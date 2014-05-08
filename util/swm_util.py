# -*- coding: utf-8 -*-
"""Utilities to depoly, configure, build, run and postprocess an experiment
of the shallow water model on the HPC Cluster of Kiel University. However,
it should work on all clusters with an PBS like batch system and Python 2.7 or
higher available on the nodes.
"""

__author__ = "Martin Claus <mclaus@geomar.de>"

import os
import shutil
import subprocess
import tempfile
from string import Template
import sys
if sys.version_info[:2] < (2, 7):
    DictClass = dict
else:
    # If Python 2.7 or higher use OrderedDict to preserve
    # the order of the Namelist
    from collections import OrderedDict as DictClass


MODULE_NAME = "swm_util"
MODEL_NL = "model_nl"
CALENDAR_NL = "calendar_nl"
DOMAIN_NL = "domain_nl"
SWM_BS_NL = "swm_bs_nl"
SWM_FORCING_NL = "swm_forcing_nl"
DIAG_NL = "diag_nl"
OUTPUT_NL = "output_nl"

NML_LINE_LENGTH = 70

class Namelist(DictClass):
    """ Class to handle Fortran Namelists
    Namelist(string) -> new namelist with fortran nml identifier string
    Namelist(string, init_val) -> new initialized namelist with nml identifier
        string and init_val beeing a valid initialisation object for the parent
        class (either OrderedDict for Python >= 2.7 or else dict).
    A fortran readable string representation of the namelist can be generated
    via str() build-in function. A string representation of the Python object
    that can be used with eval or string.Template substitution can be obtained
    by repr() build-in function. 
    """
    @property
    def name(self):
        """ Read only property name, representing the fortran namelist
        identifier.
        """
        return self._name
    
    def __init__(self, name, init_val=()):
        """x.__init__(...) initializes x; see help(type(x)) for signature"""
        self._name = name
        super(self.__class__, self).__init__(init_val)

    def __str__(self):
        """x.__str__(self) -> Fortran readable string representation of the
        namelist. If a value v is a sequence, an 1D fortran array representation
        is created using iter(v).
        """
        retstr = "&%s\n" % str(self.name)
        for k, v in self.items():
            if hasattr(v, '__iter__'):
                retstr += "%s = (/ " % k
                tmpstr = ""
                for vv in v:
                    tmpstr += "%s," % repr(vv)
                    if len(tmpstr) > NML_LINE_LENGTH:
                        if vv == v[-1]:
                            tmpstr = tmpstr[:-1]
                        retstr += tmpstr + " &\n"
                        tmpstr = ""
                retstr = retstr + tmpstr[:-1] + "/)\n"
            else:
                retstr += "%s = %s\n" % (str(k), repr(v))
        retstr += "&end\n"
        return retstr

    def __repr__(self):
        """x.__repr__(self) -> string that can be used by eval to create a copy
        of x.
        """
        retstr = "%s.%s(%s, (" % (MODULE_NAME, self.__class__.__name__,
                                 repr(self.name))
        for k, v in self.items():
            retstr += "%s, " % repr((k, v))
        retstr += "))"
        return retstr
    
    def has_name(self, name):
        """x.hasname(self, name) <==> name==x.name"""
        return name == self.name
        

class ModelController(object):
    """ Class to manage experiments on HPC Cluster
    ModelController() -> new ModelControler object
    """

    _default_namelists = [    
        Namelist(MODEL_NL, ( 
                      ("r", 0.),
                      ("k", 0.),
                      ("g", 9.80665),
                      ("gamma_new", 0.),
                      ("new_sponge_efolding", 8e-1),
                      ("gamma_new_sponge", 1.),
                      ("ah", 50.),
                      ("run_length", 100.),
                      ("dt", 1.),
                      ("meant_out", 2.628e5),
                      ("diag_start", 0.))),
        Namelist(CALENDAR_NL, (("model_start", "2000-01-01 00:00"),)),
        Namelist(DOMAIN_NL, (
                    ("a", 6371000.),
                    ("omega", 7.272205e-05),
                    ("rho0", 1024),
                    ("nx", 551),
                    ("ny", 201),
                    ("lon_s", -45.),
                    ("lon_e", 10.),
                    ("lat_s", -10.),
                    ("lat_e", 10.),
                    ("h_overwrite", 1.))),
        Namelist(SWM_BS_NL, (
                    ("bs_filename", ""),
                    ("bs_varname", "PSI"),
                    ("bs_chunksize", 0))),
        Namelist(SWM_FORCING_NL, (
                    ("filename", "tau_x.nc"),
                    ("varname", "TAUX"),
                    ("forcingtype", "W"),
                    ("component", "Z"))),
        Namelist(DIAG_NL, (
                    ("filename", "u_out"),
                    ("varname", "U"),
                    ("type", "S"),
                    ("frequency", 100))),
        Namelist(DIAG_NL, (
                    ("filename", "v_out"),
                    ("varname", "V"),
                    ("type", "S"),
                    ("frequency", 100))),
        Namelist(DIAG_NL, (
                    ("filename", "eta_out"),
                    ("varname", "ETA"),
                    ("type", "S"),
                    ("frequency", 100))),
        Namelist(OUTPUT_NL, (
                    ("oprefix", "default"),
                    ("osuffix", "")))]

    _default_header = {"freq_wind": 4.354613902181461e-08,
                       "tau_scale": 1.,
                       "fxdep": "* SIN(FREQ_WIND * dt * itt)"}

    _header_template_string = """
#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN

/* OpenMP */
#define OMPCHUNK Nx
#define OMPSCHEDULE GUIDED

/* Switch for Shallow Water Model */
#define SWM
#define SWM_TSTEP_HEAPS
/* Switches for forcing */
#define FREQ_WIND ${freq_wind}
#define TAU_SCALE ${tau_scale}
#define FXDEP ${fxdep}
/* Switches for damping */
#define LATERAL_MIXING

#define H_OVERWRITE

/* Possible calculation of Coriolis-Parameter */
#define NOF 0
#define FPLANE 1
#define BETAPLANE 2
#define SPHERICALGEOMETRY 3

/* Standard calculation of Coriolis-Parameter */
#ifndef CORIOLIS
#define CORIOLIS SPHERICALGEOMETRY
#endif

#endif
"""

    def __init__(self):
        """x.__init__(self) initialize the object.
            Reads environmental variables set by the batch system.
        """
        self._tempdirs = []
        self._dirs = {}
        self._vars = {}
        self._set_vars()

    def __del__(self):
        """ Calls x.cleanup(self) to delete temporary data on the cluster node."""
        self.cleanup()
        super(self.__class__, self).__del__()

    def cleanup(self):
        """x.cleanup(self) -> None. Deletes temporary data on the cluster node."""
        for i in reversed(range(len(self._tempdirs))):
            d = self._tempdirs[i]
            shutil.rmtree(d)
            del self._tempdirs[i]
        

    def deploy_experiment(self):
        """x.deploy_expreriment(self) -> None
            Creates temporary directories (hopefully on the local file system)
            and copies the model source and input data there.
        """
        # create directories in $TMPDIR
        self._dirs['lwrkdir'] = self._create_tempdir(prefix="swm")
        self._dirs['loutdir'] = self._create_localdir("output")
        self._dirs['lindir'] = self._create_localdir("input", path_only=True)
        self._dirs['lmodeldir'] = self._create_localdir("swm", path_only=True)

        # copy model an input files
        shutil.copytree(self._dirs['modeldir'], self._dirs['lmodeldir'])
        shutil.copytree(self._dirs['indir'], self._dirs['lindir'])

    def write_header(self, header_dic=None):
        """x.write_header(self, header_dic) -> string
            Writes header constructed from self._default_header and
            header_dic to the local include directory and returns the string.
        """
        # create header file content
        header_string = self._formatted_header(header_dic)
        # write header to local model dir
        with open(os.path.join(self._dirs['lmodeldir'],
                               'include', 'model.h'), 'w') as h_file:
            h_file.write(header_string)
        h_file.close()
        return header_string

    def compile_model(self):
        """x.compile_model(self) -> None. Compiles the deployed model."""
        pwd = os.getcwd()
        # compiling model
        os.chdir(self._dirs['lmodeldir'])
        subprocess.call('make model >/dev/null', shell=True)
        os.chdir(pwd)

    def write_namelist(self, nml=None):
        """x.write_namelist(self, nml) -> string
            Writes a fortran readable namelist file to the local work directory
            and returns its string representation.
        """
        # Setting parameters
        if not nml:
            nml = []
        nml_string = self._formatted_namelist(nml)
        # write namelist to local work dir
        with open(os.path.join(self._dirs['lwrkdir'],
                               'model.namelist'), 'w') as nl_file:
            nl_file.write(nml_string)
        nl_file.close()
        return nml_string

    def run_model(self):
        """x.run_model(self) -> None. Start the model run."""
        pwd = os.getcwd()
        os.chdir(self._dirs['lwrkdir'])
        subprocess.call('time ' + os.path.join(self._dirs['lmodeldir'],
                                               'build', 'model'),
                        shell=True)
        os.chdir(pwd)

    def post_process(self):
        """x.post_process(self) -> None.
            Moves the model output to the global file system.
        """
        # post process
        dir_list = os.listdir(self._dirs['loutdir'])
        for file in dir_list:
            print "Move: {0} -> {1}".format(os.path.join(
                                                self._dirs['loutdir'], file),
                                            self._dirs['outdir'])
            shutil.move(os.path.join(self._dirs['loutdir'], file),
                        os.path.join(self._dirs['outdir'], file))

    def _set_vars(self):
        # Setting up variables
        self._vars['jobid'] = str(os.getenv('PBS_JOBID')).split('.')[0].split(':')[-1]
        self._vars['pbs_jobname'] = str(os.getenv('PBS_JOBNAME'))
        self._dirs['workdir'] = str(os.getenv('PBS_O_WORKDIR'))
        self._dirs['modeldir'] = os.path.join(self._dirs['workdir'], "swm")
        self._dirs['outdir'] = os.path.join(self._dirs['workdir'], 'output')
        self._dirs['indir'] = os.path.join(self._dirs['workdir'], 'input')
        self._dirs.update({'lwrkdir': '', 'loutdir':'', 'lindir':'',
                           'lmodeldir':''})

    def _formatted_namelist(self, in_nml_list):
        settings = self._default_namelists
        nml_name_list = []
        nml_list = []
        for nml in settings:
            if nml.name in [SWM_FORCING_NL, DIAG_NL]:
                if nml.name in nml_name_list:
                    continue
                if nml.name in [nl.name for nl in in_nml_list]:
                    nml_list.extend(
                        [nl for nl in in_nml_list if nl.name == nml.name])
                else:
                    nml_list.append(nml)
                nml_name_list.append(nml.name)
            else:
                nml_list.append(nml)
                nml_name_list.append(nml.name)
                in_nml = [nl for nl in in_nml_list if nl.name == nml.name]
                if len(in_nml) == 1:
                    nml_list[-1].update(in_nml[0])
                elif len(in_nml) > 1:
                    raise ValueError("More than one namelist %s defined" % nml.name)

        nml_string = ""
        for nml in nml_list:
            if nml.name == OUTPUT_NL:
                nml["oprefix"] = os.path.join(self._dirs['loutdir'], nml['oprefix'])
            elif nml.name == SWM_FORCING_NL and 'filename' in nml:
                nml["filename"] = os.path.join(self._dirs['lindir'], nml['filename'])
            elif nml.name == SWM_BS_NL and "bs_filename" in nml:
                nml['bs_filename'] = os.path.join(self._dirs['lindir'],
                                               nml['bs_filename'])
            nml_string += str(nml)
        self._vars['formatted_namelist'] = nml_string
        return nml_string

    def _get_nml_from_list(self, nml_name, nml_list):
        out_list = []
        for nml in nml_list:
            if nml.getname() == nml_name:
                out_list.append(nml)
        return out_list

    def _formatted_header(self, dic):
        settings = self._default_header.copy()
        if dic:
            dict_lower = dict((k.lower(), v) for k, v in dic.iteritems())
            settings.update(dict_lower)
        header_template = Template(self._header_template_string)
        self._vars['formatted_header'] = header_template.substitute(settings)
        return self._vars['formatted_header']

    def _create_tempdir(self, prefix="tmp"):
        path = tempfile.mkdtemp(prefix=prefix)
        self._tempdirs.append(path)
        return path

    def _create_localdir(self, d, path_only=False):
        path = os.path.join(self._dirs['lwrkdir'], d)
        if not path_only:
            os.mkdir(path, 0700)
        return path
