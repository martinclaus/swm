# -*- coding: utf-8 -*-
"""Utilities to depoly, configure, build, run and postprocess an experiment
of the shallow water model on the HPC Cluster of Kiel University. However,
it should work on all clusters with an PBS like batch system and Python 2.7 or
higher available on the nodes.
"""

__author__ = "Martin Claus <mclaus@geomar.de>"

import os
import shutil
import subprocess32 as subprocess
import tempfile
from string import Template

from .namelist import Namelist

MODEL_NL = "model_nl"
CALENDAR_NL = "calendar_nl"
DOMAIN_NL = "domain_nl"
SWM_TIMESTEP_NL = "swm_timestep_nl"
SWM_FORCING_NL = "swm_forcing_nl"
DIAG_NL = "diag_nl"
OUTPUT_NL = "output_nl"

NML_LINE_LENGTH = 70

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
                    ("in_file_H", ""),
                    ("h_overwrite", 1.))),
        Namelist(SWM_TIMESTEP_NL, (
                    ("filename", ""),
                    ("varname", "PSI"),
                    ("chunksize", 0),
                    ("ab_chi", .1))),
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

    _default_header = {"tau_scale": 1.,
                       "fxdep": "* 1._8 ",
                       "velocity_sponge": "",
                       "newtonian_sponge": "",
                       "model": "FULLY_NONLINEAR"}

    _header_template_string = """
#ifndef FILE_MODEL_SEEN
#define FILE_MODEL_SEEN

/* OpenMP */
#define OMPCHUNK Nx
#define OMPSCHEDULE GUIDED

/* Switch for Shallow Water Model */
#define SWM
#define BAROTROPIC
#define ${model}
#define SWM_TSTEP_ADAMSBASHFORTH

/* Switches for forcing */
#define TAU_SCALE ${tau_scale}
#define FXDEP ${fxdep}

/* Switches for damping */
#define LATERAL_MIXING
#define NEWTONIAN_COOLING
#define LINEAR_BOTTOM_FRICTION
#define VELOCITY_SPONGE "${velocity_sponge}"
#define NEWTONIAN_SPONGE "${newtonian_sponge}"

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

    def __init__(self, workdir=None, model="swm", debug=False):
        """x.__init__(self, model="swm", debug=False) initialize the object.
            Reads environmental variables set by the batch system.
            model : directory containing the model
            debug : if true, verbose debugging information are printed to stdout
        """
        self.model = model
        self.debug = debug
        self._tempdirs = []
        self._dirs = {}
        if workdir:
            self._dirs["workdir"] = workdir
        self._vars = {}
        self._set_vars()

    def __del__(self):
        """ Calls x.cleanup(self) to delete temporary data on the cluster node."""
        self.cleanup()

    def cleanup(self):
        """x.cleanup(self) -> None. Deletes temporary data on the cluster node."""
        for i in reversed(range(len(self._tempdirs))):
            d = self._tempdirs[i]
            if self.debug: print "Remove %s" % self._tempdirs[i]
            shutil.rmtree(d)
            del self._tempdirs[i]
        

    def deploy_experiment(self):
        """x.deploy_expreriment(self) -> None
            Creates temporary directories (hopefully on the local file system)
            and copies the model source and input data there.
        """
        # create directories in $TMPDIR
        self._dirs['lwrkdir'] = self._create_tempdir(prefix=self.model)
        self._dirs['loutdir'] = self._create_localdir("output")
        self._dirs['lindir'] = self._create_localdir("input", path_only=True)
        self._dirs['lmodeldir'] = self._create_localdir(self.model,
                                                        path_only=True)

        # copy model an input files
        for dir_key in ('modeldir', 'indir'):
            shutil.copytree(self._dirs[dir_key], self._dirs['l' + dir_key])
            if self.debug:
                print "Move %s -> %s" % (self._dirs[dir_key], self._dirs['l' + dir_key])

    def write_header(self, header_dic=None):
        """x.write_header(self, header_dic) -> string
            Writes header constructed from self._default_header and
            header_dic to the local include directory and returns the string.
        """
        # create header file content
        header_string = self._formatted_header(header_dic)
        # write header to local model dir
        h_path = os.path.join(self._dirs['lmodeldir'], 'include', 'model.h')
        with open(h_path, 'w') as h_file:
            h_file.write(header_string)
        h_file.close()
        if self.debug: print "%s \n%s" % (h_path, header_string)
        return header_string

    def compile_model(self, target="model"):
        """x.compile_model(self, target="model") -> None. Compiles the deployed model with target as make target."""
        pwd = os.getcwd()
        # compiling model
        os.chdir(self._dirs['lmodeldir'])
        make_process = subprocess.Popen('make %s 2>&1' % target, shell=True,
                                        stdout=subprocess.PIPE,
                                        executable='/bin/bash')
        output = make_process.communicate()[0]
        if self.debug:
            print "Compile model:\n%s" % output
        if make_process.returncode != 0:
            raise RuntimeError("make failed with exit status %d" % (-1 * make_process.returncode))
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
        nml_path = os.path.join(self._dirs['lwrkdir'], 'model.namelist')
        with open(nml_path, 'w') as nl_file:
            nl_file.write(nml_string)
        nl_file.close()
        if self.debug: print "%s\n%s" % (nml_path, nml_string)
        return nml_string

    def run_model(self):
        """x.run_model(self) -> None. Start the model run."""
        pwd = os.getcwd()
        os.chdir(self._dirs['lwrkdir'])
        run_process = subprocess.Popen('time ' + os.path.join(self._dirs['lmodeldir'],
                                               'bin', 'model') + ' 2>&1',
                        shell=True, 
                        stdout=subprocess.PIPE,
                        executable='/bin/bash')
        output = run_process.communicate()[0]
        if self.debug: print "Run model:\n%s" % output
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
        # move profiling output to global FS
        profiling_output = os.path.join(self._dirs['lwrkdir'], "gmon.out")
        if os.path.exists(profiling_output):
            print "Move Profiling output ..."
            shutil.move(profiling_output,
                        os.path.join(self._dirs["workdir"], "gmon.out"))
            shutil.move(os.path.join(self._dirs['lmodeldir'], 'bin', 'model'),
                        os.path.join(self._dirs["workdir"], "model"))

    def _set_vars(self):
        # Setting up variables
        self._vars['jobid'] = str(os.getenv('PBS_JOBID')).split('.')[0].split(':')[-1]
        self._vars['pbs_jobname'] = str(os.getenv('PBS_JOBNAME'))
        if 'workdir' not in self._dirs:
            self._dirs['workdir'] = str(os.getenv('PBS_O_WORKDIR'))
        self._dirs['modeldir'] = os.path.join(self._dirs['workdir'], self.model)
        self._dirs['outdir'] = os.path.join(self._dirs['workdir'], 'output')
        self._dirs['indir'] = os.path.join(self._dirs['workdir'], 'input')
        if self.debug:
            for dic in (self._vars, self._dirs):
                print(str(dic))
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
        filename_map = {DOMAIN_NL: ["in_file_H"],
                        SWM_FORCING_NL: ["filename"],
                        SWM_TIMESTEP_NL: ["filename"]}
        for nml in nml_list:
            for key in nml.keys():
                if nml.name == OUTPUT_NL and key == "oprefix":
                    nml[key] = os.path.join(self._dirs['loutdir'], nml[key])
                elif nml.name in filename_map and key in filename_map[nml.name] and nml[key]:
                    nml[key] = os.path.join(self._dirs["lindir"], nml[key])
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
        if self.debug: print "Created directory %s" % path
        self._tempdirs.append(path)
        return path

    def _create_localdir(self, d, path_only=False):
        path = os.path.join(self._dirs['lwrkdir'], d)
        if not path_only:
            os.mkdir(path, 0700)
            if self.debug: print "Created directory %s" % path
        return path
