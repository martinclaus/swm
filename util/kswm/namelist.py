# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 10:27:57 2014

@author: mclaus
"""
import sys
import re

if sys.version_info[:2] < (2, 7):
    DictClass = dict
else:
    # If Python 2.7 or higher use OrderedDict to preserve
    # the order of the Namelist
    from collections import OrderedDict as DictClass

MODULE_NAME = "namelist"
NML_LINE_LENGTH = 70
# Config file parser, called from the class initialization
varname   = r'[a-zA-Z][a-zA-Z0-9_]*'
valueBool = re.compile(r"(\.(true|false|t|f)\.)",re.I)
quote = re.compile(r"([\'\"]{1}.*[\'\"]{1})")
namelistname = re.compile(r"&(" + varname + r")")
paramname = re.compile(r"^(" + varname + r")")
namlistend = re.compile(r"^\$(end)?", re.I)
comment = re.compile(r"#.*")
equalsign = re.compile(r"^=$")
computation = re.compile(r"^([0-9\.e]+\s*[\*\+\-/]{1}\s*)+[0-9\.e]+", re.I)



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
                    if isinstance(vv, bool):
                        if vv:
                            rvv = ".TRUE."
                        else:
                            rvv = ".FALSE."
                    else:
                        rvv = repr(vv)
                    tmpstr += "%s," % rvv
                    if len(tmpstr) > NML_LINE_LENGTH:
                        if vv == v[-1]:
                            tmpstr = tmpstr[:-1]
                        retstr += tmpstr + " &\n"
                        tmpstr = ""
                retstr = retstr + tmpstr[:-1] + "/)\n"
            else:
                if isinstance(v, bool):
                    if v:
                        rv = ".TRUE."
                    else:
                        rv = ".FALSE."
                else:
                    rv = repr(v)
                retstr += "%s = %s\n" % (str(k), rv)
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

def parse_namelist_file(in_file):
    """ parse_namelist_file(fobj) -> list of Namelist instances. fobj can be
    any object that implements pythons file object API, i.e. that offers a
    read() method.
    """
    retlist = []
    content = _tokenize(in_file.read())
    in_file.seek(0, 0)
    for item in content:
        match = re.match(namelistname, item)
        if match:
            nmlname = match.group(1)
            nml = Namelist(nmlname)
            retlist.append(nml)
            continue
        match = re.match(paramname, item)
        if match:
            pname = match.group(1)
            nml[pname] = []
            continue
        for pattern in (namlistend, equalsign):        
            match = re.match(pattern, item)
            if match:
                continue
        match = re.match(valueBool, item)
        if match:
            nml[pname].append(match.group(1)[1].lower()=="t")
            continue
        match = re.match(quote, item)
        if match:
            nml[pname].append(match.group(1)[1:-1])
            continue
        try:
            nml[pname].append(int(item))
        except ValueError:
            pass
        else:
            continue
        try:
            nml[pname].append(float(item))
        except ValueError:
            pass
        else:
            continue
        match = re.match(computation, item)
        if match:
            nml[pname].append(eval(item))
    for nml in retlist:
        for k, v in nml.iteritems():
            if len(v) == 1:
                nml[k] = v[0]
    return retlist

def _tokenize(text):
    fs = "$FS$"
    text = re.sub(comment, '', text)        
    for char, rep in zip((r'\n', r',', ' ', '=', ), (fs, fs, fs, fs+'='+fs)):
        text = text.replace(char, rep)
    text = text.split(fs)
    return [token.strip() for token in text if token.strip() != '']