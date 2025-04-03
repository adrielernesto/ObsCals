# cython : language_level = 3
# cython : nonecheck = True
# distutils: language = c++

from lib cimport *

cdef class TovObsCalc:
    units = { "MeV4_MeV^3" : 1, "MeVfm^3_fm^3" : 2}
    cdef double rhomin, rmax, drmax, drmin
    def __init__(self,rhomin,rmax,drmin,drmax):
        self.rhomin = rhomin
        self.rmax = rmax
        self.drmin = drmin
        self.drmax = drmax
    def Calc(self,namein, nameout, unit):
        tov_main(bytes(namein,"utf8"),bytes(nameout,"utf8"),
            self.rhomin,self.drmin,self.drmax,self.rmax,unit)


cdef class GammaObsCalc:
    units = { "MeV4_MeV^3" : 1, "MeVfm^3_fm^3" : 2}
    cdef double rhomin, rmax, drmax, drmin
    def __init__(self,rhomin,rmax,drmin,drmax):
        self.rhomin = rhomin
        self.rmax = rmax
        self.drmin = drmin
        self.drmax = drmax
    def Calc(self,namein, nameout, unit):
            gamma_main(bytes(namein,"utf8"),bytes(nameout,"utf8"),
                self.rhomin,self.drmin,self.drmax,self.rmax,unit)


cdef class CilindricObsCalc:
    units = { "MeV4_MeV^3" : 1, "MeVfm^3_fm^3" : 2}
    cdef double rhomin, rmax, drmax, drmin
    def __init__(self,rhomin,rmax,drmin,drmax):
        self.rhomin = rhomin
        self.rmax = rmax
        self.drmin = drmin
        self.drmax = drmax
    def Calc(self,namein, nameout, unit):
        cilindric_main(bytes(namein,"utf8"),bytes(nameout,"utf8"),
                self.rhomin,self.drmin,self.drmax,self.rmax,unit)


cdef class GeneralCilindricObsCalc:
    units = { "MeV4_MeV^3" : 1, "MeVfm^3_fm^3" : 2}
    cdef double rhomin, rmax, drmax, drmin
    def __init__(self,rhomin,rmax,drmin,drmax):
        self.rhomin = rhomin
        self.rmax = rmax
        self.drmin = drmin
        self.drmax = drmax
    def Calc(self,namein, nameout, unit):
        general_cilindric_main(bytes(namein,"utf8"),bytes(nameout,"utf8"),
                self.rhomin,self.drmin,self.drmax,self.rmax,unit)

















