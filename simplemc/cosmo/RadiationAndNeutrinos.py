##
# This is a class that provides relevant support
# for treating radiation and neutrinos.
# Much of this functionality was in BasicCosmo,
# but it became clutterish there.
##


from simplemc.cosmo.NuDensity import ZeroNuDensity, NuDensity
from simplemc.cosmo import cosmoApprox as CA
from simplemc.cosmo.paramDefs import mnu_par, Nnu_par
import sys


class RadiationAndNeutrinos:
    """
    Class for neutrinos, Mnu and Nnu

    Parameters
    ----------
    mnu : float
        Standard value for the mass of neutrinos.
    Nnu : float
        Standard value for the number of families.
    varyMnu : bool
        Whether varying mass for neutrinos.
    varyNnu : bool
        Whether varying families for neutrinos.
    degenerate : bool
        Degenerate neutrinos.
    disable : bool
        Set radiation to zero.
    """
    # Type this into google
    # 8*pi^5*(boltzmann constant)^4/(15*(h*c)^3))*(1 Kelvin)**4/(3*(100 km/s/Mpc)^2/(8*Pi*G)*(speed of light)^2)
    # to get
    omrad_fac = 4.48130979e-7
    omrad_pref_ = omrad_fac*CA.Tcmb**4

    # check  4.48130979e-7*(2.725**4)=2.47099e-05 c.f. Doddelson p. 41
    # now Omrad= omread_pref_/h**2

    def __init__(self, mnu=mnu_par.value, Nnu=Nnu_par.value,
                 varyMnu=False, varyNnu=False, degenerate=False, disable=False):
        
        self.disabled = disable
        if (self.disabled):
            self.Omrad = 0
            self.Omnuh2 = 0
            self.Omnu = 0
            self.NuDensity = ZeroNuDensity()
            return

        self.varyMnu = False
        self.varyNnu = False
        self.NuDensity = NuDensity(CA.Tcmb, Nnu, mnu, degenerate, fact=self.omrad_fac)
        # print "Relic neutrino density:",self.NuDensity.omnuh2today


    def setVaryMnu(self, T=True):
        if self.disabled:
            print("Cannot vary radiation parameter if disabled")
            sys.exit(1)
        self.varyMnu = T


    def setVaryNnu(self, T=True):
        if self.disabled:
            print("Cannot vary radiation parameter if disabled")
            sys.exit(1)
        self.varyNnu = T


    def setMnu(self, mnu):
        if self.disabled:
            print("Cannot vary radiation parameter if disabled")
            sys.exit(1)
        self.NuDensity.setMnu(mnu)


    def setNnu(self, mnu):
        if self.disabled:
            print("Cannot vary radiation parameter if disabled")
            sys.exit(1)
        self.NuDensity.setNnu(mnu)


    def mnu(self):
        if self.disabled:
            return 0
        return self.NuDensity.mnu_


    def Nnu(self):
        if self.disabled:
            return 0
        return self.NuDensity.Nnu_


    def freeParameters(self):
        if self.disabled:
            return []
        l = []
        if (self.varyMnu):
            mnu_par.setValue(self.mnu())
            l.append(mnu_par)
        if (self.varyNnu):
            Nnu_par.setValue(self.Nnu())
            l.append(Nnu_par)
        return l


    def updateParams(self, pars):
        if self.disabled:
            return True
        for p in pars:
            if p.name == "mnu":
                self.setMnu(p.value)
            elif p.name == "Nnu":
                self.setNnu(p.value)
        # This assumes we know about h.
        # This base class alone doesn't, but if inherited in say
        # LCDM, it should. If it doesn't, the best thing is to
        # crash anyway, rather then to try to do something clever
        self.Omrad = self.omrad_pref_/(self.h**2)
        self.Omnuh2 = self.NuDensity.omnuh2today
        self.Omnu = self.Omnuh2/self.h**2
        return True
