import numpy as np
from scipy import optimize
from scipy import interpolate

G_kpc = 4.302e-6#kpc/SolarMass(km/s)^2
####
data_path = "simplemc/data/Blok_McGaugh_&_Rubin_(2001)/"
data = np.loadtxt(data_path+'U4115.dat')
LBS_path = "simplemc/data/LBS/"
ditc_l = {'l0': LBS_path+'XY_l0_phi0_1.dat',
            'l1': LBS_path+'XY_l1_phi0_1.dat',
            'l2': LBS_path+'XY_l2_phi0_1.dat'}

#XY0 = np.loadtxt(LBS_path+'XY_l0_phi0_1.dat')
#XY1 = np.loadtxt(LBS_path+'XY_l1_phi0_1.dat')
vecRp_data = np.array([row[1] for row in data])# galactocentric distance [kpc]
def Mass_func(r,phi,l):
    Int = np.zeros(len(r))
    dr = np.diff(r)[0]
    for i in range(0,len(r)-1):
        Int[i+1] = dr*(phi[i+1]**2.*r[i+1]**(2.*l+2.)) + Int[i]
    return Int
def Vc2_cir(r,eps,M):
    units =8.95e10*eps**2.
    return (units*M)/r
####
##Units for r in kpc
###
def r_units(r,eps,m_a):
    return (6.39e-27*r)/(eps*m_a)
def Vc_xy(r,m_a,eps,XY,l):
    M_r0 = Mass_func(XY[0],XY[1],l)#Integrates rho(r) to obtain M(r)
    Vc2_r0 = Vc2_cir(XY[0],eps,M_r0)#Vc^2[km/s]^2 theoretical
    X0_units = r_units(XY[0],eps,m_a)#r[kpc] theoretical
    M_r0_units = M_r0*eps*1.34e-10/m_a#M(r) with Solar Mass units
    if X0_units[-1]<vecRp_data[-1]:
        #array from last element of the r[kpc] theoretical to the last element of the data array,
        # with 80 elements. It can be replaced by np.arange(X0_units[-1],vecRp_data[-1],0.1) 
        #but you have to be careful in the next function with interpolate
        r_array = np.linspace(X0_units[-1],vecRp_data[-1],80)
        Vc2_rmayor = G_kpc*M_r0_units[-1]/r_array#Computes Vc^2 with with the last result from M(r)
        Vc2_total = np.append(Vc2_r0,Vc2_rmayor)#creates an array of Vc^2 with Vc2_r0 and Vc2_rmayor
        r_total = np.append(X0_units,r_array)
        return r_total,Vc2_total
    else:
        return X0_units,Vc2_r0
def Vc_inter(r,m_a,eps,l):
    if l == 0.:
        file_XY = np.loadtxt(ditc_l['l0'])
    if l == 1.:
        file_XY = np.loadtxt(ditc_l['l1'])
    if l == 2.:
        file_XY = np.loadtxt(ditc_l['l2'])
    Vc = Vc_xy(r,m_a,eps,file_XY,l)
    #If you want to use np.arange in the previous function, It is recommended to use extrapolate
    f = interpolate.interp1d(Vc[0],Vc[1],fill_value='extrapolate')
    Vc_new = f(r)
    return Vc_new