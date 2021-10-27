import numpy as np
import scipy as sp
from scipy import integrate
from scipy.special import zeta
import matplotlib.pyplot as plt
import math
import lha

gammaE=np.euler_gamma
CF=4.0/3.0
CA=3.0
Nf=4.0
Tr=1.0/2.0

pi=np.pi

c0=0.0391
BNP=1.93

bmax=7.

beta0=11./3.*CA-4./3.*Tr*Nf
beta1=34./3.*CA**2.-20./3.*CA*Tr*Nf-4.*CF*Tr*Nf
beta2=2857./54.*CA**3.+(2.*CF**2.-205./9.*CF*CA-1415./27.*CA**2.)*Tr*Nf+(44./9.*CF+158./27.*CA)*Tr**2.*Nf**2.

d_20=CF*(CA*(404.0/27.0-14.0*zeta(3))-112.0/27.0*Tr*Nf)

Gam0=4.*CF
Gam1=4.*CF*((67./9.-pi**2./3)*CA-20./9.*Tr*Nf)
Gam2=4.*CF*(CA**2*(245.0/6.0-134.0*pi**2.0/27.0+11.0*pi**4.0/45.0+22.0/3.0*zeta(3))+\
        CA*Tr*Nf*(-418.0/27.0+40.0*pi**2.0/27.0-56.0/3.0*zeta(3))+CF*Tr*Nf*(-55.0/3.0+16.0*zeta(3))\
        -16.0/27.0*Tr**2.0*Nf**2)
        
gam0=-6.*CF
gam1=CF**2.*(-3.+4.*pi**2.-48*zeta(3))+CF*CA*(-961./27.-11*pi**2/3.+52.*zeta(3))+CF*Tr*Nf*(260./27.+4.*pi**2./3.)
gam2=CF**3.0*(-29.0-6.0*pi**2.0-16.0*pi**4.0/5.0-136.0*zeta(3)+32.0*pi**2.0/3.0*zeta(3)+480.0*zeta(5))\
        +CF**2.0*CA*(-151.0/2.0+410.0*pi**2.0/9.0+494.0*pi**4.0/135.0-1688.0/3.0*zeta(3)-16*pi**2.0/3.0*zeta(3)-240.0*zeta(5))\
        +CF*CA**2.0*(-139345.0/1458.0-7163*pi**2.0/243.0-83.0*pi**4.0/45.0+7052.0/9.0*zeta(3)-88.0*pi**2.0/9.0*zeta(3)-272.0*zeta(5))\
        +CF**2.0*Tr*Nf*(5906.0/27.0-52.0*pi**2.0/9.0-56.0*pi**4.0/27.0+1024.0/9.0*zeta(3))\
        +CF*CA*Tr*Nf*(-34636.0/729.0+5188.0*pi**2.0/243.0+44.0*pi**4.0/45.0-3856.0/27.0*zeta(3))\
        +CF*Tr**2.0*Nf**2.0*(19336.0/729.0-80.0*pi**2.0/27.0-64.0/27.0*zeta(3))

def lmu(mu,b):
    lmu=np.log(mu**2.*b**2./4.)+2.*gammaE
    return lmu
  
def v(mu,b):
    v=gam0/Gam0\
    +lha.a_s(mu)*(beta0/12.*lmu(mu,b)**2.-gam0*Gam1/Gam0**2.+(gam1+d_20)/Gam0)
    return v

def rad_resum(var):
    [mu, b]=var
    X=beta0* lha.a_s(mu)* lmu(mu, b)
    rad_resum=-Gam0/2./beta0*np.log(1.-X)+\
    lha.a_s(mu)*(1./2./beta0/(1.-X))*(-beta1*Gam0/beta0*(np.log(1.-X)+X)+Gam1*X)+\
    lha.a_s(mu)**2.*(Gam0*beta1**2.*(np.log(1-X)**2.-X**2.)/4./beta0**3.+(X**2.-2.*X-2.*np.log(1.-X))*beta1*Gam1/4./beta0**2.+\
    Gam0*beta2*X**2./4./beta0**2.-Gam2/4./beta0*X*(X-2.)+d_20)/(1.-X)**2.
    return rad_resum
    
def rad(var):
    [mu, b]=var
    bstar=b/(1.+b**2./BNP**2.)**0.5
    rad=rad_resum([mu, bstar])+c0*b*bstar
    return rad

def p_mub(mu,b):
    p_mub=2.*beta0*rad([mu,b])/Gam0
    return p_mub

def g(mu,b):
    #print(mu, b)
    p=p_mub(mu, b)
    
    g=Gam0/2./beta0**2.*(np.exp(-p)-1.+p)/lha.a_s(mu)+\
    Gam0/2./beta0**2.*(beta1/beta0*(np.exp(-p)-1.+p-p**2./2.)-Gam1/Gam0*(np.exp(-p)-1.+p)+\
    beta0*gam1/gam0*p)+\
    Gam0/2./beta0**2.*lha.a_s(mu)*((Gam1**2./Gam0**2.-Gam2/Gam0)*(np.cosh(p)-1.)+(beta1*Gam1/beta0/Gam0-beta2/beta0)*(np.sinh(p)-p)+\
    (beta0*gam2/Gam0-beta0*gam1*Gam1/Gam0**2.)*(np.exp(p)-1.))
   
    #print(g)
    return g

def zeta_pert(var):
    [mu, b]=var
    #print(v(mu, b))
    zeta_pert=mu*2.*np.exp(-gammaE)/b*np.exp(-v(mu,b))
    return zeta_pert

def zeta_exact(var):
    [mu, b]=var
    zeta_exact=mu**2.*np.exp(-g(mu, b)/rad([mu, b]))
    return zeta_exact

def zeta_NP(var):
    [mu, b]=var
    zeta_NP=zeta_pert(var)*np.exp(-b**2/BNP**2.)+zeta_exact(var)*(1.-np.exp(-b**2/BNP**2.))
    return zeta_NP


def evo(var):
    [mu, b]=var
    if b>bmax:
        return 0
    evo=(mu**2./zeta_NP(var))**(-rad(var))
    return evo**2.


print(Gam0,Gam1,Gam2)
print(gam0,gam1,gam2)
print(beta0,beta1,beta2)
print(lha.a_s(2.))
print(gammaE)












