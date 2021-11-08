import numpy as np
import scipy as sp
from scipy import integrate
from scipy.special import zeta
import matplotlib.pyplot as plt
import math
import time
import lhapdf
import lha


gammaE=np.euler_gamma
CF=4.0/3.0
CA=3.0
Nf=4.0
Tr=1.0/2.0

pi=math.pi
c0=0.0391
BNP=1.93
#BNPf=4.0

zeta_2=zeta(2.0)
zeta_3=zeta(3.0)
zeta_4=zeta(4.0)
zeta_5=zeta(5.0)

beta0=11.0/3.0*CA-4.0/3.0*Nf*Tr
beta1=34.0/3.0*CA**2.0-20.0/3.0*CA*Tr*Nf-4.0*CF*Tr*Nf
beta2=2857./54.*CA**3.+(2.*CF**2.-205./9.*CF*CA-1415./27.*CA**2.)*Tr*Nf+(44./9.*CF+158./27.*CA)*Tr**2.*Nf**2.

n_cut=1#(n_cut=1,2,3)

def gam_n(n):
    gam_n=0
    if n==0:
        gam_n=1.0
    elif n==1:
        gam_n=(67.0/9.0-pi**2.0/3.0)*CA-20.0/9.0*Tr*Nf
    elif n==2:
        gam_n=CA**2*(245.0/6.0-134.0*pi**2.0/27.0+11.0*pi**4.0/45.0+22.0/3.0*zeta_3)+\
        CA*Tr*Nf*(-418.0/27.0+40.0*pi**2.0/27.0-56.0/3.0*zeta_3)+CF*Tr*Nf*(-55.0/3.0+16.0*zeta_3)\
        -16.0/27.0*Tr**2.0*Nf**2
    else:
        print('error in gam_n: n out of number')
        return 0;
    return 4.0*CF*gam_n

#Gam0=gam_n(0)/4.0/CF
#Gam1=gam_n(1)/4.0/CF
#Gam2=gam_n(2)/4.0/CF

def gamv_n(n):
    gamv_n=0
    if n==0:
        gamv_n=-6.0*CF
    elif n==1:
        gamv_n=CF**2.0*(-3.00+4*pi**2.0-48.0*zeta_3)+CF*CA*(-961.0/27.0-11.0*pi**2.0/3.0+52.0*zeta_3)\
        +CF*Tr*Nf*(260.0/27.0+4.0*pi**2.0/3.0)
    elif n==2:
        gamv_n=CF**3.0*(-29.0-6.0*pi**2.0-16.0*pi**4.0/5.0-136.0*zeta_3+32.0*pi**2.0/3.0*zeta_3+480.0*zeta_5)\
        +CF**2.0*CA*(-151.0/2.0+410.0*pi**2.0/9.0+494.0*pi**4.0/135.0-1688.0/3.0*zeta_3-16*pi**2.0/3.0*zeta_3-240.0*zeta_5)\
        +CF*CA**2.0*(-139345.0/1458.0-7163*pi**2.0/243.0-83.0*pi**4.0/45.0+7052.0/9.0*zeta_3-88.0*pi**2.0/9.0*zeta_3-272.0*zeta_5)\
        +CF**2.0*Tr*Nf*(5906.0/27.0-52.0*pi**2.0/9.0-56.0*pi**4.0/27.0+1024.0/9.0*zeta_3)\
        +CF*CA*Tr*Nf*(-34636.0/729.0+5188.0*pi**2.0/243.0+44.0*pi**4.0/45.0-3856.0/27.0*zeta_3)\
        +CF*Tr**2.0*Nf**2.0*(19336.0/729.0-80.0*pi**2.0/27.0-64.0/27.0*zeta_3)
    else:
        print('error in gamv_n: n out of number')
        return 0;
    return gamv_n

def d_nk(n,k):
    d_nk=0
    if n==1 and k==0:
        d_nk=0
    elif n==1 and k==1:
        d_nk=2.0
    elif n==2 and k==0:
        d_nk=CA*(404.0/27.0-14.0*zeta_3)-112.0/27.0*Tr*Nf
    elif n==2 and k==1:
        d_nk=2.0*Gam1
    elif n==2 and k==2:
        d_nk=beta0
    elif n==3 and k==0:
        d_nk=CA**2.0*(297029.0/1458.0-3196.0/81.0*zeta_2-6164.0/27.0*zeta_3-77.0/3.0*zeta_4+88.0/3.0*zeta_2*zeta_3+96.0*zeta_5)\
        +CA*Nf*(-31313.0/729.0+412.0/81.0*zeta_2+452.0/27.0*zeta_3-10.0/3.0*zeta_4)\
        +CF*Nf*(-1711.0/54.0+152.0/9.0*zeta_3+8.0*zeta_4)+Nf**2.0*(928.0/729.0+16.0/9.0*zeta_3)
    elif n==3 and k==1:
        d_nk=2.0*beta0*(CA*(404.0/27.0-14.0*zeta_3)-112.0/27.0*Tr*Nf)+2.0*Gam2
    elif n==3 and k==2:
        d_nk=2.0*Gam1*beta0+beta1
    elif n==3 and k==3:
        d_nk=2.0/3.0*beta0**2.0
    else:
        print('error in d_nk')
    return CF*d_nk

def l_mu(mu,b):
    l_mu=math.log(mu**2.0*b**2.0/4.0/math.exp(-2.0*gammaE))
    return l_mu

def r3(b):
    r3=-2.*d_nk(2,0)**2./gam_n(0)**2.
    return r3
    
def p_mub(mu,b):
    p_mub=2.*beta0*rad([mu,b])/gam_n(0)
    return p_mub
    
def g_0(mu,b):

    g_0=(math.exp(-p_mub(mu,b))+p_mub(mu,b)-1.)/beta0/p_mub(mu,b)
    return g_0

def g(mu,b):
    p=p_mub(mu, b)
    g0=g_0(mu, b)
    a_s=lha.a_s(mu)
    
    g=g0/a_s+g0*(beta1/beta0-gam_n(1)/gam_n(0))+gamv_n(1)/gamv_n(0)-beta1*p/2./beta0**2.+\
    a_s*(g0*(gam_n(0)*beta2-beta1*gam_n(1))/beta0/gam_n(0)+(np.cosh(p)-1.)/p*(beta0*gam_n(1)**2.-beta0*gam_n(0)*gam_n(2)+\
    beta1*gam_n(0)*gam_n(1)-beta2*gam_n(0)**2.)/beta0**2./gam_n(0)**2.+(np.exp(p)-1.)/p*(gam_n(0)*gamv_n(2)-gam_n(1)*gamv_n(1))/gam_n(0)**2.)
    
    
    return g

#...................evolution........................

def evo_fac(var):#演化因子，由函数rad,zeta决定
    mu=var[0]
    b=var[1]
    #print(var)
    
    muf=mu
    ztf=muf**2.0
    bstar=b/math.sqrt(1.0+b**2.0/BNP**2.)
    mui=2.0*math.exp(-gammaE)/b+2.
    zti=4.0*math.exp(-2.0*gammaE+3.0/2.0)/b**2.0
    
    
    if  b>b_x(mu):
        return 0
    
    #path1=integrate.quad(gam_F, mui, muf, args=ztf)[0]-rad([mui, b])*math.log(ztf/zti)
    #path2=integrate.quad(gam_F, mui, muf, args=zti)[0]-rad([muf, b])*math.log(ztf/zti)
    #print(zt_NP([mu, b]))
    #print(rad(var))
    path3=(ztf/zt_NP([mu, b]) )**(-rad(var))
    #print(var)
    #print(path3)
    #path4=ad.in_R(mui, zti, muf, b)
    #path=[math.exp(path1),math.exp(path2),path3]
    evofac=path3
    #evofac=math.exp(path4)
    #print(path3)
    return evofac**2.

def b_x(mu):
    #print(lha.a_s(mu))
    b_x=2.*math.exp(-gammaE)*math.exp(1./beta0/lha.a_s(mu)/2.)/mu
    return b_x
#..................................................


def v_mub(mu,b):
    #print((gamv_n(1)))
    v_mub=gamv_n(0)/gam_n(0)+\
    lha.a_s(mu)*(beta0/12.*l_mu(mu,b)**2.-gamv_n(0)*gam_n(1)/gam_n(0)**2.+(gamv_n(1)+d_nk(2,0))/gam_n(0))#+\
    #lha.a_s(mu)**2.*(beta0**2.*l_mu(mu,b)**3./24.+(beta1+beta0*gam_n(1)/gam_n(0))*l_mu(mu,b)**2./12.+(-beta0*\
    #gamv_n(0)*gam_n(1)/gam_n(0)**2.+(8.*d_nk(2,0)+3.*gam_n(1))*beta0/3./gam_n(0))*l_mu(mu,b)/2.+\
    #gamv_n(0)*gam_n(1)**2./gam_n(0)**3.-(gam_n(1)*(d_nk(2,0)+gamv_n(1))+gamv_n(0)*gam_n(2))/gam_n(0)**2.+(d_nk(3,0)+gamv_n(2))/gam_n(0)+\
    #r3(b)/l_mu(mu,b))
    #if v_mub<-700:
    #    v_mub=-700
    return v_mub

def zt_pert(var):
    mu=var[0]
    b=var[1]
    #print(v_mub(mu, b))
    zt_pert=2.*mu*math.exp(-gammaE)*math.exp(-v_mub(mu, b))/b
    #print(zt_pert)
    return zt_pert

def zt_exact(var):
    mu=var[0]
    b=var[1]
    zt_exact=mu**2.*math.exp(-g(mu,b))
    return zt_exact

def zt_NP(var):
    [mu,b]=var
    zt_NP=zt_pert(var)*math.exp(-b**2./BNP**2.)+zt_exact(var)*(1.-math.exp(-b**2./BNP**2.))
    #print(zt_NP)
    return zt_NP
    
def cusp(var):
    mu=var
    cusp=0
    for i in range(n_cut):
        cusp=lha.a_s(mu)**(i+1.0)*gam_n(i)+cusp
    return cusp

def cusp_mu(var):
    cusp_mu=cusp(var)/var
    return cusp_mu

def gamv(var):
    mu=var
    gamv=0
    for i in range(n_cut):
        gamv=lha.a_s(mu)**(i+1.0)*gamv_n(i)+gamv
    return gamv

def rad_per_resum(var):
    [mu, b]=var
    #print(var)
    X=beta0*lha.a_s(mu)*l_mu(mu,b)
    #print(X)
    rad_per_resum=-gam_n(0)/2./beta0*math.log(1.-X)+\
    lha.a_s(mu)*(-beta1*gam_n(0)/beta0*(math.log(1.-X)+X)+gam_n(1)*X)/2./beta0/(1.-X)#+\
    lha.a_s(mu)**2.*(gam_n(0)*beta1**2.*(math.log(1-X)**2.-X**2.)/4./beta0**3.+(X**2.-2.*X-2.*math.log(1.-X))*beta1*gam_n(1)/4./beta0**2.+\
    gam_n(0)*beta2*X**2./4./beta0**2.-gam_n(2)/4./beta0*X*(X-2.)+d_nk(2,0))/(1.-X)**2.

    #print(rad_per_resum)
    return rad_per_resum

def rad_per(var):
    mu, b=var
    rad_per=0
    for i in range(n_cut):
        n=i+1
        for k in range(n+1):
            rad_per=lha.a_s(mu)**n*l_mu(mu,b)**k*d_nk(n,k)+rad_per
    return rad_per

def diff_rad(var):
    mu, b=var
    prec=0.0001
    
    muu=mu+prec/2.0
    mud=mu-prec/2.0
    diff_rad=(rad([muu, b])-rad([mud, b]))/prec
    return diff_rad
    
def del_cusp(var):
    mu, b=var
    del_cusp=cusp(mu)-mu*diff_rad([mu,b])
    return del_cusp

def a_int(mu, ztf):#for in_R
    a_int=(16.0/3.0*lha.a_s(mu)*math.log(mu**2.0/ztf)+24.0/3.0*lha.a_s(mu))/mu
    return a_int
    
def in_R(mui, zti, muf, b):
    ztf=muf**2.0
    bstar=b/math.sqrt(1.0+b**2.0/BNPf)
    in_R=integrate.quad(a_int, mui, muf, args=ztf)[0]-(8.0/3.0*lha.a_s(mui)*math.log(mui**2*b**2/4.0/math.exp(-2.0*gammaE))+\
    c0*b*b)*math.log(ztf/zti)
    return in_R

def gam_F(mu, zt):
    gam_F=(cusp(mu)*math.log(mu**2.0/zt)-gamv(mu))/mu
    return gam_F

def rad(var):
    mu=var[0]
    b=var[1]
    bstar=b/math.sqrt(1.0+b**2.0/BNP**2.)
    mu0=2.0*math.exp(-gammaE)/b+2.0
    #rad=integrate.quad(cusp_mu,mu0,mu)[0]+rad_per([mu0,b])+c0*b*b
    rad=rad_per_resum([mu,bstar])+c0*bstar*b
    return rad

def zt(var):
    mu=var[0]
    b=var[1]
    
    bstar=b/math.sqrt(1.0+b**2.0/BNPf)
    
    mub=2.0*math.exp(-gammaE)/b
    #ob=mub**2
    ob=4.0*math.exp(-2.0*gammaE+3.0/2.0)/b**2.0
    
    ztup=math.log(ob)*rad([mub, b])+integrate.quad(asint, mub, mu)[0]
    ztdown=rad([mu, b])
    zt=ztup/ztdown
    if zt>706:
        zt=math.inf
    elif zt<-700:
        zt=math.exp(-700)
    else:
        zt=math.exp(zt)
    return zt

def asint(var):#zeta依赖函数，用于积分运算
    asint=(2.0*cusp(var)*math.log(var)-gamv(var))/var 
    return asint






