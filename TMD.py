import numpy as np
import scipy as sp
import math
import lhapdf
from scipy import integrate
from scipy.special import gamma
from scipy.special import jv
#............................constant..........................

gammaE=np.euler_gamma
CF=4.0/3.0
p = lhapdf.mkPDF("CT14lo",0)

def a_s(var):
    mu=var
    a_s=p.alphasQ(mu)
    return a_s

def f_1(var):
    x=var[0]
    mu=var[1]
    f_1=p.xfxQ(21, x, mu)/x
    return f_1

def d_1(var):
    z=var[0]
    mu=var[1]
    d_1=z+mu
    return d_1

#...........................TMD evolution.........................

def evo_fac(var):
    mu=var[0]
    b=var[1]

    evofac=(abs( mu**2/zeta(var)) )**(-2.0*rad(var))
    return evofac

# rapidity anomalous dimension
def rad(var):
    mu=var[0]
    b=var[1]
    rad=rad_pert(var)+d_np(b)
    return rad

def rad_pert(var):
    mu=var[0]
    b=var[1]
    radpert=a_s(mu)*l_mu(var)*2*CF
    return radpert

def b_star(var):
    BNPf=2.29
    b=var
    bstar=b/math.sqrt(1+b**2/BNPf)
    return bstar

def d_np(var):
    c0=0.022
    b=var
    dnp=c0*b*b_star(b)
    return dnp

def l_mu(var):
    mu=var[0]
    b=var[1]
    lmu=math.log(mu**2*b**2/4/math.exp(-2*gammaE))
    return lmu

# zeta
def zeta(var):
    mu=var[0]
    b=var[1]

    tb=2.0*math.exp(-gammaE)/b
    ob=4.0*math.exp(-2.0*gammaE+3.0/2.0)/b**2
    
    zetaa=(ob*rad([tb,b])+integrate.quad(plus,tb,mu)[0])/rad(var)
    return zetaa

def plus(var):
    mu=var
    pluss=(cad(mu)*mu-gamma_v(mu))/2
    return pluss

def cad(var):
    mu=var
    cadd=a_s(mu)*4*CF
    return cadd

def gamma_v(var):
    mu=var
    gammav=a_s(mu)*(-6)*CF
    return gammav

#..........................TMD distribution...............
def f1(var):
    x=var[0]
    mu=var[1]
    b=var[2]
    f1=f_1([x,mu])*f_1NP([x,b])
    return f1

def d1(var):
    z=var[0]
    mu=var[1]
    b=var[2]
    d1=d_1([z,mu])*d_1NP([z,b])
    return d1

def f_1NP(var):
    lam=0.86

    x=var[0]
    b=var[1]
    f_1NP=lam*g1(x)**2*b**2/4/(1+lam*g1(x))
    f_1NP=(1-f_1NP)/2/math.pi
    f_1NP=f_1NP*math.exp(-g1(x)*b**2/4)
    return f_1NP

def d_1NP(var):
    lamF=5.5

    z=var[0]
    b=var[1]
    d_1NP=(lamF/z**2)*g4(z)**2*(1-g4(z)*b**2/4/z**2)*math.exp(-g4(z)**2*b**2/4/z**2)
    d_1NP=d_1NP+g3(z)*math.exp(-g3(z)*b**2/4/z**2)
    d_1NP=d_1NP/(2*math.pi*z**2*(g3(z)+(lamF/z**2)*g4(z)**2))
    return d_1NP

def g1(var):
    N1=0.28
    alpha=2.95
    sig=0.17

    x=var
    g1=N1*(1-x)**alpha*x**sig/0.9**alpha/0.1**sig
    return g1

def g3(var):
    N3=0.21
    beta=1.65
    delta=2.28
    gam=0.14

    z=var
    g3=N3*(z**beta+delta)*(1-z)**gam/(0.5**beta+delta)/0.5**gam
    return g3

def g4(var):
    N4=0.13
    beta=1.65
    delta=2.28
    gam=0.14

    z=var
    g4=N4*(z**beta+delta)*(1-z)**gam/(0.5**beta+delta)/0.5**gam
    return g4

#..............................sivers fit..................................

def fsiver(var,par):
    x=var[0]
    b=var[1]
    beta=par[0]
    eps=par[1]
    Nq=par[2]
    r0=par[3]
    r1=par[4]
    r2=par[5]

    fsiver=math.exp(-(r0+x*r1)*b**2/math.sqrt(1+r2*x**2*b**2))
    fsiver=Nq*(1-x)*x**beta*(1+eps*x)/n([beta,eps])*fsiver
    return fsiver

def n(par):
    beta=par[0]
    eps=par[1]
    n=(3+beta+eps+beta*eps)*gamma(beta+1)/gamma(beta+4)
    return n

#..................................AUT for theory..........................
def aut_th(var,par):
    x=var[0]
    Qf=var[1]
    z=var[2]
    ph=var[3]

    beta=par[0]
    eps=par[1]
    Np=par[2]
    r0=par[3]
    r1=par[4]
    r2=par[5]

    allpar=(x,Qf,z,ph,beta,eps,Np,r0,r1,r2)
    aut_th=-integrate.quad(aut_inte1,0.1,5.0,args=allpar)[0]/integrate.quad(aut_inte2,0,5.0,args=allpar)[0]
    return aut_th

def aut_inte1(b,x,Qf,z,ph,beta,eps,Np,r0,r1,r2):
    Q=math.sqrt(Qf)
    par=[beta,eps,Np,r0,r1,r2]
    aut1=b**2*jv(1,b*ph/z)*evo_fac([Q,b])*fsiver([x,b],par)*d1([z,Q,b])
    return aut1

def aut_inte2(b,x,Qf,z,ph,beta,eps,Np,r0,r1,r2):
    Q=math.sqrt(Qf)
    par=[beta, eps, Np, r0, r1, r2]
    aut2=b*jv(0,b*ph/z)*evo_fac([Q,b])*f1([x,Q,b])*d1([z,Q,b])
    return aut2


'''
var=[0.10,0.1,0.1,0.1]
par=[1,1,1,1,1,1]

x=var[0]
Qf=var[1]
z=var[2]
ph=var[3]

beta=par[0]
eps=par[1]
Np=par[2]
r0=par[3]
r1=par[4]
r2=par[5]

#print(4.0**(-0.2))
for i in range(1000):
	b=0.001+1.0*i
	print(zeta([Qf,b]))
'''

