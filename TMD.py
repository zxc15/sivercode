import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
import lhapdf
import lha
from scipy import integrate
from scipy.special import gamma
from scipy.special import jv
#程序中涉及的常数
gammaE=np.euler_gamma
CF=4.0/3.0
c0=0.022
BNPf=2.29

#...........................TMD evolution.........................

def evo_fac(var):#演化因子，由函数rad,zeta决定
    mu=var[0]
    b=var[1]

    evofac=( mu**2/zeta(var) )**(-2.0*rad(var))
    return evofac
    
# rapidity anomalous dimension

def rad(var):
    mu=var[0]
    b=var[1]
    rad=lha.a_s(mu)*2*CF*math.log(mu**2*b**2/4.0/math.exp(-2.0*gammaE))+c0*b**2*(1.0+b**2/BNPf)**(-0.5)
    return rad


# zeta

def zeta(var):
    mu=var[0]
    b=var[1]

    tb=2.0*math.exp(-gammaE)/b
    ob=4.0*math.exp(-2.0*gammaE+3.0/2.0)/b**2
    
    zetaup=4.0*math.exp(-2.0*gammaE+3.0/2.0)*c0*(1.0+b**2/BNPf)**(-0.5)+CF*integrate.quad(asint,tb,mu)[0]
    zetadown=lha.a_s(mu)*2.0*CF*math.log(mu**2 *b**2/4.0/math.exp(-2.0*gammaE))+c0*b**2*(1+b**2/BNPf)**(-0.5)
    
    zeta=zetaup/zetadown
    return abs(zeta)

def asint(var):#zeta依赖函数，用于积分运算
	asint=2.0*var*lha.a_s(var)+3.0*lha.a_s(var)
	return asint




#..........................TMD distribution...............
def f1(var):#TMDpdf
    
    x=var[0]
    mu=var[1]
    b=var[2]
    target=var[3]
    
    fp=lha.f_1([x,mu,target])
    fnp=f_1NP([x,b])
    f1={i:fp[i]*fnp for i in fp}
    return f1

def d1(var):#TMDff
    
    z=var[0]
    mu=var[1]
    b=var[2]
    hardon=var[3]
    charge=var[4]
    
    dp=lha.d_1([z,mu,hardon,charge])
    dnp=d_1NP([z,b])
    d1={i:dp[i]*dnp for i in dp}
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

def fsiver(var,par):#sivers函数，其中有六个参数
    x=var[0]
    b=var[1]
    beta=par[0]
    eps=par[1]
    Nq=par[2]
    r0=par[3]
    r1=par[4]
    r2=par[5]

    fsiver=math.exp(-(r0+x*r1)*b**2/math.sqrt(1.0+r2*x**2*b**2))
    fsiver=Nq*(1-x)*x**beta*(1+eps*x)/n([beta,eps])*fsiver
    return fsiver

def n(par):
    beta=par[0]
    eps=par[1]
    n=(3+beta+eps+beta*eps)*gamma(beta+1)/gamma(beta+4)
    return n

#.............................AUT for theory..........................

def aut_th(var,par):#最终得到的AUT,由演化因子(evo_fac),TMDpdf(f1,fsiver),以及TMDff(d1)决定
    x=var[0]
    Qf=var[1]
    z=var[2]
    ph=var[3]
    target=var[4]
    hadron=var[5]
    charge=var[6]

    beta=par[0]
    eps=par[1]
    Np=par[2]
    r0=par[3]
    r1=par[4]
    r2=par[5]
    
    Q=math.sqrt(Qf)
    flavor_charge={2:2.0/3.0, 1:-1.0/3.0, 3:-1.0/3.0, 4:2.0/3.0, 5:-1.0/3.0, -2:-2.0/3.0, -1:1.0/3.0, -3:1.0/3.0, -4:-2.0/3.0, -5:1.0/3.0}
    aut_th_up=0
    aut_th_down=0
    m_target={'deuteron':1.876, 'proton':0.938, 'neutron':0.940, '3He':2.808}

    for i in m_target:
        if i==target:
            M=m_target[i]
    for i in flavor_charge:
        allpar=(x,Qf,z,ph,beta,eps,Np,r0,r1,r2,target,hadron,charge,i)
        pole=2*math.exp(-gammaE)/Q#被积函数有极点，积分在极点处分段
        inte1=integrate.quad(aut_inte1,0,pole,args=allpar)[0]+integrate.quad(aut_inte1,pole,math.inf,args=allpar)[0]
        inte2=integrate.quad(aut_inte2,0,pole,args=allpar)[0]+integrate.quad(aut_inte2,pole,math.inf,args=allpar)[0]
        
        aut_th_up=flavor_charge[i]**2*inte1+aut_th_up
        aut_th_down=flavor_charge[i]**2*inte2+aut_th_down
        
    aut_th=-M*aut_th_up/aut_th_down
    return aut_th

def aut_inte1(b,x,Qf,z,ph,beta,eps,Np,r0,r1,r2,target,hadron,charge,flavor):
    Q=math.sqrt(Qf)
    par=[beta,eps,Np,r0,r1,r2]
    ff=d1([z,Q,b,hadron,charge])
    aut1=b**2*jv(1,b*ph/z)*evo_fac([Q,b])*fsiver([x,b],par)*ff[flavor]
    return aut1



def aut_inte2(b,x,Qf,z,ph,beta,eps,Np,r0,r1,r2,target,hadron,charge,flavor):
    Q=math.sqrt(Qf)
    par=[beta, eps, Np, r0, r1, r2]
    pdf1=f1([x,Q,b,target])
    ff=d1([z,Q,b,hadron,charge])
    aut2=b*jv(0,b*ph/z)*evo_fac([Q,b])*pdf1[flavor]*ff[flavor]
    return aut2


