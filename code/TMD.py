import numpy as np
import scipy as sp
import readdate as ex
import evolution as evl
import anomalous_dimension as ad
import matplotlib.pyplot as plt
import math
import time
import lhapdf
import lha
from scipy import integrate
from scipy.special import gamma
from scipy.special import jv
#程序中涉及的常数

gammaE=np.euler_gamma
pi=math.pi
PpD=0.935
PnD=0.935
Pp=-0.028
Pn=0.86
#lam1=0.14
#lam2=0.075

lad1=0.198
lad2=9.3
lad3=431.
lad4=2.12
lad5=-4.44

eta1=0.26
eta2=0.476
eta3=0.478
eta4=0.483



flavor_use=[2, 1, 3, -2, -1, -3]

#m_target={'deuteron':1.876, 'proton':0.938, 'neutron':0.940, '3He':2.808}
m_target={'deuteron':0.939, 'proton':0.939, 'neutron':0.939, '3He':0.939}

#..........................TMD distribution...............
def mu0(b):
    mu0=2.*math.exp(-gammaE)/b+2.
    return mu0

def f1(var):#TMDpdf
    
    x=var[0]
    b=var[1]
    target=var[2]
    #print([x, b])
    #f1_high_order=f1_c3(var)+f1_c4(var)+f1_c2(var)
    #print(f1_high_order)
    fp=lha.f_1([x,mu0(b),target])
    fnp=f_np([x,b])
    #print(f1_high_order)
    f1={i:fp[i]*fnp for i in flavor_use}
    #f1={i:fp[i]*fnp+f1_high_order for i in flavor_use}
    return f1

def f1_c3(var):
    [x, b, target]=var
    f1_c3=0
    for i in flavor_use:
        allpar=(x, b, target, i)
        f1_c3=integrate.quad(inte_f1_c3, x, 1., args=allpar)[0]+f1_c3
    return f1_c3

def inte_f1_c3(y, x, b, target, flavor):
    mui=mu0(b)
    fp=lha.f_1([x/y, mui, target])[flavor]
    fnp=f_np([x,b])
    inte=2.*lha.a_s(mui)*evl.CF*(1.-y)*fp*fnp/y
    return inte

def f1_c4(var):
    [x, b, target]=var
    mui=mu0(b)
    fnp=f_np([x,b])
    f1_c4=0
    for i in flavor_use:
        fp=lha.f_1([x, mui, target])[i]
        f1_c4=-lha.a_s(mui)*evl.CF*fp*fnp*pi**2./6.+f1_c4
    return f1_c4

def f1_c2(var):
    [x, b, target]=var
    mui=mu0(b)
    fnp=f_np([x,b])
    f1_c2=0
    L=evl.lmu(mui, b)
    for i in flavor_use:
        allpar=(x, b, target, i)
        fp=lha.f_1([x, mui, target])[i]
        f1_c2=integrate.quad(inte_f1_c2, x, 1, args=allpar)[0]\
        +integrate.quad(inte_f1_c2_p1, x, 1, args=allpar)[0]\
        +integrate.quad(inte_f1_c2_p2, 0, x, args=allpar)[0]+3.*fp*fnp+f1_c2
    f1_c2=-L*lha.a_s(mui)*evl.CF*f1_c2
    return f1_c2
    
def inte_f1_c2(y, x, b, target, flavor):
    mui=mu0(b)

    fp=lha.f_1([x/y, mui, target])[flavor]
    fnp=f_np([x,b])
  
    inte=2.*(-1-y)*fp*fnp/y
    return inte
    

def inte_f1_c2_p1(y, x, b, target, flavor):
    mui=mu0(b)

    fp=lha.f_1([x/y, mui, target])[flavor]
    fnp=f_np([x,b])
    
    fp1=lha.f_1([x, mui, target])[flavor]
  
    inte=(4./y*fp*fnp-4.*fp1*fnp)/(1.-y)
    return inte

def inte_f1_c2_p2(y, x, b, target, flavor):
    mui=mu0(b)
    fnp=f_np([x,b])
    
    fp1=lha.f_1([x, mui, target])[flavor]
  
    inte=-4.*fp1*fnp/(1.-y)
    return inte

def d1(var):#TMDff
    
    z=var[0]
    b=var[1]
    hardon=var[2]
    charge=var[3]
    
    dp=lha.d_1([z,mu0(b),hardon,charge])
    dnp=d_np([z,b])
    d1={i:dp[i]*dnp/z**2. for i in flavor_use}
    return d1

def f_np(var):
    [x, b]=var
    f_np=lad1*(1.-x)+lad2*x+x*(1-x)*lad5
    f_np=-f_np*b**2./(1.+lad3*x**(lad4)*b**2.)**0.5
    f_np=math.exp(f_np)
    return f_np

def d_np(var):
    [z, b]=var
    d_np=(eta1*z+eta2*(1.-z))*b**2./(1.+eta3*(b/z)**2.)**0.5/z**2.
    d_np=math.exp(-d_np)*(1.+eta4*b**2./z**2.)
    return d_np

#print(f1_c2(  [.101, 0.11, 'proton'] ))

'''
def f_np(var):
    [z, b]=var
    f_np=math.exp(-lam2*z*b**2./(1.+z**2.*b**2.*lam2**2./lam1**2.))
    return f_np

def f_1NP(var):
    lam=0.86

    x=var[0]
    b=var[1]
    f_1NP=lam*g1(x)**2*b**2/4.0/(1.0+lam*g1(x))
    f_1NP=(1.0-f_1NP)/2.0/math.pi
    f_1NP=f_1NP*math.exp(-g1(x)*b**2.0/4.0)
    return f_1NP

def d_1NP(var):
    lamF=5.5

    z=var[0]
    b=var[1]
    d_1NP=(lamF/z**2)*g4(z)**2*(1.0-g4(z)*b**2/4.0/z**2)*math.exp(-g4(z)**2*b**2/4.0/z**2)
    d_1NP=d_1NP+g3(z)*math.exp(-g3(z)*b**2/4.0/z**2)
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
'''

def evo(var):
    [mu, b]=var
    evo=evl.evo(var)
    #evo=ad.evo_fac(var)
    return evo

#..............................sivers fit..................................
def fsiver_all(var,par):#共21个参数
    x=var[0]
    b=var[1]
    target=var[2]
    #其中参数r0,r1,r2与味道无关，并分别定义了轻味夸克的参数beta,eps,Nq
    [beta_sea, N_s, N_sea]=par[3]
    [r0, r1, r2] =par[0]
    [beta_u, eps_u, Nq_u]=par[1]
    [beta_d, eps_d, Nq_u]=par[2]
    par_s=[beta_sea, 0, N_s]
    par_sea=[beta_sea, 0, N_sea]
    
    #本例暂不考虑重味夸克的siverPDF，所以，将重味夸克的参数设置为0
    par_cb=[0,0,0,0,0,0]
    #定义质子和中子的sivers参数
    flavor_p_par={2:par[1]+par[0], 1:par[2]+par[0], 3:par_s+par[0], 4:par_cb, 5:par_cb, -2:par_sea+par[0], -1:par_sea+par[0], -3:par_sea+par[0], -4:par_cb, -5:par_cb}
    flavor_n_par={2:par[2]+par[0], 1:par[1]+par[0], 3:par_s+par[0], 4:par_cb, 5:par_cb, -2:par_sea+par[0], -1:par_sea+par[0], -3:par_sea+par[0], -4:par_cb, -5:par_cb}
    #不同靶核的siversPDF
    fsiver_p={i:fsiver([x,b],flavor_p_par[i]) for i in flavor_use}
    fsiver_n={i:fsiver([x,b],flavor_n_par[i]) for i in flavor_use}
    
    if target=='proton':   
        fsiver_all={i:fsiver([x,b],flavor_p_par[i]) for i in flavor_use}
    elif target=='neutron':
        fsiver_all={i:fsiver([x,b],flavor_n_par[i]) for i in flavor_use}
    elif target=='deuteron':
        fsiver_all={i:(PnD*fsiver_n[i]+PpD*fsiver_p[i])/2.0 for i in flavor_use}
    elif target=='3He':
        fsiver_all={i:(Pn*fsiver_n[i]+Pp*2.0*fsiver_p[i])/3.0 for i in flavor_use}
    else:
        print('error:target input error')
    return fsiver_all

def fourier_mean_std(var):
    [x, kt, mu, target, flavor]=var
    random_par=ex.data_par()
    f_fsi=[]
    for i in random_par:
        ffa=fourier_fsiver_all([x, kt, mu, target],i)[flavor]
        f_fsi.append(ffa)
    fourier_mean_std=[np.mean(f_fsi), np.std(f_fsi)]
    return fourier_mean_std

def fourier_fsiver_all(var,par):
    [x, kt, mu, target]=var
    fourier_fsiver_all={}
    
    for j in flavor_use:
        allpar=(x, kt, mu, target, j, par)
        fourier_fsiver_all[j]=integrate.quad(fourier_inte, 0., evl.bmax, args=allpar)[0]
        
    return fourier_fsiver_all
    
def fourier_f1_all(var):
    [x, kt, mu, target]=var
    fourier_f1_all={}
    
    for j in flavor_use:
        allpar=(x, kt, mu, target, j)
        fourier_f1_all[j]=integrate.quad(fourier_f1_inte, 0., evl.bmax, args=allpar)[0]
        
    return fourier_f1_all

def fourier_inte(b, x, kt, mu, target, flavor, par):
    M=0.
    for i in m_target:
        if i==target:
            M=m_target[i]
    if M==0.:
        print('error in fourier_inte')
           
    fsi=fsiver_all([x, b, target],par)[flavor]
    fourier_inte=M**2.*b**2./(2.*pi*kt)*jv(1,b*kt)*fsi*evo([mu, b])**.5
    return fourier_inte

def fourier_f1_inte(b, x, kt, mu, target, flavor):
    fsi=f1([x, b, target])[flavor]
    fourier_inte=b/(2.*pi)*jv(0,b*kt)*fsi*evo([mu, b])**.5
    return fourier_inte

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

    
    Q=math.sqrt(Qf)
    flavor_charge={2:2.0/3.0, 1:-1.0/3.0, 3:-1.0/3.0, 4:2.0/3.0, 5:-1.0/3.0, -2:-2.0/3.0, -1:1.0/3.0, -3:1.0/3.0, -4:-2.0/3.0, -5:1.0/3.0}
    
    aut_th_up=0
    aut_th_down=0
    

    for i in m_target:#判断靶例子的质量
        if i==target:
            M=m_target[i]
    for i in flavor_use:#对不同味道的积分叠加
        allpar=(x,Qf,z,ph,par,target,hadron,charge,i)
        pole=2*math.exp(-gammaE)/Q#被积函数有极点，积分在极点处分段
        inte1=integrate.quad(aut_inte1,0,pole,args=allpar)[0]+integrate.quad(aut_inte1, pole,evl.bmax ,args=allpar)[0]
        inte2=integrate.quad(aut_inte2,0,pole,args=allpar)[0]+integrate.quad(aut_inte2, pole,evl.bmax ,args=allpar)[0]
        
        aut_th_up=flavor_charge[i]**2*inte1+aut_th_up
        aut_th_down=flavor_charge[i]**2*inte2+aut_th_down
 
    aut_th=-M*aut_th_up/aut_th_down
    
    return aut_th

def aut_inte1(b,x,Qf,z,ph,par,target,hadron,charge,flavor):
    Q=math.sqrt(Qf)
    fsi=fsiver_all([x,b,target],par)
    ff=d1([z,b,hadron,charge])
    aut1=b**2*jv(1,b*ph/z)*evo([Q,b])*fsi[flavor]*ff[flavor]
    return aut1

def aut_inte2(b,x,Qf,z,ph,par,target,hadron,charge,flavor):
    Q=math.sqrt(Qf)
    pdf1=f1([x,b,target])
    ff=d1([z,b,hadron,charge])
    aut2=b*jv(0,b*ph/z)*evo([Q,b])*pdf1[flavor]*ff[flavor]
    return aut2



'''
par=[[0.75, 2.8, 203.0],[ -0.34, -3.9, -0.016],[ -0.7, 5.9, 0.36],[ 2.4, 0.64, -0.36]]

target='proton'
hadron='pi'
charge='+'
flavor=1
b=1.0
x=0.1
kt=2.
mu=2.0
print(fourier_mean_std([x, kt, mu, target, flavor]))
'''

'''
var1=[0.1, 3.18, 0.35, 0.141, target, hadron, charge]
var2=[0.043, 3.18, 0.35, 0.158, target, hadron, charge]

x=0.53
mu=18
fp=lha.f_1([x,mu,target])
print(fp)
'''


