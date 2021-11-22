import numpy as np
import scipy as sp
import readdate as ex
import evolution as evl
#import anomalous_dimension as ad
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

inte_epsabs=1.49e-9
inte_epsrel=1.49e-5
inte_limit=1000
inte_limlst=50


lad1=0.198
lad2=9.3
lad3=431.
lad4=2.12
lad5=-4.44

eta1=0.26
eta2=0.476
eta3=0.478
eta4=0.483

Nf=evl.Nf
CF=evl.CF
BNP=evl.BNP

flavor_use=[2, 1, 3, -2, -1, -3]

qq=2.
qg=1.

#m_target={'deuteron':1.876, 'proton':0.938, 'neutron':0.940, '3He':2.808}
m_target={'deuteron':0.939, 'proton':0.939, 'neutron':0.939, '3He':0.939}

def inte(fun, a, b, par):
    if par==0:
        inte=integrate.quad(fun, a, b, epsabs=inte_epsabs, epsrel=inte_epsrel, limit=inte_limit, limlst=inte_limlst)[0]
    else:
        inte=integrate.quad(fun, a, b, args=par, epsabs=inte_epsabs, epsrel=inte_epsrel, limit=inte_limit, limlst=inte_limlst)[0]
    return inte

#..........................TMD distribution...............
def mu0(b):
    #bstar=b/(1.+b**2./BNP**2.)**0.5
    mu0=2.*math.exp(-gammaE)/b+2.
    return mu0

def f1(var):#TMDpdf
    
    x=var[0]
    b=var[1]
    target=var[2]
    fp=lha.f_1([x,mu0(b),target])
    fnp=f_np([x,b])
    #c1=fc1(var)
    c1=0
    f1={i:fp[i]*fnp+c1 for i in flavor_use}
    return f1

def f_np(var):
    [x, b]=var
    f_np=lad1*(1.-x)+lad2*x+x*(1-x)*lad5
    f_np=-f_np*b**2./(1.+lad3*x**(lad4)*b**2.)**0.5
    f_np=math.exp(f_np)
    return f_np

def fc1(var):
    [x, b, target]=var
    '''
    fc1=c1_p_qg(var)+c1_c_qg(var)
    for i in flavor_use:
        par=[x, b, target, i]
        fc1=p1_q(par)+p2_q(par)+c1_c_qq(par)+fc1
    '''
    fc1=c1_p_qg(var)+c1_c_qg(var)
    for i in flavor_use:
        par=[x, b, target, i]
        fc1=p1_q(par)+p2_q(par)+c1_c_qq(par)+fc1

    return fc1

def c1_p_qg(var):
    [x, b, target]=var
    mui=mu0(b)
    a_s=lha.a_s(mui)
    bstar=b/(1.+b**2./BNP**2.)**0.5
    L=evl.lmu(mui, b)
    allpar=(x, b, target)
    c1=-a_s*L*inte(c1_p_qg_inte, x, 1., allpar)
    return c1

def c1_c_qg(var):
    [x, b, target]=var
    mui=mu0(b)
    a_s=lha.a_s(mui)
    allpar=(x, b, target)
    c1=a_s*inte(c1_c_qg_inte, x, 1., allpar)
    return c1

def p1_q(var):
    [x, b, target, flavor]=var
    mui=mu0(b)
    a_s=lha.a_s(mui)
    bstar=b/(1.+b**2./BNP**2.)**0.5
    L=evl.lmu(mui, b)
    allpar=(x, b, target, flavor)
    fp2=lha.f_1([x, mui, target])[flavor]
    fnp=f_np([x,b])
    c1=-qq*a_s*L*CF*(inte(p1_q_inte, x, 1., allpar)+4.*np.log(1.-x)*fp2*fnp)
    return c1

def p2_q(var):
    [x, b, target, flavor]=var
    mui=mu0(b)
    a_s=lha.a_s(mui)
    bstar=b/(1.+b**2./BNP**2.)**0.5
    L=evl.lmu(mui, b)
    allpar=(x, b, target, flavor)
    fp=lha.f_1([x, mui, target])[flavor]
    fnp=f_np([x,b])
    c1=-qq*a_s*L*CF*(inte(p2_q_inte, x, 1., allpar)+3.*fp*fnp)
    return c1

def c1_c_qq(var):
    [x, b, target, flavor]=var
    mui=mu0(b)
    a_s=lha.a_s(mui)
    allpar=(x, b, target, flavor)
    fp=lha.f_1([x, mui, target])[flavor]
    fnp=f_np([x,b])
    c1=a_s*CF*(inte(c1_c_qq_inte, x, 1., allpar)-pi**2./6.*fp*fnp)
    return c1

def c1_p_qg_inte(y, x, b, target):
    mui=mu0(b)
    fp=lha.f_1([x/y, mui, target])[0]
    fnp=f_np([x,b])
    c1=2.*qg*(1.-2.*y+2.*y**2.)*fp*fnp/y
    return c1

def c1_c_qg_inte(y, x, b, target):
    mui=mu0(b)
    fp=lha.f_1([x/y, mui, target])[0]
    fnp=f_np([x,b])
    c1=2.*(1.-y)*fp*fnp
    return c1

def p1_q_inte(y, x, b, target, flavor):
    mui=mu0(b)
    fp1=lha.f_1([x/y, mui, target])[flavor]
    fp2=lha.f_1([x, mui, target])[flavor]
    fnp=f_np([x,b])
    c1=4.*fnp*(fp1/y-fp2)/(1.-y)
    return c1
    
def p2_q_inte(y, x, b, target, flavor):
    mui=mu0(b)
    fp1=lha.f_1([x/y, mui, target])[flavor]
    fp2=lha.f_1([x, mui, target])[flavor]
    fnp=f_np([x,b])
    c1=2.*(-1.-y)*fp1*fnp/y
    return c1

def c1_c_qq_inte(y, x, b, target, flavor):
    mui=mu0(b)
    fp=lha.f_1([x/y, mui, target])[flavor]
    fnp=f_np([x,b])
    c1=2.*(1.-y)*fp*fnp/y
    return c1


def dp_fourier_f1(var):
    [x, kt, mu, target, flavor]=var
    #....................fc0 傅里叶变换....................
    def fc0b_inte(b):
        fp=lha.f_1([x,mu0(b),target])
        fnp=f_np([x,b])
        fc0b=fp[flavor]*fnp 
        fc0b_inte=b/(2.*pi)*jv(0,b*kt)*fc0b*evo([mu, b])**.5
        return fc0b_inte
    fc0=inte(fc0b_inte, 0, evl.bmax, 0)
    
    #....................fc1 傅里叶变换.....................
    
    def fc1b_g_inte(y, b):
        mui=mu0(b)
        a_s=lha.a_s(mui)
        L=evl.lmu(mui, b)
        fp=lha.f_1([x/y, mui, target])[0]
        fnp=f_np([x,b])
        
        P_inte=-a_s*L*2.*Nf*(1.-2.*y+2.*y**2.)*fp*fnp/y
        C_inte=a_s*2.*(1.-y)*fp*fnp
        PC=P_inte+C_inte
        fc1b_g_inte=b/(2.*pi)*jv(0,b*kt)*PC*evo([mu, b])**.5
        return fc1b_g_inte
    fc1_g=integrate.nquad(fc1b_g_inte, [[x, 1.], [0, evl.bmax]])[0]
    
    def fc1_q(fla):
        def fc1b_q_inte(y, b):
            mui=mu0(b)
            a_s=lha.a_s(mui)
            L=evl.lmu(mui, b)
            fp1=lha.f_1([x/y, mui, target])[fla]
            fp2=lha.f_1([x, mui, target])[fla]
            fnp=f_np([x,b])
            
            P1_inte=-a_s*L*CF*4.*fnp*(fp1/y-fp2)/(1.-y)
            P2_inte=-a_s*L*CF*(2.*(-1.-y)*fp1*fnp/y+3.*fp2*fnp)
            C_inte=a_s*CF*(2.*(1.-y)/y*fp1*fnp-pi**2./6.*fp2*fnp)
            PC=P1_inte+P2_inte+C_inte
            fc1_q_inte=b/(2.*pi)*jv(0,b*kt)*PC*evo([mu, b])**.5
            return fc1_q_inte
        
        def fc1b_q_inte0(b):
            mui=mu0(b)
            a_s=lha.a_s(mui)
            L=evl.lmu(mui, b)
            fp2=lha.f_1([x, mui, target])[fla]
            fnp=f_np([x,b])
            
            P1_inte0=-a_s*L*CF*4.*np.log(1.-x)*fp2*fnp
            fc1_q_inte0=b/(2.*pi)*jv(0,b*kt)*P1_inte0*evo([mu, b])**.5
            return fc1_q_inte0
            
        fc1_q=integrate.nquad(fc1b_q_inte, [[x, 1.], [0, evl.bmax]])[0]
        fc1_q=fc1_q+inte(fc1b_q_inte0, 0, evl.bmax, 0)
        return fc1_q
    
    fc1=fc1_g
    for i in flavor_use:
        fc1=fc1_q(i)+fc1
   
        
    fc=fc0+fc1
    
    print(fc)
    return fc


def d1(var):#TMDff
    
    z=var[0]
    b=var[1]
    hardon=var[2]
    charge=var[3]
    
    dp=lha.d_1([z,mu0(b),hardon,charge])
    dnp=d_np([z,b])
    d1={i:dp[i]*dnp/z**2. for i in flavor_use}
    return d1

def d_np(var):
    [z, b]=var
    d_np=(eta1*z+eta2*(1.-z))*b**2./(1.+eta3*(b/z)**2.)**0.5/z**2.
    d_np=math.exp(-d_np)*(1.+eta4*b**2./z**2.)
    return d_np

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
        fourier_fsiver_all[j]=inte(fourier_inte, 0., evl.bmax, allpar)
        
    return fourier_fsiver_all
    
def fourier_f1_all(var):
    [x, kt, mu, target]=var
    #print(var)
    fourier_f1_all={}
    
    for j in flavor_use:
        allpar=(x, kt, mu, target, j)
        fourier_f1_all[j]=inte(fourier_f1_inte, 0, evl.bmax, allpar)
        #print(fourier_f1_all[j])
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
        inte1=inte(aut_inte1,0,pole,allpar)+inte(aut_inte1, pole,evl.bmax ,allpar)
        inte2=inte(aut_inte2,0,pole,allpar)+inte(aut_inte2, pole,evl.bmax ,allpar)
        
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

#print(p1_q([0.2, 1, 'proton', 3]))

x, kt, mu, target, flavor=0.1, 0.001, 2., 'proton', 2
#var=[x,kt,mu,target, flavor]
#dp_fourier_f1(var)
#allpar=(x, kt, mu, target, flavor)
#print(inte(fourier_f1_inte, 0, evl.bmax, allpar))
#print(fourier_f1_all([x,kt,mu,target]))


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


