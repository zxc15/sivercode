import sys
import math
import readdate as ex
import TMD as th
import lha
from iminuit import Minuit
import time


def fit(var):
    par_all=[var[0],var[1],var[2]]
    par_u=[var[3],var[4],var[5]]
    par_d=[var[6],var[7],var[8]]
    par_sea=[var[9],var[10],var[11]]
    
    print(var)
    par=[par_all,par_u,par_d,par_sea]
    experiment=['compass2009','compass2015','hermers','JLab2011','JLab2014']
    hadron=['K','pi','p']
    charge=['+','-','0']
    fit=0
    data_number=0
    for l in experiment:
        for j in hadron:
            for k in charge:
                ex_var=[l,j,k]
                test_ff=lha.d_1([0.2,1.2,j,k])
                if test_ff==0:#检验lha中是否有符合输入条件的碎裂函数
                    continue
                try:
                    exdata=ex.aut_ex(ex_var)#检验数据库中是否有符合输入条件的实验数据
                except:
                    continue
                #print(ex_var)
                target=' '
                if ex_var[0]=='compass2009':
                    target='deuteron'
                elif ex_var[0]=='compass2015' or ex_var[0]=='hermers':
                    target='proton'
                elif ex_var[0]=='JLab2011':
                    target='neutron'
                elif ex_var[0]=='JLab2014':
                    target='3He'
                else:
                    print('error: don not have '+ex_var[0]+' data')
                
                for i in exdata:
                    th_var=[i[0],i[1],i[2],i[3],target,ex_var[1],ex_var[2]]
                    th_par=par
                    err_ex=i[5]**2+i[6]**2+i[7]**2
                    fit=(th.aut_th(th_var,th_par)-i[4])**2/err_ex+fit
                    data_number=data_number+1
                    #print(fit)
                #print(data_number)   
    print(fit/(data_number-12))
    return fit

def fitcut(var):
    par_all=[var[0],var[1],var[2]]
    par_u=[var[3],var[4],var[5]]
    par_d=[var[6],var[7],var[8]]
    par_sea=[var[9],var[10],var[11]]
    print(var)
    par=[par_all,par_u,par_d,par_sea]
    experiment=['compass2009','JLab2011','JLab2014','hermers']
    hadron=['K','pi']
    charge=['+','-']

    fitcut=0
    data_number=0
    for l in experiment:
        for j in hadron:
            for k in charge:
                target=' '
                binn=''
                if l=='compass2009':
                    target='deuteron'
                    binn='ph'
                elif l=='hermers':
                    target='proton'
                    binn=''
                elif l=='JLab2011':
                    target='neutron'
                    binn='x'
                elif l=='JLab2014':
                    target='3He'
                    binn='x'
                else:
                    print('error: don not have '+ex_var[0]+' data')
                ex_var=[l,j,k,binn]
                test_ff=lha.d_1([0.2,1.2,j,k])
                if test_ff==0:#检验lha中是否有符合输入条件的碎裂函数
                    continue
                try:
                    exdata=ex.aut_excut(ex_var)#检验数据库中是否有符合输入条件的实验数据
                except:
                    continue
                if len(exdata)==0:
                    continue
                    
                #print(ex_var)
                
                
                for i in exdata:
                    th_var=[i[0],i[1],i[2],i[3],target,ex_var[1],ex_var[2]]
                    th_par=par
                    err_ex=i[5]**2+i[6]**2+i[7]**2
                    fitcut=(th.aut_th(th_var,th_par)-i[4])**2/err_ex+fitcut
                    data_number=data_number+1
                    #print(fitcut)
                #print(data_number)   
    print(fitcut/(data_number))
    return fitcut

def minuit(par):
    par_name=('r0', 'r1', 'r2', 'beta_u', 'eps_u', 'N_u', 'beta_d', 'eps_d', 'N_d', 'beta_sea', 'N_s', 'N_sea')
    m = Minuit(fitcut,par,name=par_name)
    m.errordef=1
    #m.fixed['r0', 'r1', 'r2', 'beta_u', 'eps_u', 'N_u', 'beta_d', 'eps_d', 'N_d', 'beta_sea', 'N_s', 'N_sea']=True
    m.limits=[(0, 14.7), (0, 62.8), (0, 4200), (-4.0, 4.0), (-44, 40), (-0.3, 3.0), (-8.0, 8.0), (-60, 65), (-4.0, 5.0), (-24.0, 27.0), (-7.0, 7.0), (-5.0, 4.0)]
    m.migrad()
    m.hesse()
    print(m)
    print(m.values)
    var=[*m.values]
    print('fit= %s'%fit(var))
    return m
par=[0.75, 2.8, 203.0, -0.34, -3.9, -0.016, -0.7, 5.9, 0.36, 2.4, 0.64, -0.36]
#par9=[1.07, 0, 21.2, -0.51, -5.07, -0.05, -0.18, -2.64,  2.3, 0,  1.44, -1.2]

#par9=[ 1.07165618e+00,  6.87293007e-03,  2.12011761e+01, -5.12497682e-01,
# -5.06878162e+00, -4.62877879e-02, -1.79534249e-01, -2.64023592e+00,
#  2.27697227e+00, -8.98372104e-03,  1.43659657e+00, -1.15718137e+00]

#par9=[ 1.072,  0.007,  21.201, -0.512, -5.069, -0.046, -.180, -2.640, 2.277, -0.009,  1.437, -1.157]
#par=[0.37, 0.03, 2.51, -0.48, -4.80, -0.02, -0.62, -1.06,  1.27,   1.55,  0.99, -0.57]
#par=[0.08, 1.14, 11.1, -0.48, -4.72, -0.012, -0.87, 2.55, 0.79, 1.74, 0.49, -0.27]
#par=[1.81, 2.8, 4060.07, -0.35, -3.85, -0.00304, -0.933, -1.20, 0.127, 3.50, 0.31, -0.254]
#par=[0.26, 12.0, 29.4, -0.38, -4.21, -0.007, -0.30, -2.20, 0.16, 3.37, 1.01, -0.81]

#fitcut(par)


start =time.time()
minuit(par)
end =time.time()
hour=(end-start)/3600
print('grep_Running time: %s hour' %hour)

