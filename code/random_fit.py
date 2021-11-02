import sys
import math
import readdate as ex
import TMD as th
import lha
from iminuit import Minuit
import time

exdata=ex.data_random()
par1=[0.75, 2.8, 203.0, -0.34, -3.9, -0.016, -0.7, 5.9, 0.36, 2.4, 0.64, -0.36]
par4=[1.4, 0, 104, -0.4, -4.22, -0.017, -0.02, -2.58, 0.15, 0.04, 0.17, -0.12]

def random_fit(var):
    print(var)
    par_all=[var[0],var[1],var[2]]
    par_u=[var[3],var[4],var[5]]
    par_d=[var[6],var[7],var[8]]
    par_sea=[var[9],var[10],var[11]]
    par=[par_all,par_u,par_d,par_sea]
    data_number=0
    random_fit=0
    for i in exdata:
        if data_number==0:
            target, hadron, charge='deuteron', 'K', '+'
        elif data_number==1:
            target, hadron, charge='deuteron', 'K', '-'
        elif data_number==2:
            target, hadron, charge='deuteron', 'pi', '+'
        elif data_number==3:
            target, hadron, charge='deuteron', 'pi', '-'
        elif data_number==4:
            target, hadron, charge='neutron', 'pi', '+'
        elif data_number==5:
            target, hadron, charge='neutron', 'pi', '-'
        elif data_number==6:
            target, hadron, charge='3He', 'K', '+'
        elif data_number>6 and data_number<19:
            target, hadron, charge='proton', 'K', '+'
        elif data_number>18 and data_number<31:
            target, hadron, charge='proton', 'K', '-'
        elif data_number>30 and data_number<42:
            target, hadron, charge='proton', 'pi', '+'
        elif data_number>41 and data_number<53:
            target, hadron, charge='proton', 'pi', '-'
        else:
            print('error: in random_fit' )
        data_number=data_number+1
        var=[i[0], i[1], i[2], i[3], target, hadron, charge]
        print(target, hadron, charge)
        err_ex=i[5]**2
        random_fit=(th.aut_th(var,par)-i[4])**2/err_ex+random_fit
        print(random_fit)
    print(data_number)
    print(random_fit/(data_number))
    return random_fit

def random_minuit(par):
    par_name=('r0', 'r1', 'r2', 'beta_u', 'eps_u', 'N_u', 'beta_d', 'eps_d', 'N_d', 'beta_sea', 'N_s', 'N_sea')
    m = Minuit(random_fit,par,name=par_name)
    m.errordef=1
    #m.fixed['r0', 'r1', 'r2', 'beta_u', 'eps_u', 'N_u', 'beta_d', 'eps_d', 'N_d', 'beta_sea', 'N_s', 'N_sea']=True
    m.limits=[(0, 14.0), (0, 200.8), (0, 14200), (-4.0, 4.0), (-44, 40), (-0.3, 3.0), (-8.0, 8.0), (-60, 65), (-4.0, 5.0), (-24.0, 27.0), (-7.0, 7.0), (-5.0, 4.0)]
    m.migrad()
    m.hesse()
    print('grep_par','r0', 'r1', 'r2', 'beta_u', 'eps_u', 'N_u', 'beta_d', 'eps_d', 'N_d', 'beta_sea', 'N_s', 'N_sea')
    print('grep_val', *m.values)
    print('grep_err', *m.errors)
    var=[*m.values]
    print('grep_fit= %s'%random_fit(var))
    return m


start =time.time()
random_minuit(par1)
end =time.time()
print('grep_Running time: %s Seconds'%(end-start))

