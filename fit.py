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
    par_s=[var[9],var[10],var[11]]
    par_ub=[var[12],var[13],var[14]]
    par_db=[var[15],var[16],var[17]]
    par_sb=[var[18],var[19],var[20]]
	
    par=[par_all,par_u,par_d,par_s,par_ub,par_db,par_sb]
    experiment=['compass2009','compass2015','hermers','JLab2011']#暂不考虑3He为靶的情况
    hadron=['K','pi','p']
    charge=['+','-','0']
    fit=0
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
                    #print(fit)
    print(fit)
    return fit

def minuit(par):
    par_name=('r0','r1','r2','beta_u', 'eps_u', 'Nq_u','beta_d', 'eps_d', 'Nq_d','beta_s', 'eps_s', 'Nq_s', 'beta_ub', 'eps_ub', 'Nq_ub', 'beta_db', 'eps_db', 'Nq_db', 'beta_sb', 'eps_sb', 'Nq_sb')
    m = Minuit(fit,par,name=par_name)
    m.errordef=1
    #m.fixed['beta','eps','Np','r0','r1','r2']=True
    m.limits=[(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None)]
    m.migrad()
    m.hesse()
    print(m.values)
    var=[*m.values]
    print('fit= %s'%fit(var))
    return m


par=[2.0,2.0,2.0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
#print(fit(par))
start =time.time()
minuit(par)
end =time.time()
print('Running time: %s Seconds'%(end-start))




