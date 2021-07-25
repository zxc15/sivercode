import sys
import math
import readdate as ex
import TMD as th
from iminuit import Minuit
import time

ex_var=['JLab2014','K','-']
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
	

def fit(var):
    beta=var[0]
    eps=var[1]
    Np=var[2]
    r0=var[3]
    r1=var[4]
    r2=var[5]
    exdata=ex.aut_ex(ex_var)
    fit=0
    
    for i in exdata:
        th_var=[i[0],i[1],i[2],i[3],target,ex_var[1],ex_var[2]]
        th_par=var
        err_ex=i[5]**2+i[6]**2+i[7]**2
        fit=(th.aut_th(th_var,th_par)-i[4])**2/err_ex+fit
    return fit

def minuit(par):
    m = Minuit(fit,par,name=('beta','eps','Np','r0','r1','r2'))
    m.errordef=1
    #m.fixed['beta','Np','r0','r1']=True
    m.limits=[(0,None),(0,None),(0,None),(0,None),(0,None),(0,None)]
    m.migrad()
    m.hesse()
    print(m.values)
    return m

#par=[0.2,0.2,0.01,3.0,4.0,0.2]
#par=[15.03, 22.42, 0.11, 28.72, 148.29, 6.98]
par=[15.03, 113.29, 0.11, 28.72, 148.29, 0.33]

#par=[10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
#par=[15.03,22.42,0.11,3.0,4.0,0.2]
start =time.time()
minuit(par)
end =time.time()
print('Running time: %s Seconds'%(end-start))




