import pandas as pd
import numpy as np
import os

def aut_ex(var):#var=[experiment,hadron,charge,bin]

    skip=[0,1,2,3,4]
    binn=var[3]
    title=['x','Q^2','y','z','ph','w','clo','c1','c2','c3','siv','s1','s2','s3']
    path_head='/home/zengchunhua/workarea/python/siverdata/'
    exp=['compass2009','compass2015','JLab2011','JLab2014','hermers']
#   bin=['x', 'z', 'ph']
#   had=['K', 'pi', 'p']
#   cha=['+', '-', '0']
    
    if var[0]==exp[4]:
        path_head=path_head+'hermers/1D/'
        fil=path_head+var[1]+var[2]+'/'+binn+'.txt'
        
    else:
        for i in range(4):
            if var[0]==exp[i]:
                k=i
        
        if k==0 or k==1:
            expt='compass/'
            name1='/CS_'
            if k==1:
                name1='/S_'
                
        if k==2 or k==3:
            expt='JLab/'
            name1='/CS_'
        
        fil=path_head+expt+exp[k]+name1+var[1]+var[2]+'_'+binn+'.txt'
    
    
    
    read_data=pd.read_csv(fil, skiprows=skip, header=None, names=title,sep='\t')
    Qf_limit=0
    delta_max=0.8
    siv_data=[]
    for i in range(len(read_data)):
        data_set=[0,0,0,0,0,0,0,0]
        data_set[0]=read_data['x'][i]
        data_set[1]=read_data['Q^2'][i]
        data_set[2]=read_data['z'][i]
        data_set[3]=read_data['ph'][i]
        data_set[4]=read_data['siv'][i]
        data_set[5]=read_data['s1'][i]
        data_set[6]=read_data['s2'][i]
        data_set[7]=read_data['s3'][i]
        
        delta=data_set[3]/data_set[2]/(data_set[1]**0.5)
        Qf=data_set[1]

        if delta<delta_max and Qf>Qf_limit:
            siv_data.append(data_set)
        else:
            continue

    return siv_data

def aut_excut(var):#var=[experiment,hadron,charge,bin]

    skip=[0,1,2,3,4]
    binn=var[3]
    title=['x','Q^2','y','z','ph','w','clo','c1','c2','c3','siv','s1','s2','s3']
    path_head='/home/zengchunhua/workarea/python/siverdata/'
    exp=['compass2009','compass2015','JLab2011','JLab2014','hermers']
#   bin=['x', 'z', 'ph']
#   had=['K', 'pi', 'p']
#   cha=['+', '-', '0']

    Qf_limit=0
    delta_max=0.3
    siv_data=[]
    
    if var[0]==exp[4]:
        for k in range(4):
            for j in range(4):
                for i in range(4):
                    her3D=[i+1,j+1,k+1,var[1],var[2]]
                    data_set=read_hermers3D(her3D)
                    
                    delta=data_set[3]/data_set[2]/(data_set[1]**0.5)
                    Qf=data_set[1]
                    if delta<delta_max and Qf>Qf_limit:
                        siv_data.append(data_set)
                    else:
                        continue
                    
                            
    else:
        for i in range(4):
            if var[0]==exp[i]:
                k=i
        
        if k==0 or k==1:
            expt='compass/'
            name1='/CS_'
            if k==1:
                name1='/S_'
                
        if k==2 or k==3:
            expt='JLab/'
            name1='/CS_'
        
        fil=path_head+expt+exp[k]+name1+var[1]+var[2]+'_'+binn+'.txt'
        read_data=pd.read_csv(fil, skiprows=skip, header=None, names=title,sep='\t')        
        for i in range(len(read_data)):
            data_set=[0,0,0,0,0,0,0,0]
            data_set[0]=read_data['x'][i]
            data_set[1]=read_data['Q^2'][i]
            data_set[2]=read_data['z'][i]
            data_set[3]=read_data['ph'][i]
            data_set[4]=read_data['siv'][i]
            data_set[5]=read_data['s1'][i]
            data_set[6]=read_data['s2'][i]
            data_set[7]=read_data['s3'][i]
            
            delta=data_set[3]/data_set[2]/(data_set[1]**0.5)
            Qf=data_set[1]
            #print(delta)
            #print(Qf)
            if delta<delta_max and Qf>Qf_limit:
                siv_data.append(data_set)
            else:
                continue
    return siv_data
    
def read_hermers3D(var):#var=[x_number,z_number,ph_number,hardon,charge]
    skip=[0,1,2,3,4]
    title=['x','Q^2','y','z','ph','w','clo','c1','c2','c3','siv','s1','s2','s3']
    path_head='/home/zengchunhua/workarea/python/siverdata/hermers/3D/'
    xbin='x%d'%(var[0])
    zbin='z%d'%(var[1])
    siv_data=[]
    fil=path_head+var[3]+var[4]+'/'+xbin+'/'+zbin+'.txt'
    read_data=pd.read_csv(fil, skiprows=skip, header=None, names=title,sep='\t')
    for i in range(len(read_data)):
        data_set=[0,0,0,0,0,0,0,0]
        data_set[0]=read_data['x'][i]
        data_set[1]=read_data['Q^2'][i]
        data_set[2]=read_data['z'][i]
        data_set[3]=read_data['ph'][i]
        data_set[4]=read_data['siv'][i]
        data_set[5]=read_data['s1'][i]
        data_set[6]=read_data['s2'][i]
        data_set[7]=read_data['s3'][i]

        siv_data.append(data_set)
    return siv_data[var[2]-1] 

def data_random():
    newdata=[]
    
    experiment=['compass2009','JLab2011','JLab2014','hermers']
    hadron=['K','pi']
    charge=['+','-']


    for l in experiment:
        for j in hadron:
            for k in charge:
                binn=''
                if l=='compass2009':
                    binn='ph'
                elif l=='hermers':
                    binn=''
                elif l=='JLab2011':
                    binn='x'
                elif l=='JLab2014':
                    binn='x'
                else:
                    print('error: in aut_random')
                var=[l,j,k,binn]
                try:
                    olddata=aut_excut(var)#检验数据库中是否有符合输入条件的实验数据
                except:
                    continue
                if len(olddata)==0:
                    continue
                #print(var)
                for i in olddata:
                    data=np.random.normal(i[4],i[5])
                    ndata=[i[0],i[1],i[2],i[3],data,i[5],i[6],i[7]]
                    #print(ndata)
                    newdata.append(ndata)
                    
    #print(len(newdata))
    return newdata


def data_par():
    skip=1
    fil='/home/zengchunhua/workarea/test/mycode/use_value.txt'
    title=['r0', 'r1', 'r2', 'beta_u', 'eps_u', 'N_u', 'beta_d', 'eps_d', 'N_d', 'beta_sea', 'N_s', 'N_sea']
    read_data=pd.read_csv(fil, skiprows=skip, header=None, names=title,sep=' ')
    par=[]
    

    for i in range(len(read_data)):
        pari=[read_data[j][i] for j in title]
        par1=[pari[0], pari[1], pari[2]]
        par2=[pari[3], pari[4], pari[5]]
        par3=[pari[6], pari[7], pari[8]]
        par4=[pari[9], pari[10], pari[11]]
        parr=[par1, par2, par3, par4]
        par.append(parr)
    
    return par

 
    
    
    
