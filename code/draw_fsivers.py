import numpy as np
import readdate as ex
import time
import lha
import evolution as evl
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import math
import TMD as th


gammaE=np.euler_gamma
BNP=evl.BNP
par1=[[0.75, 2.8, 203.0],[ -0.34, -3.9, -0.016],[ -0.7, 5.9, 0.36],[ 2.4, 0.64, -0.36]]
par2=[[7.2,640,66184],[-0.248,-3.31,-0.101],[-0.925,3.4,2.22],[4.5,1.5,-3.1]]
par3=[[2.64, 33.93, 4198.16],[-0.047, -3.06, -0.0029],[-0.71, 5.06, 0.23],[4.27, 0.65, -0.72]]
par4=[[1.4, 0, 104],[-0.4, -4.22, -0.017],[-0.02, -2.58, 0.15],[0.04, 0.17, -0.12]]
par5=[[1.75, 24.9, 32.44], [-0.38, -4.21, -0.012], [-0.31, -3.66, 0.18], [0.53, -1.46, -3.33]]
par6=[[.26, 0. ,309.73],[-.364, -4.041, -.005],[-0.92, 5.53, .285],[2.27, .279, -.11]]
par7=[[0.04, 2.65, 93.69],[-0.48, -4.75, -0.011],[-0.93, 3.78, 0.74],[1.73, 0.31, -0.17]]#LL
par8=[[0.08, 1.14, 11.1],[-0.48, -4.72, -0.012],[-0.87, 2.55, 0.79],[1.74, 0.49, -0.27]]#NNLL
#par9=[[0.1, 0.04, 0.02],[-0.47, -4.77, -0.005],[-0.81, 2.58, 0.63],[1.76, 0.5, -0.29]]#NLL par1
#par9=[[1.07, 0, 21.2], [-0.51, -5.07, -0.05], [-0.18, -2.64,  2.3], [0,  1.44, -1.2]]#NLL par8
#par9=[[ 1.072,  0.007,  21.201],[ -0.512, -5.069, -0.046], [-.180, -2.640, 2.277], [-0.009,  1.437, -1.157]]#NLLpar8
par10=[[0.95, 0.09, 195.],[-0.35, -3.9, -0.013],[-0.77, 1.8, 0.34],[2.3, 0.43, -0.23]]
def draw_fsi1(par):   
    b=np.linspace(0,4,50)
    y=[]
    for i in b:
        var=[0.1,i,'proton']
        fsi=th.fsiver_all(var,par)
        y.append(fsi[1])
    
    y=[math.log(y[i]) for i in range(len(y))]
    
    
    var=[0.1,0.5,'proton']
    print(th.fsiver_all(var,par)[1])

    plt.figure()
    
    #x,y,轴刻度与说明设置
    #plt.xlim((0, 4))
    plt.ylim((-5.5, 0))
    plt.xlabel('b(GeV$^{-1}$)', fontsize=14)
    #plt.ylabel(r'$f^{\perp}_{T;q\to u}$')
    x_ticks=np.linspace(0,4,11)
    y_ticks=[0.005,0.01,0.05,0.1,0.5,1]
    ylog=[math.log(i) for i in y_ticks]
    plt.xticks(x_ticks)
    plt.yticks(ylog,y_ticks)
    
    #坐标轴边框设置(颜色和位置)
    ax=plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_position(('data',-5.5))
    ax.spines['left'].set_position(('data',0))
    #ax.xaxis.set_ticks_position('bottom')
    #ax.yaxis.set_ticks_position('left')
    
    #图例设置
    plt.plot(b,y,label=r'f$^{\perp}_{T;q\to u}$(0.1,b)')
    plt.legend(loc='best', fontsize=14)
    
    plt.show()
    return 0

def draw_fsi2(par):
    X=np.arange(0.001,1, 0.05)
    #X=np.arange(0.2,0.8, 0.05)
    B=np.arange(0.01, 3, 0.1)
    
    new_ticks=np.array([math.log(i) for i in X])

    x, b=np.meshgrid(X,B)
    x_ticks, b_ticks=np.meshgrid(new_ticks,B)
    
    (yn,xn)=x.shape

    #Zd=[]
    Zu=[]
    for i in range(yn):
        Zyd=[]
        Zyu=[]
        for j in range(xn):
            var=[x[i][j],b[i][j],'proton']
            #var=[x[i][j],b[i][j],'pi','+']
            #fsi=th.fsiver_all(var,par)
            fsi=th.f1(var)
            #Zyd.append(x[i][j]*fsi[1])
            Zyu.append(x[i][j]*fsi[2])
            #print(Zy)
        #Zd.append(Zyd)
        Zu.append(Zyu)
        
    #Zd=np.array(Zd)
    Zu=np.array(Zu)
    print(Zu)
     
    '''
    #画图fsiver_d<-p
    fig=plt.figure()
    fig.suptitle(r'f$^{\perp}_{T;q\to d}$', fontsize=14)
    ax=Axes3D(fig)
    ax.plot_surface(x_ticks, b_ticks, Zd, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
    
    #cbar=plt.contourf(x_ticks, b_ticks, Zd, cmap=plt.get_cmap('rainbow'))
    #fig.colorbar(cbar)
    #plt.contourf(x_ticks, b_ticks, Zd, cmap=plt.get_cmap('rainbow'))
    
    #ax.set_zlim(-1, 4)
    plt.xticks([math.log(0.01),math.log(0.1),math.log(1)],[0.01, 0.1, 1])
    plt.xlabel('x', fontsize=14)
    plt.ylabel('b(GeV)$^{-1}$', fontsize=14)
    '''
    
    #画图-fsiver_u<-p
    fig=plt.figure(num=2)
    #fig.suptitle(r'-f$^{\perp}_{T;q\to u}$', fontsize=14)
    fig.suptitle(r'$xF_{u\leftarrow p}$', fontsize=14)
    
    ax=Axes3D(fig)
    ax.plot_surface(x_ticks, b_ticks, Zu, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
    #ax.plot_surface(x, b, Zu, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
    '''
    cbar=plt.contourf(x_ticks, b_ticks, Zu, cmap=plt.get_cmap('rainbow'))
    fig.colorbar(cbar)
    plt.contourf(x_ticks, b_ticks, Zu, cmap=plt.get_cmap('rainbow'))
    '''
    
    #ax.set_zlim(-1, 4)
    plt.xticks([math.log(0.01),math.log(0.1),math.log(1)],[0.01, 0.1, 1])
    plt.xlabel('x', fontsize=14)
    plt.ylabel('b(GeV)$^{-1}$', fontsize=14)
    
    plt.show()
    
    return 0

def draw_f1(var):
    flavor_num=var
    X=np.arange(0.001,1, 0.05)
    #X=np.arange(0.2,0.8, 0.05)
    B=np.arange(0.01, 3, 0.1)
    
    flavor={2:'u', 1:'d', 3:'s', -2:'sea', -1:'sea', -3:'sea'}
    nam=flavor[flavor_num]
    
    new_ticks=np.array([math.log(i) for i in X])

    x, b=np.meshgrid(X,B)
    x_ticks, b_ticks=np.meshgrid(new_ticks,B)
    
    (yn,xn)=x.shape

    Z=[]
    for i in range(yn):
        Zy=[]
        for j in range(xn):
            var=[x[i][j],b[i][j],'proton']
            #var=[x[i][j],b[i][j],'pi','+']
            #f1=th.fsiver_all(var,par)
            f1=th.f1(var)
    
            Zy.append(x[i][j]*f1[flavor_num])
        Z.append(Zy)
        
    Z=np.array(Z)
    print(Z)
   
    
    #画图-fsiver_u<-p
    fig=plt.figure(num=2)
    fig.suptitle(r'$xF_{%s \leftarrow p}$' %nam, fontsize=14)
    
    ax=Axes3D(fig)
    ax.plot_surface(x_ticks, b_ticks, Z, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
    #ax.plot_surface(x, b, Zu, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
   
    
    #ax.set_zlim(-1, 4)
    plt.xticks([math.log(0.01),math.log(0.1),math.log(1)],[0.01, 0.1, 1])
    plt.xlabel('x', fontsize=14)
    plt.ylabel('b(GeV)$^{-1}$', fontsize=14)
    
    plt.show()
    
    return 0

def draw_aut(var, ax, par):
    [experiment,hadron,charge,binn]=var
    target=' '
    
    bin_num=0
    if binn=='x':
        bin_num=0
    elif binn=='z':
        bin_num=2
    elif binn=='ph':
        bin_num=3
    else:
        print('error: don not have bin '+binn)


    ylim_up=0.39
    ylim_down=-0.39
    if experiment=='compass2009':
        target='deuteron'
        ylim_up=0.14
        ylim_down=-0.14
    elif experiment=='compass2015' or experiment=='hermers':
        target='proton'
    elif experiment=='JLab2011':
        target='neutron'
    elif experiment=='JLab2014':
        target='3He'
    else:
        print('error: don not have '+experiment+' data')
    
    exdata=ex.aut_ex(var)
    
    cutdata=0

    x=[exdata[i+cutdata][bin_num] for i in range(len(exdata)-cutdata)]

    yth=[th.aut_th([exdata[i+cutdata][0],exdata[i+cutdata][1],exdata[i+cutdata][2],exdata[i+cutdata][3],target,hadron,charge],par) for i in range(len(exdata)-cutdata)]

    yex=[exdata[i+cutdata][4] for i in range(len(exdata)-cutdata)]
    err=[math.sqrt(exdata[i+cutdata][5]**2) for i in range(len(exdata)-cutdata)]

    low_error=err
    upper_error=err
    error_limit=[low_error,upper_error]
    
    final_state_hadron=' '
    if hadron+charge=='K+':
        final_state_hadron='K$^+$'
    elif hadron+charge=='K0':
        final_state_hadron='K$^0$'
    elif hadron+charge=='K-':
        final_state_hadron='K$^-$'
    elif hadron+charge=='pi+':
        final_state_hadron=r'$\pi^+$'
    elif hadron+charge=='pi-':
        final_state_hadron=r'$\pi^-$'
    elif hadron+charge=='pi0':
        final_state_hadron=r'$\pi^0$'
    elif hadron+charge=='p+':
        final_state_hadron='$p$'
    elif hadron+charge=='p-':
        final_state_hadron=r'$\bar{p}$'
    else:
        print('error: final state hadron input error')
        return 0

    tit='final state hadron: '+final_state_hadron
    if len(yth)<20:
        ax.errorbar(x,yth,xerr=0.02,fmt='.')
    else:
        ax.plot(x,yth)
    plt.axhline(y=0, color='k')
    ax.set_title(tit)
    ax.set_xlabel(binn, fontsize=12)
    ax.set_ylabel('A$_{UT}$', fontsize=12)
    delta1=0.3
    delta2=0.8
    for i in range(len(exdata)):
        delta=exdata[i][3]/exdata[i][2]/(exdata[i][1]**0.5)
        if delta<delta2:
            mfc_c='r'
            if delta<delta1:
                mfc_c='k'
        else:
            mfc_c='w'
        ax.errorbar(x[i],yex[i],yerr=err[i], fmt='o',ms=8,mfc=mfc_c, mec='k',ecolor='#FF8C00', capsize=5)
    #ax.errorbar(x,yex,yerr=error_limit,fmt='o',ms=8,mfc='c',mec='r',capsize=10)

    #ax.set_xlim(0.1,0.64)
    ax.set_ylim(ylim_down, ylim_up)
    return 0

def draw_AUT(experiment, par, binn):
    target=' '
    if experiment=='compass2009':
        target='deuteron'
        nr, nc=3, 2
        var1=['compass2009','K','+', binn]
        var2=['compass2009','pi','+', binn]
        var3=['compass2009','K','-', binn]
        var4=['compass2009','pi','-', binn]
        var5=['compass2009','K','0', binn]
        var=[var1,var2,var3,var4,var5]
    elif experiment=='compass2015':
        target='proton'
        nr, nc=2, 2
        var1=['compass2015','K','+', binn]
        var2=['compass2015','pi','+',binn]
        var3=['compass2015','K','-', binn]
        var4=['compass2015','pi','-', binn]
        var=[var1,var2,var3,var4]
    elif experiment=='hermers':
        target='proton'
        nr, nc=4, 2
        var1=['hermers','K','+', binn]
        var2=['hermers','p','+', binn]
        var3=['hermers','K','-', binn]
        var4=['hermers','p','-', binn]
        var5=['hermers','pi','+', binn]
        var6=['hermers','pi','-', binn]
        var7=['hermers','pi','0', binn]
        var=[var1,var2,var3,var4,var5,var6,var7]
    elif experiment=='JLab2011':
        target='neutron'
        nr, nc=1, 2
        var=[['JLab2011','pi','+', binn],['JLab2011','pi','-', binn]]
    elif experiment=='JLab2014':
        target='$^3$He'
        nr, nc=1, 2
        var=[['JLab2014','K','+', binn],['JLab2014','K','-', binn]]
    else:
        print('error: don not have '+experiment+' data')
  
    #fig,axs = plt.subplots(nrows=nr, ncols=nc, constrained_layout=True, num=3)
    
    fig=plt.figure(num=3)
    fig.suptitle(experiment+'    '+'target:'+target)
    
    for i in range(len(var)):
        ax=plt.subplot(nr,nc,i+1)
        if (i % 2)==0:
            xch=0
        else:
            xch=0.05
        pos=ax.get_position()
        ax.set_position([pos.x0+xch, pos.y0, pos.width, pos.height*0.8])
        draw_aut(var[i],ax,par)
    
    #plt.tight_layout(pad=0.4, w_pad=None, h_pad=None)
    plt.show()
    return 0

def draw_exdata(target, binn):
    experiment=[]
    if target=='proton':
        experiment=['compass2015','hermers']
    elif target=='neutron':
        experiment=['JLab2011']
    elif target=='deuteron':
        experiment=['compass2009']
    elif target=='3He':
        experiment=['JLab2014']
    else:
        print('error:do not have target:'+target )
        return 0

    hadron=['K','pi','p']
    charge=['+','-','0']
    
    for i in experiment:
        for j in hadron:
            for k in charge:
                var=[i,j,k,binn]

                try:
                    exdata=ex.aut_ex(var)
                except:
                    continue

                print(var)
                x=[exdata[l][0] for l in range(len(exdata))]
                y=[exdata[l][1] for l in range(len(exdata))]
                   
                #err=[math.sqrt(exdata[i][5]**2+exdata[i][4]**2+exdata[i][4]**2) for i in range(len(exdata))]
                err=[i*0 for i in range(len(exdata))]
                    
                low_error=err
                upper_error=err
                error_limit=[low_error,upper_error]
                
                mark_style=''
                if i=='compass2009':
                    mark_style='^'
                elif i=='compass2015':
                    mark_style='v'
                elif i=='hermers':
                    mark_style='o'
                elif i=='JLab2011':
                    mark_style='<'
                elif i=='JLab2014':
                    mark_style='>'
                else:
                    print('error: don not have '+i+' data')
           
                    
                mark_color=''
                if j=='K':
                    mark_color='r'
                elif j=='pi':
                    mark_color='b'
                elif j=='p':
                    mark_color='k'
                else:
                    print('error: don not have '+j+' data')
                    
                
                bar_color=''
                if k=='+':
                    bar_color='r'
                elif k=='-':
                    bar_color='g'
                elif k=='0':
                    bar_color='k'
                else:
                    print('error: don not have '+k+' data')
                
                final_state_hadron=' '
                if j+k=='K+':
                    final_state_hadron='K$^+$'
                elif j+k=='K0':
                    final_state_hadron='K$^0$'
                elif j+k=='K-':
                    final_state_hadron='K$^-$'
                elif j+k=='pi+':
                    final_state_hadron=r'$\pi^+$'
                elif j+k=='pi-':
                    final_state_hadron=r'$\pi^-$'
                elif j+k=='pi0':
                    final_state_hadron=r'$\pi^0$'
                elif j+k=='p+':
                    final_state_hadron='$p$'
                elif j+k=='p-':
                    final_state_hadron=r'$\bar{p}$'
                else:
                    print('error: final state hadron input error')

                plt.errorbar(x,y,yerr=error_limit,fmt=mark_style,ms=6,mfc=mark_color,mec=mark_color,ecolor=bar_color,capsize=10,label=i+':'+final_state_hadron)
    plt.title(target)            
    plt.legend(loc='best',fontsize=9)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('A$_{UT}$', fontsize=12)
    #plt.xlim(0.5,0.95)
    #plt.ylim(-1.8,1.8)
    plt.show()
    return 0

def draw_hermers3D(var,par):
    [hadron, charge, z_num, ph_num]=var
    exdata=[ex.read_hermers3D([i+1, z_num, ph_num, hadron, charge]) for i in range(4)]
    target='proton'
    
    x=[exdata[i][0] for i in range(len(exdata))]

    yth=[th.aut_th([exdata[i][0],exdata[i][1],exdata[i][2],exdata[i][3],target,hadron,charge],par) for i in range(len(exdata))]

    yex=[exdata[i][4] for i in range(len(exdata))]
    err=[math.sqrt(exdata[i][5]**2+exdata[i][6]**2+exdata[i][7]**2) for i in range(len(exdata))]
    #err=[math.sqrt(exdata[i][5]**2) for i in range(len(exdata))]

    final_state_hadron=' '
    ylim_up=0.1
    ylim_down=-0.1
    if hadron+charge=='K+':
        final_state_hadron='K$^+$'
        ylim_up=0.23
        ylim_down=-0.07
    elif hadron+charge=='K0':
        final_state_hadron='K$^0$'
    elif hadron+charge=='K-':
        final_state_hadron='K$^-$'
        ylim_up=0.26
        ylim_down=-0.26
    elif hadron+charge=='pi+':
        final_state_hadron=r'$\pi^+$'
        ylim_up=0.14
        ylim_down=-0.09
    elif hadron+charge=='pi-':
        final_state_hadron=r'$\pi^-$'
        ylim_up=0.14
        ylim_down=-0.09
    elif hadron+charge=='pi0':
        final_state_hadron=r'$\pi^0$'
    elif hadron+charge=='p+':
        final_state_hadron='$p$'
    elif hadron+charge=='p-':
        final_state_hadron=r'$\bar{p}$'
    else:
        print('error: final state hadron input error')
        return 0

    tit='final state hadron: '+final_state_hadron
    if len(yth)<20:
        plt.errorbar(x,yth,xerr=0.02,fmt='.')
    else:
        plt.plot(x,yth)
    plt.axhline(y=0, color='k')
    plt.title(tit)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('A$_{UT}$', fontsize=12)
    delta1=0.3
    delta2=0.8
    for i in range(4):
        delta=exdata[i][3]/exdata[i][2]/(exdata[i][1]**0.5)
        if delta<delta2:
            mfc_c='r'
            if delta<delta1:
                mfc_c='k'
        else:
            mfc_c='w'
        plt.errorbar(x[i],yex[i],yerr=err[i],fmt='o',ms=12,mfc=mfc_c, mec='k',ecolor='#FF8C00', capsize=5)

    #plt.xlim(0.1,0.64)
    plt.ylim(ylim_down, ylim_up)
    plt.show()
    return 0

def draw_authermers3D(var1):
    [hadron, charge, z_number, ph_number]=var
    X=np.arange(0.04, 0.24, 0.05)
    Qf=np.arange(1.1, 6.1, 1)

    x, qf=np.meshgrid(X,Qf)
    (h,l)=x.shape
    z=[0.239, 0.321, 0.424, 0.569]
    ph=[0.143, 0.298, 0.448, 0.770]
    par=par1
    aut=[]
    print(x)
    print('..........')
    print(qf)
    for i in range(h):
        aul=[]
        for j in range(l):
            var=[x[i][j],qf[i][j],z[z_number-1],ph[ph_number-1],'proton',hadron,charge]
            au=th.aut_th(var,par)
            aul.append(au)
        aut.append(aul)
        
    aut=np.array(aut)
    print(aut)
    title='A$_{UT}$ '+'z=%f '%z[z_number-1]+'ph=%f'%ph[ph_number-1]
    fig=plt.figure()
    fig.suptitle(title, fontsize=14)
    ax=Axes3D(fig)
    ax.plot_surface(x, qf, aut, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
    #cbar=plt.contourf(x, qf, aut, cmap=plt.get_cmap('rainbow'))
    #fig.colorbar(cbar)
    #plt.contourf(x, qf, aut, cmap=plt.get_cmap('rainbow'))
    
    #ax.set_zlim(-1, 1)
    ax.set_xlabel('x', fontsize=14)
    ax.set_ylabel('Q$^{2}$(GeV)$^{2}$', fontsize=14)
    plt.show()
    
    return 0
   
def draw_rad_per(var):
    mu=var
    
    x=np.linspace(0.001,20,50)
    #y1=[ad.b_x(i) for i in x]
    y2=[evl.g(mu, i) for i in x]
    
    #plt.plot(x,y1,label='bx')
    plt.plot(x,y2,label='y2')
    
    ax=plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_position(('data',0))
    ax.spines['left'].set_position(('data',0))
    
    
    #plt.ylim(0,2)
    plt.ylabel('evlution(Q=%s):sum to g0'%mu, fontsize=14)
    plt.xlabel('b(MeV$^{-1}$)', fontsize=14)
    #plt.legend(loc='best', fontsize=14)
    plt.show()
    return 0

def path(var):
    b=var
    muf=2.0
    ztf=muf**2.0
    
    
    bstar=b/math.sqrt(1.0+b**2.0/BNP**2.)
    mui=2.0*math.exp(-th.gammaE)/b
    zti=mui**2.0
    zti=4.0*math.exp(-2.0*gammaE+3.0/2.0)/b**2.0+10.
    
    path1=integrate.quad(ad.gam_F, mui, muf, args=ztf)[0]-ad.rad_per_resum([mui, b])*math.log(ztf/zti)
    #path2=integrate.quad(ad.gam_F, mui, muf, args=zti)[0]-ad.rad([muf, b])*math.log(ztf/zti)
    path3=ad.evo_fac([muf,b])**(0.5)
    #path4=ad.in_R(mui, zti, muf, b)
    path=[math.exp(path1),path3]
    path=[path1,math.log(path3)]
    
    return path

def draw_path():
    #pole=2*math.exp(-gammaE)/91
    #print(2*math.exp(-gammaE)/91)
    #print(path(pole))
    b=np.linspace(0.01,3.03,50)
    y1=[path(i)[0] for i in b]
    y2=[path(i)[1] for i in b]
    #y3=[path(i)[2] for i in b]
  
    plt.plot(b,y1,label='old_S', color='b')
    plt.plot(b,y2,label='new_R', color='g')
    plt.xlabel('b(GeV$^{-1}$)', fontsize=14)
    plt.ylabel('ln(evolution)')
    #plt.plot(b,y3,label='y3')
    #plt.plot(b,y,label='evolution')
    plt.legend(loc='best', fontsize=14)
    tit='evolution:Q=2'
    plt.title(tit)
    plt.show()
    return 0

def d0(mu):
    d0=(4.0*lha.a_s(mu)*math.log(mu/91.0)+3.0*lha.a_s(mu))/mu
    return d0
    
def draw_fourier_fsi(par, flavor):#par=[par1, parfit]
    kt_d=0.001
    kt_u=2.
    kt=np.linspace(kt_d, kt_u, 3)
    mu=2.
    x=.1
    target='proton'
    
    zfh=1.
    nam1=''
    nam2=''
    if flavor==1 :
        zfh=1.
        nam1=''
        nam2='d'
    elif flavor==2:
        zfh=-1.
        nam1='-'
        nam2='u'
    elif flavor==3:
        zfh=1.
        nam1=''
        nam2='s'
    else:
        zfh=-1.
        nam1='-'
        nam2='sea'
        
    y1=[]
    y2=[]
    #y3=[]
    #y3_std=[]
    for i in kt:
        var=[x, i, mu, target]
        var_random=[x, i, mu, target,flavor]
        four_fsi1=th.fourier_fsiver_all(var,par[0])[flavor]
        #four_fsi2=th.fourier_fsiver_all(var,par[1])[flavor]
        four_fsi2=th.fourier_f1_all(var)[flavor]
        #four_fsi3=th.fourier_mean_std(var_random)

        y1.append(four_fsi1)
        y2.append(four_fsi2)
        #y3.append(four_fsi3[0])
        #y3_std.append(four_fsi3[1])
        
    #print(y1)
    #print(y2)
    #print(y3)
    #print(y3_std)
    
    #y3up=np.array(y3)+np.array(y3_std)
    #y3down=np.array(y3)-np.array(y3_std)
    
    #print(y3up)
    #print(y3down)
    
    
    y1=np.log10(zfh*np.array(y1))
    y2=np.log10(np.array(y2))
    #y3=np.log10(zfh*np.array(y3))
    #y3up=np.log10(zfh*np.array(y3up))
    #y3down=np.log10(zfh*np.array(y3down))
    
    '''
    y1=zfh*np.array(y1)
    y2=zfh*np.array(y2)
    y3=zfh*np.array(y3)
    y3up=zfh*np.array(y3up)
    y3down=zfh*np.array(y3down)
    '''
    
    plt.figure()
    
    #x,y,轴刻度与说明设置
    plt.xlim((0, kt_u))
    plt.ylim((-4, 0))
    plt.xlabel('kt(GeV$^{-1}$)', fontsize=14)
    plt.ylabel(nam1+r'f$^{\perp}_{T;%s \leftarrow p }$' %nam2, fontsize=20)
    x_ticks=np.linspace(0, kt_u, 11)
    y_ticks=[ 0.0001,.001, .01, .1, 1., 10, 100 ]
    #y_ticks=[-0.5, 0, 0.5, 1., 1.5, 2., 2.5]
    ylog=[math.log10(i) for i in y_ticks]
    plt.xticks(x_ticks)
    #plt.yticks(y_ticks)
    plt.yticks(ylog,y_ticks)
    
    '''
    #坐标轴边框设置(颜色和位置)
    ax=plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_position(('data',-4))
    ax.spines['left'].set_position(('data',0))
    #ax.xaxis.set_ticks_position('bottom')
    #ax.yaxis.set_ticks_position('left')
    '''
    #图例设置
    ax=plt.gca()
    a=r'$\mu$=%s'%mu
    b=':x=%s'%x
    plt.title(a+b)
    plt.plot(kt,y1,label='par2021_fsivers' )
    plt.plot(kt,y2,label='par2021_f1')
    #plt.plot(kt,y3,label='par_random_mean')
    #ax.fill_between(kt, y3down, y3up, alpha=0.2)
    plt.legend(loc='best', fontsize=14)
    
    plt.show()
    return 0

def draw_fourier_fsi3D(var):
    [mu, flavor]=var
    
    ax = plt.figure().add_subplot(projection='3d')

    target='proton'
    verts = []

    kt_d=0.001
    kt_u=7.0

    n=4

    kt = np.linspace(kt_d, kt_u, 10)

    zs = range(n)

    for i in zs:
        
        x=0.005*2.**i
        z_mean=[]
        z_std=[]
        for j in kt:
            z_mean_std=th.fourier_mean_std([x, j, mu, target, flavor])
            print(j,z_mean_std)
            z_mean.append(z_mean_std[0])
            z_std.append(z_mean_std[1])
        
            
        xl=kt
        yl=[i*l/l for l in xl]
        zl=z_mean
        
        ax.plot(xl,  yl, zl, color='r')
        
        zme=np.array(z_mean)
        zst=np.array(z_std)
        
        zuu=zme+zst
        zdd=zme-zst
        
        zbaru=[zuu[i] for i in range(len(zuu))]
        zbard=[zdd[-(i+1)] for i in range(len(zdd))]
        zbar=zbaru+zbard
        
        ktbaru=[kt[i] for i in range(len(kt))]
        ktbard=[kt[-(i+1)] for i in range(len(kt))]
        ktbar=ktbaru+ktbard
        
            
        verts.append([*zip(ktbar, zbar)])

    poly = PolyCollection(verts, facecolors=['b', 'b', 'b', 'b', 'b', 'b'], alpha=.5)
    ax.add_collection3d(poly, zs=zs, zdir='y')

    ax.set_xlabel('kt')
    ax.set_ylabel('x')
    ax.set_zlabel('Down quark Sivers function Q=%s'%mu)
    ax.set_xlim(0, kt_u)
    ax.set_ylim(-0.2, n)
    ax.set_zlim(0, .5)

    x_ticks=np.linspace(0, kt_u, 6)
    y_ticks=[.005*2.**i for i in zs]
    ynum=[i for i in zs]
    plt.xticks(x_ticks)
    plt.yticks(ynum,y_ticks)

    plt.show()

def draw_all_par_fourier(flavor):
    kt_d=0.01
    kt_u=2.
    kt=np.linspace(kt_d, kt_u, 10)
    mu=2.
    x=.16
    target='proton'
    
    zfh=1.
    nam1=''
    nam2=''
    if flavor==1 :
        zfh=1.
        nam1=''
        nam2='d'
    elif flavor==2:
        zfh=-1.
        nam1='-'
        nam2='u'
    elif flavor==3:
        zfh=1.
        nam1=''
        nam2='s'
    else:
        zfh=-1.
        nam1='-'
        nam2='sea'
        
    y=[]

    for par in ex.data_par():
        yb=[]
        for i in kt:
            var=[x, i, mu, target]
            var_random=[x, i, mu, target,flavor]
            four_fsi1=th.fourier_fsiver_all(var,par)[flavor]
            yb.append(four_fsi1)
        yb=zfh*np.array(yb)
        y.append(yb)
 
    
    plt.figure()
    

    plt.xlim((0, kt_u))
    plt.ylim((-0.5, 4))
    plt.xlabel('kt(GeV$^{-1}$)', fontsize=14)
    plt.ylabel(nam1+r'f$^{\perp}_{T;q\to %s }$' %nam2, fontsize=14)
    x_ticks=np.linspace(0, kt_u, 11)
 
    plt.xticks(x_ticks)
  
    #图例设置
    plt.title('NLL:random_par x=0.16')
    for i in y:
        plt.plot(kt, i,label='random_par')
    #plt.legend(loc='best', fontsize=14)
    plt.show()
    return 0
'''
a=[1,2,3]
b=[4,5,6]
aa=np.array(a)
bb=np.array(b)

bbb=[b[-(i+1)] for i in range(len(b))]

print(a+b)
'''

def draw_c1_inte():
    x=0.2
    b=1.0
    y=np.linspace(x, 1 ,50)
   
    c1=[th.c1_p_qg_inte(i, x, b,'proton') for i in y]
  
    plt.plot(y, c1,label='c1')
    
    ax=plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_position(('data',0))
    ax.spines['left'].set_position(('data',0))
    
    #plt.ylim(0,2)
    plt.ylabel('c1', fontsize=14)
    plt.xlabel('y', fontsize=14)
    plt.show()
    return 0

def draw__fsi_f1(ax, kt_u, kt, y1, y2, title, flavor):

    nam1=''
    nam2=''
    if flavor==1 :
        nam1=''
        nam2='d'
        y_down=-3.3
        y_up=1
    elif flavor==2:
        nam1='-'
        nam2='u'
        y_down=-3.3
        y_up=2
    elif flavor==3:
        nam1=''
        nam2='s'
        y_down=-4.3
        y_up=0
    elif flavor==-1 :
        nam1='-'
        nam2='db'
        y_down=-4.3
        y_up=1
    elif flavor==-2:
        nam1='-'
        nam2='ub'
        y_down=-4.3
        y_up=1
    elif flavor==-3:
        nam1='-'
        nam2='sb'
        y_down=-4.3
        y_up=1

    ax.set_xlim((0, kt_u))
    ax.set_ylim((y_down, y_up))
    ax.set_xlabel('kt(GeV$^{-1}$)', fontsize=10)
    ax.set_ylabel(nam1+r'f$^{\perp}_{T;%s \leftarrow p }$' %nam2, fontsize=14)
    x_ticks=np.linspace(0, kt_u, 11)
    #y_ticks=[ 0.0001,.001, .01, .1, 1., 10, 100 ]
    #ylog=[math.log10(i) for i in y_ticks]
    ax.set_xticks(x_ticks)
    #ax.set_yticks(ylog,y_ticks)

    plt.plot(kt,y1,label='par2021_fsivers' )
    plt.plot(kt,y2,label='par2021_f1')
    ax.legend(loc='best', fontsize=10)

    return 0

def draw_all_fsi_f1():
    fu=th.flavor_use
    par=par1
    kt_d=0.001
    kt_u=2.
    kt=np.linspace(kt_d, kt_u, 3)
    mu=2.
    x=.1
    target='proton'

    zfh=1.
 
    title=[mu, x]
    y1=[]
    y2=[]

    for i in kt:
        var=[x, i, mu, target]
        
        four_fsi=th.fourier_fsiver_all(var,par)
        four_f1=th.fourier_f1_all(var)

        y1.append(four_fsi)
        y2.append(four_f1)
        
    
    fsi={}
    f1={}
    for i in fu:
        if i==1 :
            zfh=1.
        elif i==2:
            zfh=-1.
        elif i==3:
            zfh=1.
        else:
            zfh=-1.
        fsi[i]=[np.log10(zfh*y1[j][i]) for j in range(len(kt))]
        f1[i]=[np.log10(y2[j][i]) for j in range(len(kt))]

    #print(fsi[1])
    #print(f1[1])
    a=r'$\mu$=%s'%mu
    b=':x=%s'%x
    fig=plt.figure(num=3)
    fig.suptitle(a+b)
    
    for i in range(6):
        ax=plt.subplot(3,2,i+1)
        flavor=fu[i]
        pos=ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width, pos.height*0.8])
        draw__fsi_f1(ax, kt_u, kt, fsi[flavor], f1[flavor], title, flavor)

    plt.show()
    return 0

draw_all_fsi_f1()








'''
draw_AUT('JLab2011', par, 'x')
draw_AUT('JLab2014', par, 'x')
draw_AUT('compass2009', par, 'ph')


hadron='pi'
charge='+'
draw_hermers3D([hadron,charge,1,1],par)
draw_hermers3D([hadron,charge,2,1],par)
draw_hermers3D([hadron,charge,3,1],par)
draw_hermers3D([hadron,charge,4,1],par)
'''



