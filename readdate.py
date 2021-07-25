
import pandas as pd
import os

def aut_ex(var):#var=[experiment,hadron,charge]

	skip=[0,1,2,3,4]
	title=['x','Q^2','y','z','ph','w','clo','c1','c2','c3','siv','s1','s2','s3']
	path_head='/home/zengchunhua/workarea/python/siverdata/'
	exp=['compass2009','compass2015','JLab2011','JLab2014','hermers']
#	had=['K','pi','p']
#	cha=['+','-','0']
	
	if var[0]==exp[4]:
		path_head=path_head+'hermers/1D/'
		fil=path_head+var[1]+var[2]+'/x.txt'
		
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
		
		fil=path_head+expt+exp[k]+name1+var[1]+var[2]+'_x.txt'
	
	
	
	read_data=pd.read_csv(fil, skiprows=skip, header=None, names=title,sep='\t')
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
		siv_data.append(data_set)
	return siv_data



