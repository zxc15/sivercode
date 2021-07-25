import lhapdf
import numpy as np



pdf=lhapdf.mkPDF("CT14lo",0)
aphas=lhapdf.mkPDF("CT14lo",0)
ff_pion=lhapdf.mkPDF("JAM20-SIDIS_FF_pion_nlo",0)
ff_kaon=lhapdf.mkPDF("JAM20-SIDIS_FF_kaon_nlo",0)
ff_hadron=lhapdf.mkPDF("JAM20-SIDIS_FF_hadron_nlo",0)

def a_s(var):
    mu=var
    a_s=aphas.alphasQ(mu)
    return a_s

def f_1(var):
	x, mu, target=var[0], var[1], var[2]
	if   target=='proton':
		u, d=pdf.xfxQ(2, x, mu)/x, pdf.xfxQ(1, x, mu)/x
		ub, db=pdf.xfxQ(-2, x, mu)/x, pdf.xfxQ(-1, x, mu)/x
	elif target=='neutron':
		u, d=pdf.xfxQ(1, x, mu), pdf.xfxQ(2, x, mu)
		ub, db=pdf.xfxQ(-1, x, mu), pdf.xfxQ(-2, x, mu)
	elif target=='deuteron':
		u=(pdf.xfxQ(2, x, mu)/x + pdf.xfxQ(1, x, mu)/x)/2.0
		d=(pdf.xfxQ(1, x, mu)/x + pdf.xfxQ(2, x, mu)/x)/2.0
		ub=(pdf.xfxQ(-2, x, mu)/x + pdf.xfxQ(-1, x, mu)/x)/2.0
		db=(pdf.xfxQ(-1, x, mu)/x + pdf.xfxQ(-2, x, mu)/x)/2.0
	elif target=='3He':
		u=(2.0*pdf.xfxQ(2, x, mu)/x + pdf.xfxQ(1, x, mu))/3.0
		d=(2.0*pdf.xfxQ(1, x, mu)/x + pdf.xfxQ(2, x, mu))/3.0
		ub=(2.0*pdf.xfxQ(-2, x, mu)/x + pdf.xfxQ(-1, x, mu))/3.0
		db=(2.0*pdf.xfxQ(-1, x, mu)/x + pdf.xfxQ(-2, x, mu))/3.0
	else :
		print('Fail to match target!')
		return 0
		
	s, sb=pdf.xfxQ(3, x, mu)/x, pdf.xfxQ(-3, x, mu)/x
	c, cb=pdf.xfxQ(4, x, mu)/x, pdf.xfxQ(-4, x, mu)/x
	b, bb=pdf.xfxQ(5, x, mu)/x, pdf.xfxQ(-5, x, mu)/x
	
	f_1={2:u, 1:d, 3:s, 4:c, 5:b, -2:ub, -1:db, -3:sb, -4:cb, -5:bb}
	  
	return f_1

def d_1(var):
	z, mu, hadron, charge=var[0], var[1] ,var[2], var[3]

	if hadron=='pi' and charge=='+':
		u, d, s, c, b=ff_pion.xfxQ(2,z,mu)/z, ff_pion.xfxQ(1,z,mu)/z, ff_pion.xfxQ(3,z,mu)/z, ff_pion.xfxQ(4,z,mu)/z, ff_pion.xfxQ(5,z,mu)/z
		ub, db, sb, cb, bb=ff_pion.xfxQ(-2,z,mu)/z, ff_pion.xfxQ(-1,z,mu)/z, ff_pion.xfxQ(-3,z,mu)/z, ff_pion.xfxQ(-4,z,mu)/z, ff_pion.xfxQ(-5,z,mu)/z
	elif hadron=='pi' and charge=='-':
		ub, db, sb, cb, bb=ff_pion.xfxQ(2,z,mu)/z, ff_pion.xfxQ(1,z,mu)/z, ff_pion.xfxQ(3,z,mu)/z, ff_pion.xfxQ(4,z,mu)/z, ff_pion.xfxQ(5,z,mu)/z
		u, d, s, c, b=ff_pion.xfxQ(-2,z,mu)/z, ff_pion.xfxQ(-1,z,mu)/z, ff_pion.xfxQ(-3,z,mu)/z, ff_pion.xfxQ(-4,z,mu)/z, ff_pion.xfxQ(-5,z,mu)/z
	elif hadron=='K' and charge=='+':
		u, d, s, c, b=ff_kaon.xfxQ(2,z,mu)/z, ff_kaon.xfxQ(1,z,mu)/z, ff_kaon.xfxQ(3,z,mu)/z, ff_kaon.xfxQ(4,z,mu)/z, ff_kaon.xfxQ(5,z,mu)/z
		ub, db, sb, cb, bb=ff_kaon.xfxQ(-2,z,mu)/z, ff_kaon.xfxQ(-1,z,mu)/z, ff_kaon.xfxQ(-3,z,mu)/z, ff_kaon.xfxQ(-4,z,mu)/z, ff_kaon.xfxQ(-5,z,mu)/z
	elif hadron=='K' and charge=='-':
		ub, db, sb, cb, bb=ff_kaon.xfxQ(2,z,mu)/z, ff_kaon.xfxQ(1,z,mu)/z, ff_kaon.xfxQ(3,z,mu)/z, ff_kaon.xfxQ(4,z,mu)/z, ff_kaon.xfxQ(5,z,mu)/z
		u, d, s, c, b=ff_kaon.xfxQ(-2,z,mu)/z, ff_kaon.xfxQ(-1,z,mu)/z, ff_kaon.xfxQ(-3,z,mu)/z, ff_kaon.xfxQ(-4,z,mu)/z, ff_kaon.xfxQ(-5,z,mu)/z
	elif hadron=='p' and charge=='+':
		u, d, s, c, b=ff_hadron.xfxQ(2,z,mu)/z, ff_hadron.xfxQ(1,z,mu)/z, ff_hadron.xfxQ(3,z,mu)/z, ff_hadron.xfxQ(4,z,mu)/z, ff_hadron.xfxQ(5,z,mu)/z
		ub, db, sb, cb, bb=ff_hadron.xfxQ(-2,z,mu)/z, ff_hadron.xfxQ(-1,z,mu)/z, ff_hadron.xfxQ(-3,z,mu)/z, ff_hadron.xfxQ(-4,z,mu)/z, ff_hadron.xfxQ(-5,z,mu)/z
	elif hadron=='p' and charge=='-':
		ub, db, sb, cb, bb=ff_hadron.xfxQ(2,z,mu)/z, ff_hadron.xfxQ(1,z,mu)/z, ff_hadron.xfxQ(3,z,mu)/z, ff_hadron.xfxQ(4,z,mu)/z, ff_hadron.xfxQ(5,z,mu)/z
		u, d, s, c, b=ff_hadron.xfxQ(-2,z,mu)/z, ff_hadron.xfxQ(-1,z,mu)/z, ff_hadron.xfxQ(-3,z,mu)/z, ff_hadron.xfxQ(-4,z,mu)/z, ff_hadron.xfxQ(-5,z,mu)/z
	else :
		print('error: Fail to match '+var[2]+var[3])
		return 0
	d_1={2:u, 1:d, 3:s, 4:c, 5:b, -2:ub, -1:db, -3:sb, -4:cb, -5:bb}
	return d_1
 
 

