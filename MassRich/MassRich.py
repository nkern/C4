'''
Analysis of Mass--Richness relationship for Millennium Simulation Henriques Light Cone
and SDSS DR12 Spectra with C4-identified Clusters
Below analysis relies on richnesses, phase spaces etc. data from Chris found in
	/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH*	or
	/nfs/christoq_ls/C4/sdssdr12
'''

## Import Modules
import numpy as np
import matplotlib.pyplot as mp
import astropy.io.fits as fits
import astropy.stats as astats
from sklearn import linear_model
from scipy import stats
from scipy.optimize import leastsq
import os
import sys
from fits_table import fits_table
from c4_richness import RICHNESS
import astropy.cosmology as cosmo
import DictEZ as ez
from stack_class_c4 import CFOUR
from causticpy import Caustic



## Constants
c = 2.99792458e5		#km/s
H0 = 70


## Flags



## Functions ##
def calc_scat(self,mass,observable,loglog=True,binobs=True,step=0.05,data_cut=None):
	''' Does Mass Mixing Statistics '''
	## First Calculate Scatter of Richness Relationship
	# Get Data Cut
	cut = np.where(observable > 0)
	if obs_cut == None:
		print ' Find observable value above which vertical scatter is fully encompassed'
		mp.plot(np.log(observable[cut]),np.log(true_mass[cut]),'ko',alpha=.75)
		mp.xlabel('Observable',fontsize=16)
		mp.ylabel('M200',fontsize=16)
		mp.show()
		obs_cut = input('What is that value? (Enter as a float, ex. 1.5) : ')

	# Correct obs_cut
	if type(obs_cut) == int or type(obs_cut) == float: obs_cut = np.array([obs_cut])

	# Now Calculate Scatter
	SCAT = []
	for i in range(len(obs_cut)):
		data_cut = np.where(np.log(observable) > obs_cut[i])[0]
		print 'Number of Data Points:',len(data_cut)
		step = step
		obs_max = np.log(observable[data_cut]).max()
		bins = np.arange(obs_cut[i],obs_max,step)

		mass_bins = np.array( map(lambda x: (true_mass[data_cut])[np.where((np.log(observable[data_cut])>x)&(np.log(observable[data_cut])<(x+step)))], bins))
		mass_avgs = np.array(map(lambda x:np.median(x),mass_bins))
		mass_fracs = []
		for j in range(len(mass_avgs)):
			mass_fracs.extend(np.log(mass_bins[j] / mass_avgs[j]))
		mass_fracs = np.array(mass_fracs)

		scat = astats.biweight_midvariance(mass_fracs)
		SCAT.append(scat)

	return np.array(SCAT)

def residual(vars, x, data):
	""" model for data, subtract data """
	a = vars[0]
	b = vars[1]

	model = a*x + b

	return (data-model)**2

def gauss_resid(vars,x,data):
	a = vars[0]
	b = vars[1]
	model = (1/a*np.sqrt(2*np.pi)) * np.exp(-(x-b)/(2*a**2))
	return (data-model)**2

def stats(X,Y,X_low,X_high,X_step,fitline=False,log10=False):
	''' Get Scatter for X Y relationship where X and Y are in linear space'''

	# Put data into log space
	if log10 == True: X = np.log10(X); Y = np.log10(Y)
	else: X = np.log(X); Y = np.log(Y)

	# Make Bins
	bins = np.arange(X_low,X_high-0.001,X_step)
	bin_cent = (bins + np.arange(X_low+X_step,X_high-0.001+X_step,X_step))/2.
	cut = np.where(X>X_low)[0]
	Xbins = np.array(map(lambda x: X[np.where((X<x+X_step)&(X>=x))],bins))
	Ybins = np.array(map(lambda x: Y[np.where((X<x+X_step)&(X>=x))],bins))

	# Take average of each bin
	Yavgs = np.array(map(np.median,Ybins))

	# Get deviation of each point from bin avg
	if fitline == True:
		model = linear_model.LinearRegression()
		model.fit(X[cut].reshape(X[cut].size,1),Y[cut])
		Ybins_del = Ybins - (Xbins*model.coef_[0] + model.intercept_)
	else:
		Ybins_del = Ybins-Yavgs

	# Get Scatter
	try:
		Ybins_scat = np.array(map(astats.biweight_midvariance,Ybins_del))
	except:
		Ybins_scat = np.array(map(np.std,Ybins_del))

	Ybins_scat_err = np.array(map(lambda x: 1/np.sqrt(x.size),Ybins))

	if fitline == True:
		return dict(X=X,Y=Y,bins=bins,bin_cent=bin_cent,cut=cut,Xbins=Xbins,Ybins=Ybins,Yavgs=Yavgs,Ybins_del=Ybins_del,Ybins_scat=Ybins_scat,Ybins_scat_err=Ybins_scat_err,model=model)
	else:
		return dict(X=X,Y=Y,bins=bins,bin_cent=bin_cent,cut=cut,Xbins=Xbins,Ybins=Ybins,Yavgs=Yavgs,Ybins_del=Ybins_del,Ybins_scat=Ybins_scat,Ybins_scat_err=Ybins_scat_err)

## Truth Table (M-R for Henriques using 3D members and True Mass) ##
panel1 = False
if panel1 == True:
	# Constants
	data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH'
	root = '/nfs/christoq_ls/nkern/C4/MassRich/TRUTH'

	# Load Data
	truth = fits.open(data_root+'/c4_cluster_truth_revH100_rev5.fits')[1].data
	truth['m200mean']*=1e10
	truth['m200crit']*=1e10

	Rich = truth['n10virm19']
	Mass = truth['m200crit']

	rich_low,rich_high,rich_step = 3.3,4.8,0.3
	richbins = np.arange(rich_low,rich_high,rich_step)
	mass_low,mass_high,mass_step = 14.25,15.15,0.2
	massbins = np.arange(mass_low,mass_high,mass_step)

	d = stats(Rich,Mass,rich_low,rich_high,rich_step,log10=False,fitline=True)
	globals().update(d)

	plot = True
	if plot == True:
		fig,ax = mp.subplots(1,2,figsize=(13,7))
		## ax[0] is MassRich
		p1, = ax[0].plot(X,Y,'ko',alpha=.2)
		ax[0].set_xlim(X.min()-.1,X.max()+.1)
		ax[0].set_ylim(Y.min()-.1,Y.max()+.5)
		ax[0].set_xlabel('log( N200 )',fontsize=15)
		ax[0].set_ylabel('log( Mass )',fontsize=15)
		ax[0].fill_between(richbins,Y.min()-.1,Y.max()+.5,alpha=.1)
		ax[0].plot(X,X*model.coef_+model.intercept_,'r')

		## ax[1] is Scatter
		p2 = ax[1].errorbar(bin_cent,Ybins_scat,yerr=Ybins_scat*Ybins_scat_err,fmt='o',color='b',markersize=8)	
		ax[1].set_xlim(bin_cent.min()-.2,bin_cent.max()+.2)
		ax[1].set_ylim(Ybins_scat.min()-.025,Ybins_scat.max()+.05)
		ax[1].set_xlabel('log( N200 )',fontsize=15)
		ax[1].set_ylabel('absolute scatter',fontsize=15)

		# Fit scatter points to polynomial
		def residual(vars, x, data, err):
			""" model for data, subtract data """
			a = vars[0]
			b = vars[1]
			c = vars[2]
			d = vars[3]
			model = a*(x-b)**c + d
			return ((data-model)/err)**2

		xmodel = np.arange(bin_cent[0]-.1,bin_cent[-1]+.3,0.01)
		vars = [0.125,0.1,-1.5,0.01]
		#vars = [0.015,4,2,0]
		out = leastsq(residual,vars,args=(bin_cent,Ybins_scat,Ybins_scat_err))
	#	p3, = ax[1].plot(xmodel,out[0][0]*(xmodel-out[0][1])**(out[0][2])+out[0][3],'g')
	#	ax[1].legend([p3],[str(np.round(out[0][0],3))+'(x-'+str(np.round(out[0][1],3))+')^('+str(np.round(out[0][2],3))+') + '+str(np.round(out[0][3],3))],fontsize=15)

		savefig = True
		if savefig == True:
			fig.savefig('TRUTH_scatter.eps',bbox_inches='tight')

## Truth Caustic Table (M-R for Henriques using Projected Members and Caustic Mass) ##
panel2 =False
if panel2 == True:
	# Constants
	data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC'
	root = '/nfs/christoq_ls/nkern/C4/MassRich/TRUTH_CAUSTIC'
        H0 = 70
        h = H0 / 100.

	# Load Data
	halos = fits.open(data_root+'/c4_cluster_single_whalo.fits')[1].data
#	halo_avg_data = fits.open(root+'/halo_avg_data.fits')[1].data

	# Table quoted values from Henriques
	M200 = halos['M200CRIT']
	R200 = halos['R200CRIT']
	RA = halos['HALO_RA']
	DEC = halos['HALO_DEC']
	Z = halos['HALO_Z']
	HVD = halos['VELDISP']

	# Richnesses
        new = fits.open(data_root+'/richness_cr200_bcg/new.fits')[0].data
        newb = fits.open(data_root+'/richness_cr200_bcg/newb.fits')[0].data

	# Calculate N200 and Nspec
	from causticpy import *
	C = Caustic()

        ### Run with Halo values from either Henriques or C4 algorithm
	henriq = True
	if henriq == True:
		# Halo data from Henriques

		# Good Halos
		good = np.where((halos['halo_z']>0.03))[0]

		RA = halos['halo_ra']
		DEC = halos['halo_dec']
		Z = halos['halo_z']
		RVIR = halos['r200crit']

	else:
		# Halo data from C4

		# Good Halos
		good = np.where(halos['z_biwt']>0.03)[0]

		RA = halo_avg_data['ra_avg']
		DEC = halo_avg_data['dec_avg']
		Z = halo_avg_data['z_avg']
		RVIR = halos['rho_rv1500']*0.7**2
		RVIR[np.where(RVIR<0.5)]=0.5

	raise NameError

	Nspec = []
	N200 = []
	HVD = []

	good2 = []
	for i in good:
		try:
			gal_ra, gal_dec, gal_umag, gal_gmag, gal_rmag, gal_imag, gal_zmag, gal_z = np.loadtxt(data_root+'/proj_data/C4_dr12_'+str(i)+'_data.csv',dtype='float',delimiter=',',skiprows=1,usecols=(1,2,3,4,5,6,7,13),unpack=True)
		except: continue

		clus_ra = RA[i]
		clus_dec = DEC[i]
		clus_z = Z[i]
		rvir = RVIR[i]

		ang_d,lum_d = C.zdistance(clus_z,H0)
		angles = C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
		rdata = angles * ang_d
		vdata = c * (gal_z - clus_z) / (1 + clus_z)

		# Calculate HVD
		rlim = rvir * 1.25   # Mpc
		rad_cut = np.where(rdata < rlim)[0]
		if len(rad_cut) == 0:
			continue
		elif len(rad_cut) > 2:
			hvd = astats.biweight_midvariance(vdata[np.where((rdata < rlim)&(np.abs(vdata)<3000))])
		else:
			hvd = np.std(vdata[np.where(rdata < rlim)])

		Nspec.append( np.where((rdata < rvir) & (np.abs(vdata) < hvd * 1.5))[0].size )
		HVD.append(hvd)
		good2.append(i)
		if i % 200 == 0: print i

	Nspec,HVD,good2 = np.array(Nspec),np.array(HVD),np.array(good2)
	N200 = []
	for i in good2:
		N200.append(new[i][2,1,4,9] - newb[i][2,1,4,1])
	N200 = np.array(N200)
	new_halos = halos[good2]

        write_new = True
        if write_new == True:
                keynames = new_halos.names
                keynames.extend(['N200','Nspec','HVD','orig_order'])
                d = {}
                for k in range(len(new_halos.dtype)):
                        d[keynames[k]] = new_halos[keynames[k]]
                d['N200'] = N200
                d['Nspec'] = Nspec
                d['HVD'] = HVD
                d['orig_order'] = good2
                fits_table(d,keynames,root+'/halos.fits')


## Truth Caustic C4 Table (M-R for C4 Catalogue in Henriques using Projected Members and Caustic Masses) ##
panel3 = False
if panel3 == True:
	# Constants
	data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC_C4'
	root = '/nfs/christoq_ls/nkern/C4/MassRich/TRUTH_CAUSTIC_C4'
        H0 = 70
        h = H0 / 100.

	# Load Data
	halos = fits.open(data_root+'/c4_cluster_single_whalo.fit')[1].data

	M200_true = halos['M200CRIT']
	R200_true = halos['R200CRIT']

## SDSS C4 (M-R for C4 Catalogue in Real Universe using SDSS DR12 spectra / photometry and Caustic Masses) ##
panel4 = False
if panel4 == True:
	# Constants
	data_root = '/nfs/christoq_ls/C4/sdssdr12'
	root = '/nfs/christoq_ls/nkern/C4/MassRich/DR12'
	H0 = 70
	h = H0 / 100.

	# Load Data
	halos = fits.open(data_root+'/c4_cluster_single.fits')[1].data
	new = fits.open(data_root+'/richness_cr200/new.fits')[0].data
	newb = fits.open(data_root+'/richness_cr200/newb.fits')[0].data

	# Select Out Good Data
	good = np.where((halos['z_biwt']>0.003)&(halos['single']==9999))[0]

	#raise NameError
	# Calculate N200 and Nspec
	from causticpy import *
	C = Caustic()

	RA = halos['ra_bcg']
	DEC = halos['dec_bcg']
	Z = halos['z_biwt']
	RVIR = halos['RHO200_RV1500']

	Nspec = []
	HVD = []
	good2 = []
	for i in good:
		try: 
			gal_ra, gal_dec, gal_umag, gal_gmag, gal_rmag, gal_imag, gal_zmag, gal_z = np.loadtxt(data_root+'/proj_data/C4_dr12_'+str(i)+'_data.csv',dtype='float',delimiter=',',skiprows=1,usecols=(1,2,3,4,5,6,7,13),unpack=True)
		except: continue

		ang_d,lum_d = C.zdistance(Z[i],H0)
		angles = C.findangle(gal_ra,gal_dec,RA[i],DEC[i])
		rdata = angles * ang_d
		vdata = c * (gal_z - Z[i]) / (1 + Z[i])

		# Calculate HVD
		rlim = RVIR[i] * 1.25	# Mpc
		rad_cut = np.where(rdata < rlim)[0]
		if len(rad_cut) == 0:
			continue
		elif len(rad_cut) > 2:
			hvd = astats.biweight_midvariance(vdata[np.where((rdata < rlim)&(np.abs(vdata)<3000))])
		else:
			hvd = np.std(vdata[np.where(rdata < rlim)])

		Nspec.append( np.where((rdata < rlim/1.25) & (np.abs(vdata) < hvd * 1.5))[0].size )
		HVD.append(hvd)
		good2.append(i)
		if i % 200 == 0: print i	

	Nspec,HVD,good2 = np.array(Nspec),np.array(HVD),np.array(good2)

	N200 = []
	for i in good2:
		N200.append(new[i][6,1,9,9] - newb[i][6,1,9,1])

	N200 = np.array(N200)

	new_halos = halos[good2]

	write_new = False
	if write_new == True:
		keynames = new_halos.names
		keynames.extend(['N200','Nspec','HVD','orig_order'])
		d = {}
		for k in range(len(new_halos.dtype)):
			d[keynames[k]] = new_halos[keynames[k]]
		d['N200'] = N200
		d['Nspec'] = Nspec
		d['HVD'] = HVD
		d['orig_order'] = good2
		fits_table(d,keynames,root+'/halos.fits')






## Sub Sampled Mass Bias Trend ##
'''
Run the caustic technique over individual clusters with high sampling (also stack them)
Run the caustic tehcnique over individual clusters with low sampling (also stack them)
We expect to see little change in the stacked results, large change in individual results

	At first try doing this over panel3 (truth_caustic_c4)
'''
subsample = False
if subsample == True:
	# Set Constants and Flags
	root = '/nfs/christoq_ls/'
	data_loc = 'MILLENNIUM/Henriques/TRUTH_CAUSTIC_C4/'

	def fits_data(fits_data):
		d = {}
		for i in range(len(fits_data.dtype)):
			d[fits_data.columns[i].name] = fits_data[fits_data.columns[i].name][np.where(fits_data[fits_data.columns[i].name]!=0)]
		return d

	data = []
	data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal100_Nclus20_Nens20.fits')[1].data))
        data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal75_Nclus20_Nens20.fits')[1].data))
        data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal50_Nclus20_Nens20.fits')[1].data))
        data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal25_Nclus20_Nens20.fits')[1].data))
        data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal15_Nclus20_Nens20.fits')[1].data))
	data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal10_Nclus20_Nens20.fits')[1].data))

	def residual(vars, x, data):
		a = vars[0]
		b = vars[1]
		model = a*x + b
		return (model-data)**2

	def residual2(vars, x, data):
		b = vars[0]
		model = 0.832 * x + b
		return (model-data)**2

	vars = [1,14]
	out = leastsq(residual, vars, args=(xdata,ydata))[0]

	plot = False
	if plot == True:
		fig,ax = mp.subplots(1,2,figsize=(14,7))
		p0, = ax[0].plot(xdata,np.log10(data[0]['EDGEMASS_EST']),'mo',markersize=9)
		p1, = ax[0].plot(xdata,np.log10(data[1]['EDGEMASS_EST']),'bo',markersize=9)
		p2, = ax[0].plot(xdata,np.log10(data[2]['EDGEMASS_EST']),'go',markersize=9)
		p3, = ax[0].plot(xdata,np.log10(data[3]['EDGEMASS_EST']),'co',markersize=9)
		p4, = ax[0].plot(xdata,np.log10(data[4]['EDGEMASS_EST']),'yo',markersize=9)
		p5, = ax[0].plot(xdata,np.log10(data[5]['EDGEMASS_EST']),'ro',markersize=9)
		ax[0].legend([p0,p1,p2,p3,p4,p5],['Ngal=100','Ngal=75','Ngal=50','Ngal=25','Ngal=15','Ngal=10'],loc=2)

		out = leastsq(residual, vars, args=(np.log10(data[0]['BINN200']),np.log10(data[0]['EDGEMASS_EST'])))[0]
		out1 = leastsq(residual2, [13], args=(xdata,np.log10(data[1]['EDGEMASS_EST'])))[0]
		out2 = leastsq(residual2, [13], args=(xdata,np.log10(data[2]['EDGEMASS_EST'])))[0]
		out3 = leastsq(residual2, [13], args=(xdata,np.log10(data[3]['EDGEMASS_EST'])))[0]
		out4 = leastsq(residual2, [13], args=(xdata,np.log10(data[4]['EDGEMASS_EST'])))[0]
		out5 = leastsq(residual2, [13], args=(xdata,np.log10(data[5]['EDGEMASS_EST'])))[0]
		out_params = [out[1],out1[0],out2[0],out3[0],out4[0],out5[0]]
		x = np.arange(1.4,2.6,.1)
		xdata = np.log10(data[0]['BINN200'])
	
		p0, = ax[1].plot(x,out[0]*x + out[1],'m')
		p1, = ax[1].plot(x,out[0]*x + out1[0],'b')
		p2, = ax[1].plot(x,out[0]*x + out2[0],'g')
		p3, = ax[1].plot(x,out[0]*x + out3[0],'c')
		p4, = ax[1].plot(x,out[0]*x + out4[0],'y')
		p5, = ax[1].plot(x,out[0]*x + out5[0],'r')

		bias = [np.median(np.log( 10**(out[0]*x+out[1]) / 10**(out[0]*x+i) )) for i in out_params]

		ax[1].legend([p0,p1,p2,p3,p4,p5],['-'+str(np.round(bias[0]*100,1))+'% bias','-'+str(np.round(bias[1]*100,1))+'% bias','-'+str(np.round(bias[2]*100,1))+'% bias','-'+str(np.round(bias[3]*100,1))+'% bias','-'+str(np.round(bias[4]*100,1))+'% bias','-'+str(np.round(bias[5]*100,1))+'% bias'],loc=2)
		
		ax[0].set_xlim(1,3)
		ax[0].set_ylim(13,16)
		ax[1].set_xlim(1,3)
		ax[1].set_ylim(13,16)
		ax[0].set_xlabel('log10 N200',fontsize=15)
		ax[0].set_ylabel('log10 M200',fontsize=15)

		savefig = False
		if savefig == True: fig.savefig('gal_subsample.pdf',bbox_inches='tight')
	

## Calculate Nspec and N200 from Chris' data ##
calc_nspec = False
if calc_nspec == True:

	data_loc = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC_C4'

	halos = fits.open(data_loc+'/c4_cluster_single_whalo.fit')[1].data

	good = np.where((halos['z_biwt']>0.003)&(halos['haloid']!=0)&(halos['single']==9999))[0]

	new = fits.open(data_loc+'/old/richness_cr200/new.fits')[0].data
	newb = fits.open(data_loc+'/old/richness_cr200/newb.fits')[0].data

	N200 = []
	for i in range(5020):
		N200.append(new[i][6,1,9,9] - newb[i][6,1,9,1])

	N200 = np.array(N200)[good]

	Nspec = []
	new_good = []
	for i in range(len(good)):
		try:
			f = open(data_loc+'/caustics_geom_q10/caustic_C4_'+str(i)+'.results')
			new_good.append(i)
		except:
			continue
		for j in range(14):
			f.readline()
			if j == 13: b = f.readline()

		Nspec.append(int(b[-4:]))

	new_good = np.array(new_good)
	Nspec = np.array(Nspec)
	N200 = N200[new_good]


## Calculate Average Halo Positions from Galaxies
calc_avg = False
if calc_avg == True:
	from stack_class_c4 import CFOUR
	from causticpy import Caustic
	
	root = '/nfs/christoq_ls/nkern'
	data_loc = 'MassRich/TRUTH_CAUSTIC'
	write_loc = 'individual'

	C4 = CFOUR({'H0':70,'chris_data_root':'/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC'})
	C = Caustic()

	# Load Halos
        halos = fits.open(root+'/C4/'+data_loc+'/halos.fits')[1].data
        HaloID = halos['orig_order']
        RA = halos['halo_ra']
        DEC = halos['halo_dec']
        Z = halos['halo_z']
        Nspec = halos['Nspec']
        N200 = halos['N200']
        HVD = halos['HVD']
        RVIR = halos['RVIR']
        SINGLE = halos['single']
        SUB = halos['sub']
        NC4 = halos['nc4']

	RA_AVG,DEC_AVG,Z_AVG = [],[],[]
	# Loop Over Halos
	for i in range(len(halos)):
		if i % 100 == 0: print i

		# Assign Halo Properties
		clus_ra = RA[i]
		clus_dec = DEC[i]
		clus_z = Z[i]

		# Load Galaxies
		galdata = C4.load_chris_gals(HaloID[i])	
		gal_ra,gal_dec,gal_z,gal_gmags,gal_rmags,gal_imags = galdata

		# Project Galaxies
		ang_d,lum_d = C.zdistance(clus_z,H0)
		angles = C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
		rdata = angles * ang_d
		vdata = c * (gal_z - clus_z) / (1 + clus_z)

		# Take Average Three times
		cut1 = np.where((np.abs(vdata)<1000)&(rdata<1.5))[0]
		clus_ra = astats.biweight_location(gal_ra[cut1])
		clus_dec = astats.biweight_location(gal_dec[cut1])
		clus_z = astats.biweight_location(gal_z[cut1])

		ang_d,lum_d = C.zdistance(clus_z,H0)
		angles = C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
		rdata = angles * ang_d
		vdata = c * (gal_z - clus_z) / (1 + clus_z)

                cut2 = np.where((np.abs(vdata)<500)&(rdata<1.5))[0]
                clus_ra = astats.biweight_location(gal_ra[cut2])
                clus_dec = astats.biweight_location(gal_dec[cut2])
                clus_z = astats.biweight_location(gal_z[cut2])

                ang_d,lum_d = C.zdistance(clus_z,H0)
                angles = C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
                rdata = angles * ang_d
                vdata = c * (gal_z - clus_z) / (1 + clus_z)

		cut3 = np.where((np.abs(vdata)<500)&(rdata<1.5))[0]
		RA_AVG.append(astats.biweight_location(gal_ra[cut3]))
		DEC_AVG.append(astats.biweight_location(gal_dec[cut3]))
		Z_AVG.append(astats.biweight_location(gal_z[cut3]))

	RA_AVG,DEC_AVG,Z_AVG = np.array(RA_AVG),np.array(DEC_AVG),np.array(Z_AVG)

	write_file = True
	if write_file == True:
		d = {}
		for j in halos.names:
			d[j] = halos[j]
		d['RA_AVG'] = RA_AVG
		d['DEC_AVG'] = DEC_AVG
		d['Z_AVG'] = Z_AVG

		fits_table(d,halos.names+['RA_AVG','DEC_AVG','Z_AVG'],data_loc+'/halos2.fits')

## Calculate Cluster Richnesses
cluster_rich = False
if cluster_rich == True:

        root = '/nfs/christoq_ls/nkern'
        data_loc = 'MassRich/TRUTH_CAUSTIC'
        write_loc = 'individual'

        C4 = CFOUR({'H0':70,'chris_data_root':'/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC'})
	C = Caustic()
	H0 = 72.0
	c = 2.99792e5
	Cosmo = cosmo.LambdaCDM(H0,0.3,0.7)
	keys = ['C','H0','c','Cosmo','C4']
	varib = ez.create(keys,locals())
	R = RICHNESS(varib)


        # Load Halos
        halos = fits.open(root+'/C4/'+data_loc+'/halos.fits')[1].data
        HaloID = halos['orig_order']
        RA = halos['ra_avg']
        DEC = halos['dec_avg']
        Z = halos['z_avg']
        RVIR = halos['RVIR']

	# Load Galaxy Data
	gals = fits.open('/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC/m19.1_allgals_wsdss_specerrs_abs.fits')[1].data

	# Derive Gal Cut Out Parameters
	arcs = np.array(Cosmo.arcsec_per_kpc_proper(Z))*15000. / 3600.    # 15 Mpc in degrees

	richness = []
	# Loop over clusters
	for i in range(len(HaloID)):

		if i % 100 == 0: print i

		# Define Cluster parameters
		clus_ra = RA[i]
		clus_dec = DEC[i]
		clus_z = Z[i]
		clus_rvir = RVIR[i]
		haloid = HaloID[i]

		d_dec = arcs[i]
		d_ra = d_dec/np.cos(clus_dec*np.pi/180)
		d_z = 0.03	

		# Cut Out Galaxy Data Around Cluster
		cut = np.where( (np.abs(gals['gal_ra']-clus_ra)<d_ra) & (np.abs(gals['gal_dec']-clus_dec)<d_dec) & (np.abs(gals['gal_z_app']-clus_z) < d_z))[0]
		gal_ra = gals['gal_ra'][cut];gal_dec = gals['gal_dec'][cut];gal_z=gals['gal_z_app'][cut];gal_gmags=gals['gal_sdss_g'][cut];gal_rmags=gals['gal_sdss_r'][cut];gal_imags=gals['gal_sdss_i'][cut];gal_absr=gals['r_absmag'][cut]

		# Run Richness Estimator
		rich = R.richness_est(gal_ra,gal_dec,gal_z,np.zeros(len(gal_z)),gal_gmags,gal_rmags,gal_imags,gal_absr,haloid,clus_ra,clus_dec,clus_z,clus_rvir=clus_rvir,spec_list=None,use_specs=False,use_bcg=False,fit_rs=False,fixed_vdisp=False,fixed_aperture=False,plot_sky=False,plot_gr=False,plot_phase=False)

		richness.append(rich)


	richness = np.array(richness)

	write_data = True
	if write_data == True:
		d = {}
		d['orig_order'] = HaloID
		d['N200'] = richness
		d['RVIR'] = RVIR
		d['M200'] = halos['m200crit']
		fits_table(d,['orig_order','N200','RVIR','M200'],data_loc+'/kern_richnesses.fits')



## Do Cluster Pair Elimination Technique
def pair_find(rdata,vdata):
	# Set Binning Constants
	rlim = 3		# Mpc
	vlim = 4000		# km/s
	vstep = 1000		# km/s
	nbin = vlim/vstep*2

	# Do binning and find number of points in each bin
	bins_temp = np.zeros(nbin)
	for i in range(nbin):
		bins_temp[i] = np.where((rdata<rlim)&(vdata<(vlim-i*vstep))&(vdata>(vlim-(i+1)*vstep)))[0].size

	# Make bin gridding denser by factor of 2, interpolate for higher resolution
	bins = np.zeros(nbin*2+1)
	j = 0
	for i in range(len(bins)):
		if i == 0 or i == range(len(bins))[-1]:
			bins[i] = 0
		elif i % 2 != 0:
			bins[i] = bins_temp[j]
		else:
			bins[i] = np.mean([bins_temp[j],bins_temp[j+1]])
			j += 1

	# Create 3 double peaked models and one single peaked model
	v_range = np.linspace(vlim,-vlim,nbin*2+1)
	peak = np.max(bins)
	double1 = np.zeros(nbin*2+1)		# clusters at +/- 1500 km/s
	double1[nbin/2:nbin/2+3] = peak
	double1[-nbin/2-3:-nbin/2] = peak
	double2 = np.zeros(nbin*2+1)		# clusters at 0km/s and -3000 km/s
	double2[nbin-1:nbin+2] = peak
	double2[-4:-1] = peak
	double3 = np.zeros(nbin*2+1)		# clusters at 0 km/s and 3000 km/s
	double3[nbin-1:nbin+2] = peak
	double3[1:4] = peak
	single = np.zeros(nbin*2+1)
	single[nbin/2+3:nbin/2+6] = peak
	single[[nbin/2+2,nbin/2+6]] = peak/2
	
	# Do Simple Least Squares Fit
	d1_chi = sum((bins-double1)**2)
	d2_chi = sum((bins-double2)**2)
	d3_chi = sum((bins-double3)**2)
	s_chi = sum((bins-single)**2)

	# If any double model chi_sq is < s_chi*3, it is likely a pair
	thresh = 2
	if d1_chi < s_chi*thresh or d2_chi < s_chi*thresh or d3_chi < s_chi*thresh:
		pair = True
	else:
		pair = False	

	return pair, d1_chi, d2_chi, d3_chi, s_chi, double1, double2, double3, single, v_range, bins




