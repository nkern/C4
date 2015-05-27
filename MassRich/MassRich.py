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

## Constants




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
	truth = fits.open(data_root+'/richness/c4_cluster_truth_revH70_rev5.fits')[1].data
	truth['m200mean']*=1e10
	truth['m200crit']*=1e10

	Rich = truth['n10virm2025']
	Mass = truth['m200crit']

	rich_low,rich_high,rich_step = 2.8,4.4,0.3
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


## Truth Caustic Table (M-R for Henriques using Projected Members and Caustic Mass) ##
panel2 = False
if panel2 == True:
	# Constants
	data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC'
	root = '/nfs/christoq_ls/nkern/C4/MassRich/TRUTH_CAUSTIC'

	# Load Data
	halos = fits.open(data_root+'/c4_cluster_single_whalo.fits')[1].data

	M200_true = halos['M200CRIT']
	R200_true = halos['R200CRIT']

	Rich = truth['n10virm2025']
	Mass = truth['m200crit']

	rich_low,rich_high,rich_step = 2.8,4.4,0.3
	richbins = np.arange(rich_low,rich_high,rich_step)
	mass_low,mass_high,mass_step = 14.25,15.15,0.2
	massbins = np.arange(mass_low,mass_high,mass_step)

	d = stats(Rich,Mass,rich_low,rich_high,rich_step,log10=False,fitline=True)
	globals().update(d)



## Truth Caustic C4 Table (M-R for C4 Catalogue in Henriques using Projected Members and Caustic Masses) ##
panel3 = False
if panel3 == True:
	# Constants
	data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC_C4'
	root = '/nfs/christoq_ls/nkern/C4/MassRich/TRUTH_CAUSTIC_C4'
        H0 = 70
        h = H0 / 100.
        c = 2.99e5

	# Load Data
	halos = fits.open(data_root+'/c4_cluster_single_whalo.fit')[1].data

	M200_true = halos['M200CRIT']
	R200_true = halos['R200CRIT']





## SDSS C4 (M-R for C4 Catalogue in Real Universe using SDSS DR12 spectra / photometry and Caustic Masses) ##
panel4 = True
if panel4 == True:
	# Constants
	data_root = '/nfs/christoq_ls/C4/sdssdr12'
	root = '/nfs/christoq_ls/nkern/C4/MassRich/DR12'
	H0 = 70
	h = H0 / 100.
	c = 2.99e5

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













