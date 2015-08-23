'''
CatPrep.py
Construction / Preparation of Halo Catalogue for all four panels for stacking
Individual M200 estimate is done after catalogue is constructed, and lives in data_loc/individual
For stacking analysis, both halo catalogue and individual mass estimates (if they exist) are loaded

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
from fits_table import *
from c4_richness import RICHNESS
import astropy.cosmology as cosmo
import DictEZ as ez
from stack_class_c4 import CFOUR
from causticpy import Caustic
import cPickle as pkl
import subprocess

## Constants
c = 2.99792458e5		#km/s
H0 = 70


## Flags ##
# Which Panel to run on ?
# panel1 = TRUTH, panel2 = TRUTH_CAUSTIC, panel3 = TRUTH_CAUSTIC_C4, panel4 = DR12
panel1 = False
panel2 = True
panel3 = False
panel4 = False

# What calculations to make ?
init_cats = False		# Make initial catalogues with full Henriques data (full_halos.fits) and one with basic clean cuts (halos.fits)
calc_avg = False		# Calculate average halo centers using spectra
calc_rich = False		# Calculate miller_N200 and kern_N200 and cluster_pairs etc.
calc_cau = False		# Run caustic technique over individual clusters
write_cau = False		# Load caustic masses and write out to halos.fits
sel_func = False		# Perform cleaning selection based on redshifts, cluster pairs, etc.

# After catalogue is made
subsample = False


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

## Projection and Averaging Code ##
def proj_avg(clus_ra,clus_dec,clus_z,gal_ra,gal_dec,gal_z,vlim,rlim,C):

	# Project Galaxies
	ang_d,lum_d = C.zdistance(clus_z,H0)
	angles = C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
	rdata = angles * ang_d
	vdata = c * (gal_z - clus_z) / (1 + clus_z)

	# Take Average
	cut = np.where((np.abs(vdata)<vlim)&(rdata<1.5))[0]
	if len(cut) < 2: cut = np.where((np.abs(vdata)<vlim+500)&(rdata<rlim))[0]
	if len(cut) < 2:
		clus_ra = np.median(gal_ra[cut])
		clus_dec = np.median(gal_dec[cut])
		clus_z = np.median(gal_z[cut])
	else:
		clus_ra = astats.biweight_location(gal_ra[cut])
		clus_dec = astats.biweight_location(gal_dec[cut])
		clus_z = astats.biweight_location(gal_z[cut])

	return clus_ra, clus_dec, clus_z


## Calculate Average Halo Positions from Galaxies ##
def clus_avg(data_loc,halo_file,chris_data_root,newfilename,write_data=True,clobber=True):

	C4 = CFOUR({'H0':70,'chris_data_root':chris_data_root})
	C = Caustic()

	# Load Halos
        halos = fits.open(data_loc+'/'+halo_file)[1].data
        HaloID = halos['orig_order']
        RA = halos['ra_bcg']
        DEC = halos['dec_bcg']
        Z = halos['z_biwt']
        RVIR = halos['RVIR']
        SINGLE = halos['single']
        SUB = halos['sub']
        NC4 = halos['nc4']

	RA_AVG,DEC_AVG,Z_AVG = [],[],[]

	# Loop Over Halos
	print ''
	print '-'*40
	print '...running average cluster center code'
	for i in range(len(halos)):
		if i % 100 == 0: print '...working on cluster '+str(i)+' out of '+str(len(halos))
		try:
			# Assign Halo Properties
			clus_ra = RA[i]
			clus_dec = DEC[i]
			clus_z = Z[i]

			# Load Galaxies
			galdata = C4.load_chris_gals(HaloID[i])
			gal_ra,gal_dec,gal_z,gal_gmags,gal_rmags,gal_imags = galdata

			# Take Iterative Average, four times
			# vlim = 1500, rlim = 1.5
			clus_ra,clus_dec,clus_z = proj_avg(clus_ra,clus_dec,clus_z,gal_ra,gal_dec,gal_z,1500,1.5,C)
			# vlim = 1000, rlim = 1.5
			clus_ra,clus_dec,clus_z = proj_avg(clus_ra,clus_dec,clus_z,gal_ra,gal_dec,gal_z,1000,1.5,C)
			# vlim = 1500, rlim = 1.5
			clus_ra,clus_dec,clus_z = proj_avg(clus_ra,clus_dec,clus_z,gal_ra,gal_dec,gal_z,1000,1.5,C)
			# vlim = 2000, rlim = 1.5
			clus_ra,clus_dec,clus_z = proj_avg(clus_ra,clus_dec,clus_z,gal_ra,gal_dec,gal_z,2000,1.5,C)

		except:
			print i
			clus_ra,clus_dec,clus_z = 0, 0, 0

		RA_AVG.append(clus_ra)
		DEC_AVG.append(clus_dec)
		Z_AVG.append(clus_z)

	RA_AVG,DEC_AVG,Z_AVG = np.array(RA_AVG),np.array(DEC_AVG),np.array(Z_AVG)

	print '...finished average cluster-center calculations'

	## Write Data Out
	if write_data == True:
		print '...writing out cluster catalgoue with average centers included'
		# Dictionary of new columns
		new_keys = ['RA_AVG','DEC_AVG','Z_AVG']
		new_dic = ez.create(new_keys,locals())

		# Original fits record file
		orig_table = halos

		# Write own fits file
		keys = ['HaloID','RA','DEC','Z','RVIR','RA_AVG','DEC_AVG','Z_AVG']
		dic = ez.create(keys,locals())
		fits_table(dic,keys,data_loc+'/avg_centers.fits',clobber=True)

		# Append new columns
		fits_append(orig_table,new_dic,new_keys,filename=data_loc+'/'+newfilename,clobber=clobber)
		print '-'*40
		print ''


## Cluster Richness ##
def cluster_rich(data_loc,halo_file,chris_data_root,chris_data_file,newfilename,write_data=True,clobber=True):
	""" Calculate a Cluster's richness via Miller N200 and Kern N200
		data_loc : e.g. MassRich/TRUTH_CAUSTIC
		halo_file : e.g. halos.fits
		chris_data_root : e.g. /nfs/christoq_ls/C4/sdssdr12
		chris_data_file : e.g. DR12_GalaxyPhotoData_wabs_wedges.fits or m19.1_allgals_wsdss_specerrs_abs.fits
 """

	# Run through Kern richness estimator, mainly to get cluster pairs

        root = '/nfs/christoq_ls/nkern'

        C4 = CFOUR({'H0':70,'chris_data_root':chris_data_root})
	C = Caustic()
	H0 = 70.0
	c = 2.99792e5
	Cosmo = cosmo.LambdaCDM(H0,0.3,0.7)
	keys = ['C','H0','c','Cosmo','C4']
	varib = ez.create(keys,locals())
	R = RICHNESS(varib)

        # Load Halos
        halos = fits.open(data_loc+'/'+halo_file)[1].data
        HaloID = halos['orig_order']
        RA = halos['ra_avg']
        DEC = halos['dec_avg']
        Z = halos['z_avg']
        RVIR = halos['RVIR']

	# Load Galaxy Data
	gals = fits.open(chris_data_root+'/'+chris_data_file)[1].data
	# Change gals keys according to SDSSDR12 galaxy file
	if data_loc[-4:] != 'DR12':
		gals = dict(map(lambda x: (x,gals[x]),gals.names))
		gals['objid'] = gals.pop('GAL_HALOID')
		gals['ra'] = gals.pop('GAL_RA')
		gals['dec'] = gals.pop('GAL_DEC')
		gals['z'] = gals.pop('GAL_Z_APP')
		gals['u_mag'] = gals.pop('GAL_SDSS_U')
		gals['g_mag'] = gals.pop('GAL_SDSS_G')
		gals['r_mag'] = gals.pop('GAL_SDSS_R')
		gals['i_mag'] = gals.pop('GAL_SDSS_I')
		gals['z_mag'] = gals.pop('GAL_SDSS_Z')
		gals['r_absmag'] = gals.pop('R_ABSMAG')

	# Derive Gal Cut Out Parameters
	arcs = np.array(Cosmo.arcsec_per_kpc_proper(Z))*15000. / 3600.    # 15 Mpc in degrees

	# Kern richness arrays
	kern_N200 = []
	HVD = []
	pair_avg = []
	Nspec = []
	kern_obs_tot = []
	kern_obs_back = []		# obs_back is number of non-member galaxies in central aperture around cluster (aka. already scaled to inner aperture)

	# Miller Richness Arrays
	new = fits.open(chris_data_root+'/richness_cr200_bcg/new.fits')[0].data
	newb = fits.open(chris_data_root+'/richness_cr200_bcg/newb.fits')[0].data
	newb *= 0.2
	miller_N200 = []
	miller_obs_tot = []
	miller_obs_back = []
	colfac = 9
	bakfac = 1
	mag3 = 4
	v = 2
	radfac = 1

	# Loop over clusters
	print ''
	print '-'*40
	print '...calculating cluster richnesses'
	for i in range(len(HaloID)):

		if i % 100 == 0: print '...working on cluster '+str(i)+' out of '+str(len(HaloID))
		# Define Cluster parameters
		clus_ra = RA[i]
		clus_dec = DEC[i]
		clus_z = Z[i]
		clus_rvir = RVIR[i]
		haloid = HaloID[i]

		if np.isnan(clus_ra) == True or np.isnan(clus_dec) == True or np.isnan(clus_z) == True:
			richness.append(0)
			HVD.append(0)
			pair_avg.append(False)
			Nspec.append(0)
			continue

		# 15 Mpc in degrees of declination and degrees of RA
		d_dec = arcs[i]
		d_ra = d_dec/np.cos(clus_dec*np.pi/180)
		d_z = 0.04

		# Cut Out Galaxy Data Around Cluster
		cut = np.where( (np.abs(gals['ra']-clus_ra)<d_ra) & (np.abs(gals['dec']-clus_dec)<d_dec) & (np.abs(gals['z']-clus_z) < d_z))[0]
		gal_ra = gals['ra'][cut];gal_dec = gals['dec'][cut];gal_z=gals['z'][cut];gal_gmags=gals['g_mag'][cut];gal_rmags=gals['r_mag'][cut];gal_imags=gals['i_mag'][cut];gal_absr=gals['r_absmag'][cut]

		# Run Kern Richness Estimator
		rich = R.richness_est(gal_ra,gal_dec,gal_z,np.zeros(len(gal_z)),gal_gmags,gal_rmags,gal_imags,gal_absr,haloid,clus_ra,clus_dec,clus_z,clus_rvir=clus_rvir,spec_list=None,use_specs=False,use_bcg=False,fit_rs=False,fixed_vdisp=False,fixed_aperture=False,plot_sky=False,plot_gr=False,plot_phase=False,find_pairs=True)

		kern_N200.append(rich)
		HVD.append(R.vel_disp)
		pair_avg.append(R.pair)
		Nspec.append(R.Nspec)
		kern_obs_tot.append(R.obs_tot)
		kern_obs_back.append(R.obs_back_scaled)

		# Append Miller Richness Values
		k = halos['orig_order'][i]
		miller_N200.append(new[k,colfac,mag3,v,radfac] - newb[k,bakfac,mag3,v,radfac])
		miller_obs_tot.append(new[k,colfac,mag3,v,radfac])
		miller_obs_back.append(newb[k,bakfac,mag3,v,radfac])


	kern_N200 = np.array(kern_N200)
	HVD = np.array(HVD)
	pair_avg = np.array(pair_avg)
	Nspec = np.array(Nspec)
	kern_obs_tot = np.array(kern_obs_tot)
	kern_obs_back = np.array(kern_obs_back)

	miller_N200 = np.array(miller_N200)
	miller_obs_tot = np.array(miller_obs_tot)
	miller_obs_back = np.array(miller_obs_back)

	print '...finished calculating richnesses'
	## Write Data Out
	if write_data == True:
		print '...writing out halos.fits file'

		# Dictionary of new columns
		new_keys = ['kern_N200','HVD','pair_avg','Nspec','kern_obs_tot','kern_obs_back','miller_N200','miller_obs_tot','miller_obs_back']
		new_dic = ez.create(new_keys,locals())

		# Original fits record file
		orig_table = halos

		# Write out own fits file
		keys = ['HaloID','RVIR'] + new_keys
		dic = ez.create(keys,locals())
		fits_table(dic,keys,data_loc+'/richnesses.fits',clobber=True)

		# Append new columns
		fits_append(orig_table,new_dic,new_keys,filename=data_loc+'/'+newfilename,clobber=clobber)
		print '-'*40
		print ''


def ind_caustic(params):
	''' 
		- run caustic technique on individual clusters from halo catalogue
		params : caustic_params.py values, in a dictionary
	'''

	print ''
	print '-'*40
	print '...running caustic technique over halos.fits'

	# Create caustic_params.py file in /C4/ directory using bash command 'sed' and pre-defined parameters
	sed_string = 'sed -e "'
	for i in params.keys():
		sed_string += 's:@@'+i+'@@:'+str(params[i])+':g;'

	sed_string += '" < caustic_params_pbs.py > caustic_params.py'
	
	# Run sed command from bash shell
	print '...creating caustic_params.py file'
	os.system(sed_string)
	print '...C4/caustic_params.py file overwritten and created'

	# Run c4_stack.py
	print ''
	print ''
	print '...running c4_stack.py over '+str(params['halo_num'])+' clusters'
	os.system('python c4_stack.py 0')
	print ''
	print ''
	print '...ran c4_stack.py through!'
	print ''

	file_num = int(subprocess.check_output('ls -1 '+str(params['data_loc'])+'/'+str(params['write_loc'])+' | wc -l', shell=True))
	print '...number of files in '+str(params['data_loc'])+'/'+str(params['write_loc'])+' is '+str(file_num)


def write_m200(params,newfilename,clobber=True):
	''' load individual caustic mass estimates from individual/ and add to halos.fits '''
	# Add M200 estimates to halos.fits
	print '...loading in files from '+str(params['data_loc'])+'/'+str(params['write_loc'])
	CAUM200 = []
	CAUHVD = []
	CAUM200_EST = []
	CAUR200_EST = []
	CAUM500_EST = []
	CAUR500_EST = []
	for i in range(params['halo_num']):
		sys.stdout.write("Progress... "+str(i)+" out of "+str(params['halo_num'])+"\r")
		sys.stdout.flush()
		try:
			f = open(data_loc+'/'+write_loc+'/Ensemble_'+str(i)+'_Data.pkl','rb')
			input = pkl.Unpickler(f)
			data = input.load()
			CAUM200.append(float(data['ens_caumass']))
			CAUHVD.append(float(data['ens_hvd']))
			CAUM200_EST.append(float(data['ens_caumass_est']))
			CAUR200_EST.append(float(data['ens_r200_est']))
			CAUM500_EST.append(float(data['ens_caumass500_est']))
			CAUR500_EST.append(float(data['ens_r500_est']))
		except:
			CAUM200.append(0)
			CAUHVD.append(0)
			CAUM200_EST.append(0)
			CAUR200_EST.append(0)
			CAUM500_EST.append(0)
			CAUR500_EST.append(0)

	CAUM200 = np.array(CAUM200)
	CAUHVD = np.array(CAUHVD)
	CAUM200_EST = np.array(CAUM200_EST)
	CAUR200_EST = np.array(CAUR200_EST)
	CAUM500_EST = np.array(CAUM500_EST)
	CAUR500_EST = np.array(CAUR500_EST)


	# Define Catastrophic Failures as M200 < 1e12
	CATA_FAIL = np.array([False]*len(params['halo_num']))
	CATA_FAIL[np.where((CAUM200<1e12)|(np.isnan(CAUM200)==True))] = True

	print '...finished load'
	## Write Data Out
	if write_data == True:
		print '...writing out halos.fits file'

		# Dictionary of new columns
		new_keys = ['CAUM200','CAUHVD','CAUM200_EST','CAUR200_EST','CAUM500_EST','CAUR500_EST','CATA_FAIL']
		new_dic = ez.create(new_keys,locals())

		# Original fits record file
		orig_table = halos

		# Append new columns
		fits_append(orig_table,new_dic,new_keys,filename=data_loc+'/'+newfilename,clobber=clobber)
		print '-'*40
		print ''


def select_func(params):
	''' Takes cluster catalogue and makes cleaning selection cuts based on redshift, N200, Nspec, cluster pairs etc. '''

	print '...beginning clean selection cuts'	

	# Load Halo Catalgoue
	halos = fits.open(data_loc+'/halos.fits')[1].data

	if params['data_loc'][-2:] == 'TH':		# TRUTH
		cut = np.where(	
				(halos['Z_BCG']>0.03)&
				(halos['Z_BCG']<0.12)&
				(halos['n10virm19']>5)
			)[0]

	if params['data_loc'][-2:] == 'IC':		# TRUTH_CAUSTIC
		cut = np.where(	
				(halos['miller_N200']>5)&
				(np.isnan(halos['miller_N200'])==False)&
				(halos['Z_AVG']>0.03)&
				(halos['Z_AVG']<0.12)&
				(halos['Nspec']>5)&
				(halos['pair_avg']!=True)&
				(halos['cata_fail']!=True)
			)[0]

	elif params['data_loc'][-2:] == 'C4':		# TRUTH_CAUSTIC_C4
		cut = np.where(	
				(halos['miller_N200']>5)&
				(np.isnan(halos['miller_N200'])==False)&
				(halos['Z_AVG']>0.03)&
				(halos['Z_AVG']<0.12)&
				(halos['Nspec']>5)&
				(halos['pair_avg']!=True)&
				(halos['cata_fail']!=True)
			)[0]

	elif params['data_loc'][-2:] == '12':		# DR12
		cut = np.where(	
				(halos['miller_N200']>=5)&
				(np.isnan(halos['miller_N200'])==False)&
				(halos['Z_AVG']>0.03)&
				(halos['Z_AVG']<0.12)&
				(halos['Nspec']>5)&
				(halos['pair_avg']!=True)&
				(halos['cata_fail']!=True)
			)[0]

	print '...size before cuts = '+str(len(halos))+', size after = '+str(len(cut))
	
	new_halos = halos[cut]
	keys = new_halos.names
	dic = dict(map(lambda x: (x,new_halos[x]),keys))
	fits_table(dic,keys,data_loc+'/cut_halos.fits',clobber=True)




#####################################################
################ P R O G R A M ######################
#####################################################



## Truth Table (M-R for Henriques using 3D members and True Mass) ##
if panel1 == True:
	# Constants
	chris_data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH'
	data_loc = 'MassRich/TRUTH'

	# Load Data
	truth = fits.open(chris_data_root+'/c4_cluster_truth_revH100_rev5.fits')[1].data
	truth['m200mean']*=1e10
	truth['m200crit']*=1e10

	# Write out halos.fits file
	if init_cats == True:
		keys = truth.names
		dic = dict(map(lambda x: (x,truth[x]),keys))
		fits_table(dic,keys,data_loc+'/truth.fits',clobber=True)

	params = {'data_loc':data_loc}
	
	## Perform clean selection function cuts ##
	if sel_func == True:
		select_func(params)
	
	Rich = truth['n10virm19']
	Mass = truth['m200crit']

	rich_low,rich_high,rich_step = 3.3,4.8,0.3
	richbins = np.arange(rich_low,rich_high,rich_step)
	mass_low,mass_high,mass_step = 14.25,15.15,0.2
	massbins = np.arange(mass_low,mass_high,mass_step)

	d = stats(Rich,Mass,rich_low,rich_high,rich_step,log10=False,fitline=True)
	globals().update(d)

	plot = False
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
if panel2 == True:

	## Create local halo catalogue straight from Henriques catalogue and C4 catalogue ##
	# Constants
	chris_data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC'
	data_loc = 'MassRich/TRUTH_CAUSTIC'
	chris_photo_file = 'm19.1_allgals_wsdss_specerrs_abs.fits'

        H0 = 70
        h = H0 / 100.

	## Make full halo catalogue from Henriques catalogue
	if init_cats == True:
		# Load Data
		halos = fits.open(chris_data_root+'/c4_cluster_single_whalo.fits')[1].data

		# Assign RVIR as C4['RHO200_RV1500']*0.7**2
		RVIR = halos['RHO200_RV1500']*0.7**2

		# Assign orig_orderID, based on inherent ordering of Chris' native c4_cluster_single.fits file
		ORIG_ORDER = np.arange(len(halos))

		# Write out full halo catalogue before basic cleaning cuts
		fits_append(halos,{'RVIR':RVIR,'ORIG_ORDER':ORIG_ORDER},['RVIR','ORIG_ORDER'],data_loc+'/full_halos.fits')

		## Make secondary halo catalogue based on basic cleaning cuts ##
		# Load Data
		halos = fits.open(data_loc+'/full_halos.fits')[1].data

		# Make Cuts on RVIR detection
		select = np.where((halos['RVIR']>0.2)&(np.isnan(halos['RVIR'])!=True))[0]
		halos = halos[select]

		# Make other possible cuts, like SINGLE, KEEP, and SUB
	#	select = np.where( halos['SINGLE']==9999)[0]
	#	halos = halos[select]
	#	select = np.where( halos['KEEP'] == 1)[0]
	#	halos = halos[select]
	#	select = np.where( halos['SUB'] < 4 )[0]
	#	halos = halos[select]

		# Change all RVIR < .4 to .4
		halos['RVIR'][np.where(halos['RVIR']<.4)] = 0.4

		# Write out halos.fits file
		d = dict(map(lambda x: (x,halos[x]), halos.names))
		fits_table(d,halos.names,data_loc+'/halos.fits',clobber=True)


	## Calculate Average Position of Clusters using surrounding spectra ##
	if calc_avg == True:
		halo_file = 'halos.fits'
		newfilename = 'halos.fits'
		clus_avg(data_loc,halo_file,chris_data_root,newfilename,write_data=True,clobber=True)

	## Calculate Richnesses ##
	if calc_rich == True:
		halo_file = 'halos.fits'
		newfilename = 'halos.fits'
		chris_data_file = chris_photo_file
		cluster_rich(data_loc,halo_file,chris_data_root,chris_data_file,newfilename,write_data=True,clobber=True)


	self_stack	= False
	scale_data	= False
	lightcone	= True
	write_data	= True
	init_shiftgap	= False
	shiftgapper	= True
	edge_int_remove	= True
	small_set	= False
	mass_mix	= False
	bootstrap	= False
	new_halo_cent	= True
	mirror		= True
	true_mems	= False
	run_los		= False
	cent_offset	= None
	
	ens_num		= len(halos)
	gal_num		= 125
	line_num	= 1
	halo_num	= len(halos)
	method_num	= 0
	cell_num	= 0
	table_num	= 0
	data_loc	= data_loc
	write_loc	= 'individual'

	mm_est		= 'richness'
	edge_perc	= 0.10
	center_scat	= None
	avg_meth	= 'median'
	bootstrap_num	= None
	bootstrap_rep	= None
	
	data_set	= ''
	h		= 0.7
	root		= '/nfs/christoq_ls/nkern'

	keys = ['self_stack','scale_data','lightcone','write_data','init_shiftgap','shiftgapper','edge_int_remove','small_set','mass_mix','bootstrap','new_halo_cent','mirror','true_mems','run_los','cent_offset','ens_num','gal_num','line_num','halo_num','method_num','cell_num','table_num','data_loc','write_loc','mm_est','edge_perc','center_scat','avg_meth','bootstrap_num','bootstrap_rep','data_set','h','root']

	params = dict(map(lambda x: (x,globals()[x]), keys))

	## Calculate Individual M200s ##
	if calc_cau == True:
		ind_caustic(params)

	## Write out M200s ##
	if write_cau == True:
		write_m200(params,data_loc+'/halos.fits',clobber=True)

	## Perform clean selection function cuts ##
	if sel_func == True:
		select_func(params)


## Truth Caustic C4 Table (M-R for C4 Catalogue in Henriques using Projected Members and Caustic Masses) ##
if panel3 == True:

	## Create local halo catalogue straight from Henriques catalogue and C4 catalogue ##
	# Constants
	chris_data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC_C4'
	data_loc = 'MassRich/TRUTH_CAUSTIC_C4'
	chris_photo_file = 'm19.1_allgals_wsdss_specerrs_abs.fits'

        H0 = 70
        h = H0 / 100.

	## Make full halo catalogue from Henriques catalogue
	if init_cats == True:
		# Load Data
		halos = fits.open(chris_data_root+'/c4_cluster_single_whalo.fit')[1].data

		# Assign RVIR as C4['RHO200_RV1500']*0.7**2
		RVIR = halos['RHO200_RV1500']*0.7**2

		# Assign orig_orderID, based on inherent ordering of Chris' native c4_cluster_single.fits file
		ORIG_ORDER = np.arange(len(halos))

		# Write out full halo catalogue before basic cleaning cuts
		fits_append(halos,{'RVIR':RVIR,'ORIG_ORDER':ORIG_ORDER},['RVIR','ORIG_ORDER'],data_loc+'/full_halos.fits')

		## Make secondary halo catalogue based on basic cleaning cuts ##
		# Load Data
		halos = fits.open(data_loc+'/full_halos.fits')[1].data

		# Make Cuts on RVIR detection
		select = np.where((halos['RVIR']>0.2)&(np.isnan(halos['RVIR'])!=True))[0]
		halos = halos[select]

		# Make other possible cuts, like SINGLE, KEEP, and SUB
	#	select = np.where( halos['SINGLE']==9999)[0]
	#	halos = halos[select]
	#	select = np.where( halos['KEEP'] == 1)[0]
	#	halos = halos[select]
	#	select = np.where( halos['SUB'] < 4 )[0]
	#	halos = halos[select]

		# Change all RVIR < .4 to .4
		halos['RVIR'][np.where(halos['RVIR']<.4)] = 0.4

		# Write out halos.fits file
		d = dict(map(lambda x: (x,halos[x]), halos.names))
		fits_table(d,halos.names,data_loc+'/halos.fits',clobber=True)


	## Calculate Average Position of Clusters using surrounding spectra ##
	if calc_avg == True:
		halo_file = 'halos.fits'
		newfilename = 'halos.fits'
		clus_avg(data_loc,halo_file,chris_data_root,newfilename,write_data=True,clobber=True)

	## Calculate Richnesses ##
	if calc_rich == True:
		halo_file = 'halos.fits'
		newfilename = 'halos.fits'
		chris_data_file = chris_photo_file
		cluster_rich(data_loc,halo_file,chris_data_root,chris_data_file,newfilename,write_data=True,clobber=True)


	self_stack	= False
	scale_data	= False
	lightcone	= True
	write_data	= True
	init_shiftgap	= False
	shiftgapper	= True
	edge_int_remove	= True
	small_set	= False
	mass_mix	= False
	bootstrap	= False
	new_halo_cent	= True
	mirror		= True
	true_mems	= False
	run_los		= False
	cent_offset	= None
	
	ens_num		= len(halos)
	gal_num		= 125
	line_num	= 1
	halo_num	= len(halos)
	method_num	= 0
	cell_num	= 0
	table_num	= 0
	data_loc	= data_loc
	write_loc	= 'individual'

	mm_est		= 'richness'
	edge_perc	= 0.10
	center_scat	= None
	avg_meth	= 'median'
	bootstrap_num	= None
	bootstrap_rep	= None
	
	data_set	= ''
	h		= 0.7
	root		= '/nfs/christoq_ls/nkern'

	keys = ['self_stack','scale_data','lightcone','write_data','init_shiftgap','shiftgapper','edge_int_remove','small_set','mass_mix','bootstrap','new_halo_cent','mirror','true_mems','run_los','cent_offset','ens_num','gal_num','line_num','halo_num','method_num','cell_num','table_num','data_loc','write_loc','mm_est','edge_perc','center_scat','avg_meth','bootstrap_num','bootstrap_rep','data_set','h','root']

	params = dict(map(lambda x: (x,globals()[x]), keys))

	## Calculate Individual M200s ##
	if calc_cau == True:
		ind_caustic(params)

	## Write out M200s ##
	if write_cau == True:
		write_m200(params,data_loc+'/halos.fits',clobber=True)

	## Perform clean selection function cuts ##
	if sel_func == True:
		select_func(params)



## SDSS C4 (M-R for C4 Catalogue in Real Universe using SDSS DR12 spectra / photometry and Caustic Masses) ##
if panel4 == True:

	## Create local halo catalogue straight from Henriques catalogue and C4 catalogue ##
	# Constants
	chris_data_root = '/nfs/christoq_ls/C4/sdssdr12'
	data_loc = 'MassRich/DR12'
	chris_photo_file = 'DR12_GalaxyPhotoData_wabs_wedges.fits'

        H0 = 70
        h = H0 / 100.

	## Make full halo catalogue from Henriques catalogue
	if init_cats == True:
		# Load Data
		halos = fits.open(chris_data_root+'/c4_cluster_single.fits')[1].data

		# Assign RVIR as C4['RHO200_RV1500']*0.7**2
		RVIR = halos['RHO200_RV1500']*0.7**2

		# Assign orig_orderID, based on inherent ordering of Chris' native c4_cluster_single.fits file
		ORIG_ORDER = np.arange(len(halos))

		# Write out full halo catalogue before basic cleaning cuts
		fits_append(halos,{'RVIR':RVIR,'ORIG_ORDER':ORIG_ORDER},['RVIR','ORIG_ORDER'],data_loc+'/full_halos.fits')

		## Make secondary halo catalogue based on basic cleaning cuts ##
		# Load Data
		halos = fits.open(data_loc+'/full_halos.fits')[1].data

		# Make Cuts on RVIR detection
		select = np.where((halos['RVIR']>0.2)&(np.isnan(halos['RVIR'])!=True))[0]
		halos = halos[select]

		# Make other possible cuts, like SINGLE, KEEP, and SUB
	#	select = np.where( halos['SINGLE']==9999)[0]
	#	halos = halos[select]
	#	select = np.where( halos['KEEP'] == 1)[0]
	#	halos = halos[select]
	#	select = np.where( halos['SUB'] < 4 )[0]
	#	halos = halos[select]

		# Change all RVIR < .4 to .4
		halos['RVIR'][np.where(halos['RVIR']<.4)] = 0.4

		# Write out halos.fits file
		d = dict(map(lambda x: (x,halos[x]), halos.names))
		fits_table(d,halos.names,data_loc+'/halos.fits',clobber=True)

	## Calculate Average Position of Clusters using surrounding spectra ##
	if calc_avg == True:
		halo_file = 'halos.fits'
		newfilename = 'halos.fits'
		clus_avg(data_loc,halo_file,chris_data_root,newfilename,write_data=True,clobber=True)

	## Calculate Richnesses ##
	if calc_rich == True:
		halo_file = 'halos.fits'
		newfilename = 'halos.fits'
		chris_data_file = chris_photo_file
		cluster_rich(data_loc,halo_file,chris_data_root,chris_data_file,newfilename,write_data=True,clobber=True)


	self_stack	= False
	scale_data	= False
	lightcone	= True
	write_data	= True
	init_shiftgap	= False
	shiftgapper	= True
	edge_int_remove	= True
	small_set	= False
	mass_mix	= False
	bootstrap	= False
	new_halo_cent	= True
	mirror		= True
	true_mems	= False
	run_los		= False
	cent_offset	= None
	
	ens_num		= len(halos)
	gal_num		= 125
	line_num	= 1
	halo_num	= len(halos)
	method_num	= 0
	cell_num	= 0
	table_num	= 0
	data_loc	= data_loc
	write_loc	= 'individual'

	mm_est		= 'richness'
	edge_perc	= 0.10
	center_scat	= None
	avg_meth	= 'median'
	bootstrap_num	= None
	bootstrap_rep	= None
	
	data_set	= ''
	h		= 0.7
	root		= '/nfs/christoq_ls/nkern'

	keys = ['self_stack','scale_data','lightcone','write_data','init_shiftgap','shiftgapper','edge_int_remove','small_set','mass_mix','bootstrap','new_halo_cent','mirror','true_mems','run_los','cent_offset','ens_num','gal_num','line_num','halo_num','method_num','cell_num','table_num','data_loc','write_loc','mm_est','edge_perc','center_scat','avg_meth','bootstrap_num','bootstrap_rep','data_set','h','root']

	params = dict(map(lambda x: (x,globals()[x]), keys))

	## Calculate Individual M200s ##
	if calc_cau == True:
		ind_caustic(params)

	## Write out M200s ##
	if write_cau == True:
		write_m200(params,data_loc+'/halos.fits',clobber=True)

	## Perform clean selection function cuts ##
	if sel_func == True:
		select_func(params)








## Sub Sampled Mass Bias Trend ##
'''
Run the caustic technique over individual clusters with high sampling (also stack them)
Run the caustic tehcnique over individual clusters with low sampling (also stack them)
We expect to see little change in the stacked results, large change in individual results

	At first try doing this over panel3 (truth_caustic_c4)
'''
if subsample == True:
	# Set Constants and Flags
	root = '/nfs/christoq_ls/'
	data_loc = 'MILLENNIUM/Henriques/TRUTH_CAUSTIC/'


	data = []
	#data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal100_Nclus20_Nens20.fits')[1].data))
	#data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal75_Nclus20_Nens20.fits')[1].data))
	#data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal50_Nclus20_Nens20.fits')[1].data))
	#data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal25_Nclus20_Nens20.fits')[1].data))
	#data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal15_Nclus20_Nens20.fits')[1].data))
	#data.append(fits_data(fits.open('MassRich/TRUTH_CAUSTIC_C4/gal_subsample/C4Stack_Ngal10_Nclus20_Nens20.fits')[1].data))

	#data.append(pkl.Unpickler(open('MassRich/TRUTH_CAUSTIC/gal_subsample/run1/C4Stack_Ngal75_Nclus50_Nens4.pkl','rb')).load())
        #data.append(pkl.Unpickler(open('MassRich/TRUTH_CAUSTIC/gal_subsample/run2/C4Stack_Ngal50_Nclus50_Nens4.pkl','rb')).load())
        #data.append(pkl.Unpickler(open('MassRich/TRUTH_CAUSTIC/gal_subsample/run3/C4Stack_Ngal25_Nclus50_Nens4.pkl','rb')).load())
        #data.append(pkl.Unpickler(open('MassRich/TRUTH_CAUSTIC/gal_subsample/run4/C4Stack_Ngal15_Nclus50_Nens4.pkl','rb')).load())
        #data.append(pkl.Unpickler(open('MassRich/TRUTH_CAUSTIC/gal_subsample/run5/C4Stack_Ngal10_Nclus50_Nens4.pkl','rb')).load())

        data.append(pkl.Unpickler(open('MassRich/TRUTH_CAUSTIC/gal_subsample/run1c/C4Stack_Ngal25_Nclus15_Nens93.pkl','rb')).load())
        data.append(pkl.Unpickler(open('MassRich/TRUTH_CAUSTIC/gal_subsample/run2c/C4Stack_Ngal15_Nclus15_Nens93.pkl','rb')).load())
        data.append(pkl.Unpickler(open('MassRich/TRUTH_CAUSTIC/gal_subsample/run3c/C4Stack_Ngal10_Nclus15_Nens93.pkl','rb')).load())


	def residual(vars, x, data):
		a = vars[0]
		b = vars[1]
		model = a*x + b
		return (model-data)**2

	def residual2(vars, x, data):
		b = vars[0]
		odel = 0.952 * x + b
		return (model-data)**2

	xdata = np.log10(data[0]['BINN200'])
	ydata = np.log10(data[0]['EDGEMASS'])
	vars = [1,14]
	out = leastsq(residual, vars, args=(xdata,ydata))[0]

	plot = False
	if plot == True:
#		fig,ax = mp.subplots(1,2,figsize=(14,7))
		fig,ax = mp.subplots()
		ax = [ax]
		p0, = ax[0].plot(xdata,np.log10(data[0]['EDGEMASS']),'m.',markersize=9)
		p1, = ax[0].plot(xdata,np.log10(data[1]['EDGEMASS']),'b.',markersize=9)
		p2, = ax[0].plot(xdata,np.log10(data[2]['EDGEMASS']),'g.',markersize=9)
		#p3, = ax[0].plot(xdata,np.log10(data[3]['EDGEMASS']),'co',markersize=9)
		#p4, = ax[0].plot(xdata,np.log10(data[4]['EDGEMASS']),'yo',markersize=9)
#		p5, = ax[0].plot(xdata,np.log10(data[5]['EDGEMASS']),'ro',markersize=9)
#		ax[0].legend([p0,p1,p2,p3,p4],['Ngal=75','Ngal=50','Ngal=25','Ngal=15','Ngal=10'],loc=2)

		out = leastsq(residual, vars, args=(np.log10(data[0]['BINN200']),np.log10(data[0]['EDGEMASS'])))[0]
#		out1 = leastsq(residual2, [13], args=(xdata,np.log10(data[1]['EDGEMASS'])))[0]
#		out2 = leastsq(residual2, [13], args=(xdata,np.log10(data[2]['EDGEMASS'])))[0]
#		out3 = leastsq(residual2, [13], args=(xdata,np.log10(data[3]['EDGEMASS'])))[0]
#		out4 = leastsq(residual2, [13], args=(xdata,np.log10(data[4]['EDGEMASS'])))[0]
##		out5 = leastsq(residual2, [13], args=(xdata,np.log10(data[5]['EDGEMASS'])))[0]
#		out_params = [out[1],out1[0],out2[0],out3[0],out4[0],]
		x = np.arange(1.4,2.6,.1)
	
#		p0, = ax[1].plot(x,out[0]*x + out[1],'m')
#		p1, = ax[1].plot(x,out[0]*x + out1[0],'b')
#		p2, = ax[1].plot(x,out[0]*x + out2[0],'g')
#		p3, = ax[1].plot(x,out[0]*x + out3[0],'c')
#		p4, = ax[1].plot(x,out[0]*x + out4[0],'y')
#		p5, = ax[1].plot(x,out[0]*x + out5[0],'r')

#		bias = [np.median(np.log( 10**(out[0]*x+out[1]) / 10**(out[0]*x+i) )) for i in out_params]

#		ax[1].legend([p0,p1,p2,p3,p4],['-'+str(np.round(bias[0]*100,1))+'% bias','-'+str(np.round(bias[1]*100,1))+'% bias','-'+str(np.round(bias[2]*100,1))+'% bias','-'+str(np.round(bias[3]*100,1))+'% bias','-'+str(np.round(bias[4]*100,1))+'% bias',],loc=2)


		di = 0.043
		l1, = ax[0].plot(x,out[0]*x + out[1],'m')
		l2, = ax[0].plot(x,out[0]*x + out[1]+di,'r')
		ax[0].plot(x,out[0]*x + out[1]-di,'r')
		i = out[1] + di
		bias = np.mean(np.log( 10**(out[0]*x+out[1]) / 10**(out[0]*x+i) ))
		bias = np.round(np.abs(bias)*100,1)
		ax[0].legend([p0,p1,p2,l1,l2],['Ngal=25','Ngal=15','Ngal=10','N=25 bestfit','+/-'+str(bias)+'% bias about bestfit'],fontsize=10,loc=2)

	
		ax[0].set_xlim(1,3)
		ax[0].set_ylim(13.5,16)
		#ax[1].set_xlim(1.5,2.5)
		#ax[1].set_ylim(14,15.5)
		ax[0].set_xlabel('log10 N200',fontsize=15)
		ax[0].set_ylabel('log10 M200',fontsize=15)

		savefig = False
		if savefig == True: fig.savefig('T_C_subsample.png',dpi=100,bbox_inches='tight')
	




