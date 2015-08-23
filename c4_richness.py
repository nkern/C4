# Import Modules
import numpy as np
import pylab as mp
import astropy.io.fits as fits
import astropy.cosmology as cosmo
from causticpy import Caustic
import astrostats as astats
from mpl_toolkits.mplot3d.axes3d import Axes3D
import DictEZ as ez
import cPickle as pkl
from sklearn import linear_model
from c4_pairfind import pair_find

# Set Constants
C = Caustic()
H0 = 72.0
c = 2.99792e5
Cosmo = cosmo.LambdaCDM(H0,0.3,0.7)
keys = ['C','H0','c','Cosmo']
varib = ez.create(keys,locals())
root = '/nfs/christoq_ls/nkern'



# Calculate Richnesses
class RICHNESS(object):

	def __init__(self,varib):
		self.__dict__.update(varib)

	def richness_est(self,ra,dec,z,pz,gmags,rmags,imags,abs_rmags,haloid,clus_ra,clus_dec,clus_z,clus_rvir=None,gr_slope=None,gr_width=None,gr_inter=None,shiftgap=True,use_specs=False,use_bcg=False,fit_rs=False,spec_list=None,fixed_aperture=False,fixed_vdisp=False,plot_gr=False,plot_sky=False,plot_phase=False,find_pairs=True):
		''' Richness estimator, magnitudes should be apparent mags except for abs_rmag
			z : spectroscopic redshift
			pz : photometric redshift
			gr_slope : if fed color-mag (g-r vs. r) slope of fit to red sequence
			gr_width : if fed color-mag (g-r vs. r) width of fit to red sequence
			gr_inter : if red color-mag (g-r vs. r) intercept of fit to red sequence
			spec_list : indexing array specifying gals w/ spectroscopic redshifts, preferably fed as boolean numpy array
			fixed_aperture : clus_rvir = 1 Mpc
			fixed_vdisp : clus_vdisp = 1000 km/s (members < clus_vdisp*2)
			use_specs : includes all spectro members in richness regardless of their color, subject to counting twice
			use_bcg : use bcg RS_color as cluster RS_color
			fit_rs : fit a line to the member galaxy RS relationship, as of now this does not work and we take a flat RS relationship
			plot_gr : plot r vs g-r color diagram with RS fits
			plot_sky : plot RA and DEC of gals with signal and background annuli
			plot_phase : plot rdata and vdata phase space
			find_pairs : run c4 cluster pair identifier on phase spaces
			- If at least one quarter of the outer annuli doesn't have any galaxies (perhaps because it has run over the observation edge),
				it will return -99 as a richness
		'''
		# Put Into Class Namespace
		keys = ['ra','dec','z','pz','gmags','rmags','imags','clus_ra','clus_dec','clus_z','clus_rvir','vel_disp','color_data','color_cut','RS_color','RS_sigma','clus_color_cut','signal','background','richness','mems','phot','spec','spec_mems','deg_per_rvir','SA_outer','SA_inner','outer_red_dense','inner_background','Nphot_inner','Nphot_outer','red_inner','red_outer','all','outer_edge','inner_edge','shift_cut','set','pair','Nspec','gr_slope','gr_inter','gr_width','obs_tot','obs_back_scaled']
		self.__dict__.update(ez.create(keys,locals()))

		# Take rough virial radius measurement, this method is BAD...
		if clus_rvir == None:
			clus_rvir = np.exp(-1.86)*len(np.where((rmags < -19.55) & (rdata < 1.0) & (np.abs(vdata) < 3500))[0])**0.51
			self.clus_rvir = clus_rvir

		# Set clus_rvir = 1.0 Mpc if fixed aperture == True
		if fixed_aperture == True:
			clus_rvir = 1.0

		# Magnitue Cut at 19.1 Apparent R Mag, then at -19 Absolute R Mag and sort by them
		bright = np.where(rmags < 19.1)[0]
		data = np.vstack([ra,dec,z,pz,gmags,rmags,imags,abs_rmags])
		data = data.T[bright].T
		ra,dec,z,pz,gmags,rmags,imags,abs_rmags = data

                bright = np.where(abs_rmags < -19.5)[0]
                data = np.vstack([ra,dec,z,pz,gmags,rmags,imags,abs_rmags])
                data = data.T[bright].T
                ra,dec,z,pz,gmags,rmags,imags,abs_rmags = data

		sorts = np.argsort(abs_rmags)
                data = np.vstack([ra,dec,z,pz,gmags,rmags,imags,abs_rmags])
		data = data.T[sorts].T
                ra,dec,z,pz,gmags,rmags,imags,abs_rmags = data
		self.__dict__.update(ez.create(keys,locals()))

		# Separate Into Spectro Z and Photo Z DataSets
		# Define the Spectro Z catalogue
		if spec_list == None:
			spec_cut = np.abs(z) > 1e-4
		else:
			if type(spec_list) == list: spec_list = np.array(spec_list)
			if spec_list.dtype != 'bool':
				spec_cut = np.array([False]*len(ra))
				spec_cut[spec_list] = True
		all = {'ra':ra,'dec':dec,'z':z,'pz':pz,'gmags':gmags,'rmags':rmags,'imags':imags,'abs_rmags':abs_rmags}
		spec = {'ra':ra[spec_cut],'dec':dec[spec_cut],'z':z[spec_cut],'pz':pz[spec_cut],'gmags':gmags[spec_cut],'rmags':rmags[spec_cut],'imags':imags[spec_cut],'abs_rmags':abs_rmags[spec_cut]}
		phot = {'ra':ra[~spec_cut],'dec':dec[~spec_cut],'z':z[~spec_cut],'pz':pz[~spec_cut],'gmags':gmags[~spec_cut],'rmags':rmags[~spec_cut],'imags':imags[~spec_cut],'abs_rmags':abs_rmags[~spec_cut]}

		# Project Spectra into radius and velocity
		ang_d,lum_d = C.zdistance(clus_z,self.H0)
		angles = C.findangle(spec['ra'],spec['dec'],clus_ra,clus_dec)
		rdata = angles * ang_d
		vdata = self.c * (spec['z'] - clus_z) / (1 + clus_z)
                spec.update({'rdata':rdata,'vdata':vdata})
                self.__dict__.update(ez.create(keys,locals()))

		# Take Hard Phasespace Limits
		limit = np.where( (rdata < 5) & (np.abs(vdata) < 5000) )[0]
		clus_data = np.vstack([rdata,vdata,spec['ra'],spec['dec'],spec['z'],spec['pz'],spec['gmags'],spec['rmags'],spec['imags'],spec['abs_rmags']])
		clus_data = clus_data.T[limit].T
		rdata,vdata,spec['ra'],spec['dec'],spec['z'],spec['pz'],spec['gmags'],spec['rmags'],spec['imags'],spec['abs_rmags'] = clus_data
                spec.update({'rdata':rdata,'vdata':vdata})
		self.__dict__.update(ez.create(keys,locals()))

		# Shiftgapper for Interlopers
		if shiftgap == True:
			len_before = np.where((rdata < clus_rvir*1.5)&(np.abs(vdata)<4000))[0].size
			clus_data = np.vstack([rdata,vdata,spec['ra'],spec['dec'],spec['z'],spec['pz'],spec['gmags'],spec['rmags'],spec['imags'],spec['abs_rmags']])
			clus_data = C.shiftgapper(clus_data.T).T
			sorts = np.argsort(clus_data[-1])
			clus_data = clus_data.T[sorts].T
			rdata,vdata,spec['ra'],spec['dec'],spec['z'],spec['pz'],spec['gmags'],spec['rmags'],spec['imags'],spec['abs_rmags'] = clus_data
			shift_cut = len_before - np.where((rdata < clus_rvir)&(np.abs(vdata)<4000))[0].size

		# Measure Velocity Dispersion of all galaxies within 1 * r_vir and np.abs(vdata) < 4000 km/s
                spec.update({'rdata':rdata,'vdata':vdata})
                self.__dict__.update(ez.create(keys,locals()))
		vel_disp = astats.biweight_midvariance(vdata[np.where((rdata < clus_rvir*1)&(np.abs(vdata) < 4000))])
		if fixed_vdisp == True: vel_disp = 1000

		# Run Cluster Pair Finder
		if find_pairs == True:
			pair, d1_chi, d2_chi, d3_chi, s_chi, double1, double2, double3, single, v_range, bins = pair_find(rdata,vdata)

		# Calculate Nspec, Nspec = # of galaxies within RVIR
		try:
			Nspec = np.where(rdata<clus_rvir)[0].size
		except:
			Nspec = 0
                self.__dict__.update(ez.create(keys,locals()))

		# Get Members from Specta, get their red sequence color
		mems = np.where((spec['rdata'] < clus_rvir)&(np.abs(spec['vdata'])<2*vel_disp))[0]
		color_data = spec['gmags'] - spec['rmags']
		color_cut = np.where((color_data[mems] < 1.2) & (color_data[mems] > 0.65))[0]
		RS_color = astats.biweight_location(color_data[mems][color_cut])
		RS_shift = color_data[mems][color_cut] - RS_color
		RS_sigma = astats.biweight_midvariance(RS_shift[np.where(np.abs(RS_shift)<.15)])

		if fit_rs == True:
			clf = linear_model.LinearRegression()
			set = np.where(np.abs(color_data[mems]-RS_color)<3*RS_sigma)[0]
			clf.fit(spec['rmags'][mems][set].reshape(set.size,1),color_data[mems][set])

		if use_bcg == True:
			bright = np.argsort(all['abs_rmags'][mems])[0]
			RS_color = color_data[mems][bright]
	                RS_shift = color_data[mems][color_cut] - RS_color
	                RS_sigma = astats.biweight_midvariance(RS_shift[np.where(np.abs(RS_shift)<.15)])
			
		# spec_mems is # of members that are within 2 sigma of cluster color
		clus_color_cut = np.where(np.abs(color_data[mems] - RS_color) < RS_sigma*2)[0]
		spec_mems = len(clus_color_cut)
		self.__dict__.update(ez.create(keys,locals()))

		# If fed gr fit
		if gr_slope != None:
			def RS_color(r_mag): return r_mag*gr_slope + gr_inter
			RS_sigma = gr_width / 2
			clus_color_cut = np.where(np.abs(color_data[mems] - RS_color(spec['rmags'])[mems]) < RS_sigma*2)[0]
			spec_mems = len(clus_color_cut)
			self.__dict__.update(ez.create(keys,locals()))

		# Get Rdata from PhotoZ & SpecZ Data Set
		angles = C.findangle(all['ra'],all['dec'],clus_ra,clus_dec)
		all['rdata'] = angles * ang_d

		# Get deg per RVIR proper
		deg_per_rvir = R.Cosmo.arcsec_per_kpc_proper(clus_z).value * 1e3 * clus_rvir / 3600

		# Get Area of Virial Circle out to 1 rvir, and Area of outer annuli, which is annuli from 4rvir < R < 6rvir, or dependent on rvir
		if clus_rvir < 2.5:
			outer_edge = 6.0
			inner_edge = 4.0
		elif clus_rvir >= 2.5 and clus_rvir < 3:
			outer_edge = 5.0
			inner_edge = 3.5
		else:
			outer_edge = 3.0
			inner_edge = 2.0	
		SA_inner = np.pi*deg_per_rvir**2
		SA_outer = np.pi * ( (outer_edge*deg_per_rvir)**2 - (inner_edge*deg_per_rvir)**2 )

		# Get Number of Cluster Color Galaxies from Photo Data Set in Inner Circle and Outer Annuli
		RS_color_sig = 2.0
		if gr_slope == None:
			red_inner = np.where(((all['gmags']-all['rmags']) < RS_color + RS_color_sig*RS_sigma)&((all['gmags']-all['rmags']) > RS_color - RS_color_sig*RS_sigma)&(all['rdata']<clus_rvir))[0]
			red_outer = np.where(((all['gmags']-all['rmags']) < RS_color + 1.5*RS_sigma)&((all['gmags']-all['rmags']) > RS_color - 1.5*RS_sigma)&(all['rdata']<outer_edge*clus_rvir)&(all['rdata']>inner_edge*clus_rvir))[0]
		else:
			red_inner = np.where(((all['gmags']-all['rmags']) < RS_color(all['rmags']) + RS_color_sig*RS_sigma)&((all['gmags']-all['rmags']) > RS_color(all['rmags']) - RS_color_sig*RS_sigma)&(all['rdata']<clus_rvir))[0]
			red_outer = np.where(((all['gmags']-all['rmags']) < RS_color(all['rmags']) + 1.5*RS_sigma)&((all['gmags']-all['rmags']) > RS_color(all['rmags']) - 1.5*RS_sigma)&(all['rdata']<outer_edge*clus_rvir)&(all['rdata']>inner_edge*clus_rvir))[0]

		Nphot_inner = len(red_inner)
		Nphot_outer = len(red_outer)

		self.__dict__.update(ez.create(keys,locals()))

		# Get Solid Angle Density of Outer Red Galaxies
		outer_red_dense = Nphot_outer / SA_outer
		inner_background = int(np.ceil(outer_red_dense * SA_inner))

		# Define obs_tot and obs_back_scaled, where obs_back_scaled is obs_back already scaled to inner aperture (aka. inner_background)
		obs_tot = np.copy(Nphot_inner)
		obs_back_scaled = np.copy(inner_background)

		# If inner_background is less than Nphot_inner, then Nphot_inner -= inner_background, otherwise Nphot_inner = 0
		if inner_background < Nphot_inner:
			Nphot_inner -= inner_background
		else:
			Nphot_inner = 0

		# Richness = spec_mems + Nphot_inner or just Nphot_inner
		if use_specs == True:
			richness = spec_mems + Nphot_inner
		else:
			richness = Nphot_inner
		self.__dict__.update(ez.create(keys,locals()))

                # Plot
                if plot_gr == True:
                        fig,ax = mp.subplots()
                        ax.plot(rmags,gmags-rmags,'ko',alpha=.8)
                        ax.plot(spec['rmags'][mems],color_data[mems],'co')
			if gr_slope == None:
                        	ax.axhline(RS_color,color='r')
                        	ax.axhline(RS_color+RS_sigma,color='b')
                        	ax.axhline(RS_color-RS_sigma,color='b')
			else:
				xdata = np.arange(spec['rmags'][mems].min(),spec['rmags'][mems].max(),.1)
				ax.plot(xdata,RS_color(xdata),color='r')
				ax.plot(xdata,RS_color(xdata)+RS_sigma,color='b')
				ax.plot(xdata,RS_color(xdata)-RS_sigma,color='b')
                        ax.set_xlim(13,19)
                        ax.set_ylim(0,1.3)
                        ax.set_xlabel('Apparent R Mag',fontsize=16)
                        ax.set_ylabel('App G Mag - App R Mag',fontsize=16)
                        ax.set_title('Color-Mag Diagram, Cluster '+str(haloid))
                        fig.savefig('colormag_'+str(haloid)+'.png',bbox_inches='tight')
                        mp.close(fig)

		if plot_sky == True:
			fig,ax = mp.subplots()
			ax.plot(all['ra'],all['dec'],'ko')
			ax.plot(all['ra'][red_inner],all['dec'][red_inner],'ro')
			ax.plot(all['ra'][red_outer],all['dec'][red_outer],'yo')
			ax.plot(clus_ra,clus_dec,'co',markersize=9)
			ax.set_xlabel('RA',fontsize=15)
			ax.set_ylabel('Dec.',fontsize=15)
			ax.set_title('Richness Annuli for Halo '+str(haloid))
			fig.savefig('skyplot_'+str(haloid)+'.png',bbox_inches='tight')
			mp.close(fig)

		if plot_phase == True:
			fig,ax = mp.subplots()
			ax.plot(spec['rdata'],spec['vdata'],'ko')
			ax.plot(spec['rdata'][mems],spec['vdata'][mems],'co')
			bcg = np.where(spec['abs_rmags']==spec['abs_rmags'][mems].min())[0][0]
			ax.plot(spec['rdata'][bcg],spec['vdata'][bcg],'ro')
			ax.set_xlim(0,5)
			ax.set_ylim(-5000,5000)
			ax.set_xlabel('Radius (Mpc)',fontsize=15)
			ax.set_ylabel('Velocity (km/s)',fontsize=15)
			ax.set_title('phasespace haloid '+str(haloid))
			fig.savefig('phasespace_'+str(haloid)+'.png',bbox_inches='tight')
			mp.close(fig)

                # Check to make sure galaxies exist (somewhat) uniformely in outer annuli (aka. cluster isn't on edge of observation strip)
                # First, project all galaxies into polar coordinates centered on cluster center
		x = (all['ra']-clus_ra)/np.cos(clus_dec*np.pi/180)		# x coordinate in arcsec (of dec) centered on cluster center
		y = (all['dec']-clus_dec)
                all['radius'] = np.sqrt( x**2 + y**2 ) / deg_per_rvir	#radius scaled by RVIR
                all['theta'] = np.arctan( y / x )
		# Add corrections to arctan function
		all['theta'][np.where( (x < 0) & (y > 0) )] += np.pi	# Quadrant II
		all['theta'][np.where( (x < 0) & (y < 0) )] += np.pi	# Quadrant III
		all['theta'][np.where( (x > 0) & (y < 0) )] += 2*np.pi	# Quadrant IV
                # Then break outer annuli into 4 sections and check if at least 1 galaxy exists in each section
                sizes1 = np.array([np.where((np.abs(all['theta']-i)<=np.pi/2)&(all['radius']>inner_edge)&(all['radius']<15))[0].size for i in np.linspace(0,2*np.pi,4)])
                # Do it again but shift theta by np.pi/8 this time
                sizes2 = np.array([np.where((np.abs(all['theta']-i)<=np.pi/2)&(all['radius']>inner_edge)&(all['radius']<15))[0].size for i in np.linspace(np.pi/8,17*np.pi/8,4)])
                if 0 in sizes1 or 0 in sizes2:
			p = False
			if p == True:
				mp.plot(all['radius'],all['theta'],'ko')
				mp.xlabel('radius')
				mp.ylabel('theta')
				mp.savefig('rth_'+str(haloid)+'.png')
				mp.close()
			print 'sizes1=',sizes1
			print 'sizes2=',sizes2
                        return -99

		return richness


R = RICHNESS(varib)







