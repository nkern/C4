'''
This file sets the cluster sample from C4 data and does a whole bunch of other stuff w/ the data
'''
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

# Load Data
# Cluster File
c4 = fits.open('C4_Cluster_Single.fits')[1].data
cut = np.where(c4['z_c4']!=0)
c4 = c4[cut]
sort = np.argsort(c4['z_c4'])
c4 = c4[sort]

# Galaxy File
photo = fits.open('DR12_Photo_AbsMag_Data.fits')[1].data
dr12 = fits.open('DR12_SpectroPhoto_BestObjID.fits')[1].data
dr12 = dr12[np.where((dr12['z'] < 0.27)&(dr12['zWarning']==0))]

# Set Constants
C = Caustic()
H0 = 72.0
c = 2.99792e5
Cosmo = cosmo.LambdaCDM(H0,0.3,0.7)
data_set = 'DataSet2/'
keys = ['C','H0','c','Cosmo','data_set']
varib = ez.create(keys,locals())
root = '/nfs/christoq_ls/nkern'

## Useful Functions
def plot3d():
	fig = mp.figure()
	ax = fig.add_subplot(111,projection='3d')
	globals().update({'ax':ax})


# Write Cluster File
cluster_file = False
if cluster_file == True:
	c4_id = range(len(c4))
	col1 = fits.Column(name='halo_id',format='J',array=c4_id)
	col2 = fits.Column(name='ra_c4',format='D',array=c4['ra_c4'])
	col3 = fits.Column(name='dec_c4',format='D',array=c4['dec_c4'])
	col4 = fits.Column(name='z_c4',format='D',array=c4['z_c4'])
	col5 = fits.Column(name='ra_bcg',format='D',array=c4['ra_bcg'])
	col6 = fits.Column(name='dec_bcg',format='D',array=c4['dec_bcg'])
	col7 = fits.Column(name='z_bcg',format='D',array=c4['z_bcg'])
	col8 = fits.Column(name='z_mean',format='D',array=c4['z_mean'])
	col9 = fits.Column(name='rvir',format='D',array=c4['rho200_rv1500'])
	col10 = fits.Column(name='nc4',format='D',array=c4['nc4'])
	col11 = fits.Column(name='lum',format='D',array=c4['lum'])
	col12 = fits.Column(name='single',format='I',array=np.array(c4['single'],int).ravel())
	col13 = fits.Column(name='keep',format='I',array=c4['keep'])
	col14 = fits.Column(name='sub',format='I',array=c4['sub'])

	cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(data_set+'C4_Halos.fits',clobber=True)
	raise NameError



## Load Previously Sorted Cluster File
c4 = fits.open(data_set+'C4_Halos.fits')[1].data
# Get Avg Halo Data
avg_haloid,avg_ra,avg_dec,avg_z = np.loadtxt('Avg_Halo_Centers.tab',delimiter='\t',unpack=True)
avg_haloid = np.array(avg_haloid,int)


## Load Previously Cut Cluster Sample
load_cut_sample = True
if load_cut_sample == True:
	data_loc = 'individual'
	print ''
	print '...Loading in Previously Cut C4 Cluster Sample'
	f = open(root+'/C4/'+data_loc+'/halo_arrays.pkl','rb')
	input = pkl.Unpickler(f)
	run_dict = input.load()
	globals().update(run_dict)
	f.close()


#blob1 = (np.log10(N200)>2.55)&(np.log10(M200)>14.8)&(np.log10(M200)<15.0)
#blob2 = (np.log10(N200)>2.4)&(np.log10(M200)>15.2)&(np.log10(M200)<15.4)
#blob2[np.where(blob2==True)][9] = False
#blob2[np.where(blob2==True)][11] = False
#cut_out = ~(blob1+blob2)


# Get Cluster Cutouts
cluster_cut = False
if cluster_cut == True:
	arcs = np.array(Cosmo.arcsec_per_kpc_proper(avg_z))*15000. / 3600.    # 15 Mpc in degrees
	for i in range(len(c4['halo_id'])):
		if np.isnan(avg_z[i]) == True: continue
		if i%500 == 0: print i

		# Clus Data
		#clus_z = c4['z_mean'][i]
		#clus_ra = c4['ra_c4'][i]
		#clus_dec = c4['dec_c4'][i]

		clus_ra = avg_ra[i]
		clus_dec = avg_dec[i]
		clus_z = avg_z[i]

		# Get RA and Dec ranges
		d_dec = arcs[i]					# 15 Mpc in Declination 
		d_ra = arcs[i]/np.cos(clus_dec*np.pi/180)	# 15 Mpc in Declination
		d_z = 15000 / c * (1+clus_z)			# 15000 km/s in z

		# Create Galaxy Data Cube out of DR12 Photometric Data
		write_photo = True
		if write_photo == True:
			cut = np.where((photo['ra']<clus_ra+d_ra)&(photo['ra']>clus_ra-d_ra)&(photo['dec']<clus_dec+d_dec)&(photo['dec']>clus_dec-d_dec)&(photo['photoz']<clus_z+0.1)&(photo['photoz']>clus_z-0.1))[0]
			id = photo['objid'][cut]
			ra = photo['ra'][cut]
			dec = photo['dec'][cut]
			z = photo['z'][cut]
			photoz = photo['photoz'][cut]
			appmag_u = photo['petromag_u'][cut]
			appmag_g = photo['petromag_g'][cut]
			appmag_r = photo['petromag_r'][cut]
			appmag_i = photo['petromag_i'][cut]
			appmag_z = photo['petromag_z'][cut]
			appmag_uerr = photo['Petromag_uerr'][cut]
			appmag_gerr = photo['Petromag_gerr'][cut]
			appmag_rerr = photo['Petromag_rerr'][cut]
			appmag_ierr = photo['Petromag_ierr'][cut]
			appmag_zerr = photo['Petromag_zerr'][cut]
			absmag_u = photo['absmag_u'][cut]
			absmag_g = photo['absmag_g'][cut]
			absmag_r = photo['absmag_r'][cut]
			absmag_i = photo['absmag_i'][cut]
			absmag_z = photo['absmag_z'][cut]
			distmod = photo['DistMod'][cut]
			kcorr_u = photo['kcorr_u'][cut]
			kcorr_g = photo['kcorr_g'][cut]
			kcorr_r = photo['kcorr_r'][cut]
			kcorr_i = photo['kcorr_i'][cut]
			kcorr_z = photo['kcorr_z'][cut]

			gal_data = np.vstack([ra,dec,z,photoz,appmag_u,appmag_g,appmag_r,appmag_i,appmag_z,appmag_uerr,appmag_gerr,appmag_rerr,appmag_ierr,appmag_zerr,absmag_u,absmag_g,absmag_r,absmag_i,absmag_z,distmod,kcorr_u,kcorr_g,kcorr_r,kcorr_i,kcorr_z])

			# Project using Photoz
			#ang_d,lum_d = C.zdistance(clus_z,H0)
			#angles = C.findangle(ra,dec,clus_ra,clus_dec)	
			#rdata = angles * ang_d
			#vdata = c * (photoz - clus_z) / (1 + clus_z)
			# RV Cut
			#cut = np.where((vdata < 8000) & (vdata > -8000) & (rdata < 15))[0]
			#id = id[cut]
			#gal_data = gal_data.T[cut].T
			#ra,dec,z,photoz,appmag_u,appmag_g,appmag_r,appmag_i,appmag_z,appmag_uerr,appmag_gerr,appmag_rerr,appmag_ierr,appmag_zerr,absmag_u,absmag_g,absmag_r,absmag_i,absmag_z,distmod,kcorr_u,kcorr_g,kcorr_r,kcorr_i,kcorr_z = gal_data

			# If len(cut) == 0 continue
			#if len(cut) == 0: continue

			# Make Final Data Array with DR12_SpecObjId and DR12_Photoz data!
			col1 = fits.Column(name='objid',format='J',array=id)
			col2 = fits.Column(name='ra',format='D',array=ra)
			col3 = fits.Column(name='dec',format='D',array=dec)
			col4 = fits.Column(name='z',format='D',array=z)
			col5 = fits.Column(name='photoz',format='D',array=photoz)
			col6 = fits.Column(name='appmag_u',format='D',array=appmag_u)
			col7 = fits.Column(name='appmag_g',format='D',array=appmag_g)
			col8 = fits.Column(name='appmag_r',format='D',array=appmag_r)
			col9 = fits.Column(name='appmag_i',format='D',array=appmag_i)
			col10 = fits.Column(name='appmag_z',format='D',array=appmag_z)
			col11 = fits.Column(name='appmag_uerr',format='D',array=appmag_uerr)
			col12 = fits.Column(name='appmag_gerr',format='D',array=appmag_gerr)
			col13 = fits.Column(name='appmag_rerr',format='D',array=appmag_rerr)
			col14 = fits.Column(name='appmag_ierr',format='D',array=appmag_ierr)
			col15 = fits.Column(name='appmag_zerr',format='D',array=appmag_zerr)
			col16 = fits.Column(name='absmag_u',format='D',array=absmag_u)
			col17 = fits.Column(name='absmag_g',format='D',array=absmag_g)
			col18 = fits.Column(name='absmag_r',format='D',array=absmag_r)
			col19 = fits.Column(name='absmag_i',format='D',array=absmag_i)
			col20 = fits.Column(name='absmag_z',format='D',array=absmag_z)
			col21 = fits.Column(name='distmod',format='D',array=distmod)
			col22 = fits.Column(name='kcorr_u',format='D',array=kcorr_u)
			col23 = fits.Column(name='kcorr_g',format='D',array=kcorr_g)
			col24 = fits.Column(name='kcorr_r',format='D',array=kcorr_r)
			col25 = fits.Column(name='kcorr_i',format='D',array=kcorr_i)
			col26 = fits.Column(name='kcorr_z',format='D',array=kcorr_z)

			cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22,col23,col24,col25,col26])
			tbhdu = fits.BinTableHDU.from_columns(cols)
			tbhdu.writeto(data_set+'Halo_'+str(c4['halo_id'][i])+'_PhotoData.fits',clobber=True)


		# Get Galaxy Cutout around Cluster for Spectroscopic Data from DR12
		write_spectro = True
		if write_spectro == True:
			cut = np.where((dr12['ra']<clus_ra+d_ra)&(dr12['ra']>clus_ra-d_ra)&(dr12['dec']<clus_dec+d_dec)&(dr12['dec']>clus_dec-d_dec)&(dr12['z']<clus_z+d_z)&(dr12['z']>clus_z-d_z))[0]
			gal_specobjid = dr12['specobjid'][cut]
			gal_objid = dr12['objid'][cut]
			gal_ra = dr12['ra'][cut]
			gal_dec = dr12['dec'][cut]
			gal_z = dr12['z'][cut]
			gal_zerr = dr12['zerr'][cut]
			gal_class = dr12['class'][cut]
			gal_subclass = dr12['subclass'][cut]
			gal_DistMod = dr12['DistMod'][cut]
			gal_absMagU = dr12['absMagU'][cut]
			gal_absMagG = dr12['absMagG'][cut]
			gal_absMagR = dr12['absMagR'][cut]
			gal_absMagI = dr12['absMagI'][cut]
			gal_absMagZ = dr12['absMagZ'][cut]
			gal_kcorrU = dr12['kcorrU'][cut]
			gal_kcorrG = dr12['kcorrG'][cut]
			gal_kcorrR = dr12['kcorrR'][cut]
			gal_kcorrI = dr12['kcorrI'][cut]
			gal_kcorrZ = dr12['kcorrZ'][cut]


			# Project	
			ang_d,lum_d = C.zdistance(clus_z,H0)
			angles = C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)	
			rdata = angles * ang_d
			vdata = c * (gal_z - clus_z) / (1 + clus_z)
			# RV Cut
			cut = np.where((vdata < 7000) & (vdata > -7000) & (rdata < 12))[0]
			gal_specobjid,gal_objid = gal_specobjid[cut],gal_objid[cut]
			gal_class,gal_subclass = gal_class[cut],gal_subclass[cut]
			gal_ra,gal_dec,gal_z,gal_zerr = gal_ra[cut],gal_dec[cut],gal_z[cut],gal_zerr[cut]
			gal_DistMod,gal_absMagU,gal_absMagG,gal_absMagR,gal_absMagI,gal_absMagZ = gal_DistMod[cut],gal_absMagU[cut],gal_absMagG[cut],gal_absMagR[cut],gal_absMagI[cut],gal_absMagZ[cut]
			gal_kcorrU,gal_kcorrG,gal_kcorrR,gal_kcorrI,gal_kcorrZ = gal_kcorrU[cut],gal_kcorrG[cut],gal_kcorrR[cut],gal_kcorrI[cut],gal_kcorrZ[cut]

			# Make FITS table
			col1 = fits.Column(name='specobjid',format='J',array=gal_specobjid)
			col2 = fits.Column(name='objid',format='J',array=gal_objid)
			col3 = fits.Column(name='ra',format='D',array=gal_ra)
			col4 = fits.Column(name='dec',format='D',array=gal_dec)
			col5 = fits.Column(name='z',format='D',array=gal_z)
			col6 = fits.Column(name='zerr',format='D',array=gal_zerr)
			col7 = fits.Column(name='class',format='A',array=gal_class)
			col8 = fits.Column(name='subclass',format='A',array=gal_subclass)
			col9 = fits.Column(name='absMagU',format='D',array=gal_absMagU)
			col10 = fits.Column(name='absMagG',format='D',array=gal_absMagG)
			col11 = fits.Column(name='absMagR',format='D',array=gal_absMagR)
			col12 = fits.Column(name='absMagI',format='D',array=gal_absMagI)
			col13 = fits.Column(name='absMagZ',format='D',array=gal_absMagZ)
			col14 = fits.Column(name='kcorrU',format='D',array=gal_kcorrU)
			col15 = fits.Column(name='kcorrG',format='D',array=gal_kcorrG)
			col16 = fits.Column(name='kcorrR',format='D',array=gal_kcorrR)
			col17 = fits.Column(name='kcorrI',format='D',array=gal_kcorrI)
			col18 = fits.Column(name='kcorrZ',format='D',array=gal_kcorrZ)
			col19 = fits.Column(name='DistMod',format='D',array=gal_DistMod)
			cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19])
			tbhdu = fits.BinTableHDU.from_columns(cols)
			tbhdu.writeto(data_set+'Halo_'+str(c4['halo_id'][i])+'_SpectroData.fits',clobber=True)


# Calculate Nspec per cluster
calc_Nspec = False
if calc_Nspec == True:
	f = open('individual/halo_arrays.pkl','rb')
	input = pkl.Unpickler(f)
	run_dict = input.load()
	globals().update(run_dict)
	f.close()

	rlim = 1.0	# R_vir
	Vlim = 5000	# km/s
	Nspec = []
	Index = []
	for i in range(len(HaloID)):
		if i%100 == 0: print i
		# Define Rlim
		Rlim = RVIR[i]*rlim
		# Get Halo Data
		clus_ra = RA[i]
		clus_dec = DEC[i]
		clus_z = Z[i]
		try:
			data = fits.open(data_set+'Halo_'+str(HaloID[i])+'_SpectroData.fits')[1].data
			gal_ra,gal_dec,gal_z = data['ra'],data['dec'],data['z']
			angles = C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
			rdata = angles * C.zdistance(clus_z,H0)[0]
			vdata = c * (gal_z - clus_z) / (1 + clus_z)

			cut = np.where((rdata<Rlim)&(np.abs(vdata)<Vlim))[0]
			Nspec.append(cut.size)
			Index.append(i)

		except IOError:
			Nspec.append(-1)
			Index.append(i)
	Nspec = np.array(Nspec)
	Index = np.array(Index)


# Calculate Richnesses
class RICHNESS(object):

	def __init__(self,varib):
		self.__dict__.update(varib)

	def richness_est(self,ra,dec,z,pz,gmags,rmags,imags,abs_rmags,haloid,clus_ra,clus_dec,clus_z,clus_rvir=None,shiftgap=True,use_specs=True,fit_rs=False):
		''' Richness estimator, magnitudes should be apparent mags except for abs_rmag
			z : spectroscopic redshift
			pz : photometric redshift
		'''
		# Put Into Class Namespace
		keys = ['ra','dec','z','pz','gmags','rmags','imags','clus_ra','clus_dec','clus_z','clus_rvir','vel_disp','color_data','color_cut','RS_color','RS_sigma','clus_color_cut','signal','background','richness','mems','phot','spec','Nspec','deg_per_rvir','SA_outer','SA_inner','outer_red_dense','inner_background','Nphot_inner','Nphot_outer','red_inner','red_outer','all','outer_edge','inner_edge','shift_cut','set']
		self.__dict__.update(ez.create(keys,locals()))

		# Take rough virial radius measurement, this method is BAD...
		if clus_rvir == None:
			clus_rvir = np.exp(-1.86)*len(np.where((rmags < -19.55) & (rdata < 1.0) & (np.abs(vdata) < 3500))[0])**0.51
			self.clus_rvir = clus_rvir

		# Magnitue Cut at -19 Absolute R Mag and sort by them
		bright = np.where(abs_rmags < -19.5)[0]
		data = np.vstack([ra,dec,z,pz,gmags,rmags,imags,abs_rmags])
		data = data.T[bright].T
		ra,dec,z,pz,gmags,rmags,imags,abs_rmags = data
                data = np.vstack([ra,dec,z,pz,gmags,rmags,imags,abs_rmags])
		sorts = np.argsort(abs_rmags)
		data = data.T[sorts].T
                ra,dec,z,pz,gmags,rmags,imags,abs_rmags = data
		self.__dict__.update(ez.create(keys,locals()))

		# Separate Into Spectro Z and Photo Z DataSets
		spec_cut = np.abs(z) > 1e-4
		phot_cut = np.where(np.abs(z) < 1e-4)[0]
		all = {'ra':ra,'dec':dec,'z':z,'pz':pz,'gmags':gmags,'rmags':rmags,'imags':imags,'abs_rmags':abs_rmags}
		spec = {'ra':ra[spec_cut],'dec':dec[spec_cut],'z':z[spec_cut],'pz':pz[spec_cut],'gmags':gmags[spec_cut],'rmags':rmags[spec_cut],'imags':imags[spec_cut],'abs_rmags':abs_rmags[spec_cut]}
		phot = {'ra':ra[~spec_cut],'dec':dec[~spec_cut],'z':z[~spec_cut],'pz':pz[~spec_cut],'gmags':gmags[~spec_cut],'rmags':rmags[~spec_cut],'imags':imags[~spec_cut],'abs_rmags':abs_rmags[~spec_cut]}

		# Project Spectra into radius and velocity
		ang_d,lum_d = C.zdistance(clus_z,self.H0)
		angles = C.findangle(spec['ra'],spec['dec'],clus_ra,clus_dec)
		rdata = angles * ang_d
		vdata = self.c * (spec['z'] - clus_z) / (1 + clus_z)

		# Shiftgapper for Interlopers
		if shiftgap == True:
			len_before = np.where((rdata < clus_rvir*1.5)&(np.abs(rdata)<4000))[0].size
			clus_data = np.vstack([rdata,vdata,spec['ra'],spec['dec'],spec['z'],spec['pz'],spec['gmags'],spec['rmags'],spec['imags'],spec['abs_rmags']])
			clus_data = C.shiftgapper(clus_data.T).T
			sorts = np.argsort(clus_data[-1])
			clus_data = clus_data.T[sorts].T
			rdata,vdata,spec['ra'],spec['dec'],spec['z'],spec['pz'],spec['gmags'],spec['rmags'],spec['imags'],spec['abs_rmags'] = clus_data
			shift_cut = len_before - np.where((rdata < clus_rvir)&(np.abs(rdata)<4000))[0].size

		# Measure Velocity Dispersion of all galaxies within 1 * r_vir and np.abs(vdata) < 4000 km/s
		vel_disp = astats.biweight_midvariance(vdata[np.where((rdata < clus_rvir*1)&(np.abs(vdata) < 4000))])
		spec.update({'rdata':rdata,'vdata':vdata})
		self.__dict__.update(ez.create(keys,locals()))


		# Get Members from Specta, get their red sequence color
		mems = np.where((spec['rdata'] < clus_rvir)&(np.abs(spec['vdata'])<2*vel_disp))[0]
		color_data = spec['gmags'] - spec['rmags']
		color_cut = np.where((color_data[mems] < 1.15) & (color_data[mems] > 0.65))[0]
		RS_color = astats.biweight_location(color_data[mems][color_cut])
		RS_sigma = astats.biweight_midvariance(color_data[mems][color_cut])
		if fit_rs == True:
			clf = linear_model.LinearRegression()
			set = np.where(np.abs(color_data[mems]-RS_color)<3*RS_sigma)[0]
			clf.fit(spec['rmags'][mems][set].reshape(set.size,1),color_data[mems][set])
			
		# Nspec is # of members that are within 2 sigma of cluster color
		clus_color_cut = np.where(np.abs(color_data[mems] - RS_color) < RS_sigma*2)[0]
		Nspec = len(clus_color_cut)
		self.__dict__.update(ez.create(keys,locals()))

		# Plot
		make_plot = False
		if make_plot == True:
			fig,ax = mp.subplots()
			ax.plot(rmags,gmags-rmags,'ko',alpha=.8)
			ax.plot(spec['rmags'][mems],color_data[mems],'co')
			ax.axhline(RS_color,color='r')	
			ax.axhline(RS_color+RS_sigma,color='b')
			ax.axhline(RS_color-RS_sigma,color='b')
			ax.set_xlim(13,19)
			ax.set_ylim(0,1.3)
			ax.set_xlabel('Apparent R Mag',fontsize=16)
			ax.set_ylabel('App G Mag - App R Mag',fontsize=16)
			ax.set_title('Color-Mag Diagram, Cluster '+str(haloid))
			fig.savefig('Halo_'+str(haloid)+'_color_mag.png',bbox_inches='tight')
			mp.close(fig)

		# Get Rdata from PhotoZ Data Set
		angles = C.findangle(all['ra'],all['dec'],clus_ra,clus_dec)
		all['rdata'] = angles * ang_d

		# Get arcsec per kpc proper
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
		red_inner = np.where(((all['gmags']-all['rmags']) < RS_color + 1.5*RS_sigma)&((all['gmags']-all['rmags']) > RS_color - 1.5*RS_sigma)&(all['rdata']<clus_rvir))[0]
		red_outer = np.where(((all['gmags']-all['rmags']) < RS_color + 1.5*RS_sigma)&((all['gmags']-all['rmags']) > RS_color - 1.5*RS_sigma)&(all['rdata']<outer_edge*clus_rvir)&(all['rdata']>inner_edge*clus_rvir))[0]
		Nphot_inner = len(red_inner)
		Nphot_outer = len(red_outer)

		self.__dict__.update(ez.create(keys,locals()))

		# Get Solid Angle Density of Outer Red Galaxies
		outer_red_dense = Nphot_outer / SA_outer
		inner_background = int(np.ceil(outer_red_dense * SA_inner))

		# If inner_background is less than Nphot_inner, then Nphot_inner -= inner_background, otherwise Nphot_inner = 0
		if inner_background < Nphot_inner:
			Nphot_inner -= inner_background
		else:
			Nphot_inner = 0

		# Richness = Nspec + Nphot_inner or just Nphot_inner
		if use_specs == True:
			richness = Nspec + Nphot_inner
		else:
			richness = Nphot_inner
		self.__dict__.update(ez.create(keys,locals()))

		return richness

R = RICHNESS(varib)
calc_rich = True
if calc_rich == True:

	ra_lim = (123,250)
	dec_lim = (-4,66)
	arcs = np.array(Cosmo.arcsec_per_kpc_proper(Z))*15000. / 3600.    # 15 Mpc in degrees
	d_dec = np.copy(arcs)
	d_ra = arcs/np.cos(DEC*np.pi/180)
        danger = np.where( (RA+d_ra > ra_lim[1]) | (RA-d_ra < ra_lim[0]) | (DEC+d_dec > dec_lim[1]) | (DEC-d_dec < dec_lim[0]))[0]

	use_specs = False
	clus_id = []
	richness  = []
	nspec = []
	vel_disp = []
	for i in range(len(HaloID)):	#len(c4['halo_id'])):
		if i%400 == 0: print i

#		if c4['sub'][i]>3 or np.isnan(avg_ra)[i] == True: continue

		# Load Halo Data
		#haloid = c4['halo_id'][i]
		#clus_ra = c4['ra_bcg'][i]
		#clus_dec = c4['dec_bcg'][i]
		#clus_z = c4['z_mean'][i]
		#clus_rvir = c4['rvir'][i]
		#if clus_rvir < 0.5: clus_rvir = 0.5
		#elif clus_rvir > 2: clus_rvir = 2

		# Use Avg Halo Data
		#use_avg = True
		#if use_avg == True and ~np.isnan(avg_ra[i]):
		#	clus_ra = avg_ra[i]
		#	clus_dec = avg_dec[i]
		#	clus_z = avg_z[i]

		haloid = HaloID[i]
		clus_ra = RA[i]
		clus_dec = DEC[i]
		clus_z = Z[i]
		clus_rvir = RVIR[i]
		clus_hvd = HVD[i]

		# Get Other nearby Clusters
		nearby = np.where((avg_ra < clus_ra + 0.5)&(avg_ra > clus_ra - 0.5)&(avg_dec < clus_dec + 0.5)&(avg_dec > clus_dec - 0.5)&(avg_ra != clus_ra))[0]
		near_ra = avg_ra[nearby]
		near_dec = avg_dec[nearby]
		near_z = avg_z[nearby]

		# Load Galaxy Data
		galdata = fits.open(data_set+'Halo_'+str(HaloID[i])+'_PhotoData.fits')[1].data
		gal_ra = galdata['ra']
		gal_dec = galdata['dec']
		gal_z = galdata['z']
		gal_photoz = galdata['photoz']
		gal_absg = galdata['absmag_g']
		gal_absr = galdata['absmag_r']
		gal_absi = galdata['absmag_i']
		gal_appg = galdata['appmag_g']
		gal_appr = galdata['appmag_r']
		gal_appi = galdata['appmag_i']

		rich = R.richness_est(gal_ra,gal_dec,gal_z,gal_photoz,gal_appg,gal_appr,gal_appi,gal_absr,haloid,clus_ra,clus_dec,clus_z,clus_rvir=clus_rvir,use_specs=use_specs)
		richness.append(rich)
		clus_id.append(haloid)
		vel_disp.append(R.vel_disp)
		nspec.append(R.Nspec)

	richness = np.array(richness)
	clus_id = np.array(clus_id)
	vel_disp = np.array(vel_disp)
	nspec = np.array(nspec)

	write = True
	if write == True:
		f = open('C4_Richnesses.tab','w')
		f.write('# Richnesses and Velocity Dispersions for C4 Clusters\n')
		f.write('# HaloID, Richness, Vel_Disp, Nspec\n')
		for i in range(len(richness)):
			f.write(str(clus_id[i])+'\t'+str(richness[i])+'\t'+str(vel_disp[i])+'\t'+str(nspec[i])+'\n')
		f.close()


## Positional Recentering
recenter = False
if recenter == True:

	AVG_RA,AVG_DEC,AVG_Z = [],[],[]

	for i in range(len(c4['halo_id'])):
		if i%500 == 0: print i
		# Load Halo Data
		haloid = c4['halo_id'][i]
		clus_ra = c4['ra_c4'][i]
		clus_dec = c4['dec_c4'][i]
		clus_z = c4['z_mean'][i]
		clus_rvir = c4['rvir'][i]
		if clus_rvir < 0.5:
			clus_rvir = 0.5
		elif clus_rvir > 2:
			clus_rvir = 2

		# Load Galaxy Data
		galdata = fits.open(data_set+'Halo_'+str(i)+'_PhotoData.fits')[1].data
		gal_ra = galdata['ra']
		gal_dec = galdata['dec']
		gal_z = galdata['z']
		gal_photoz = galdata['photoz']
		gal_absg = galdata['absmag_g']
		gal_absr = galdata['absmag_r']
		gal_absi = galdata['absmag_i']
		gal_appg = galdata['appmag_g']
		gal_appr = galdata['appmag_r']
		gal_appi = galdata['appmag_i']
		gal_data = np.vstack([gal_ra,gal_dec,gal_z,gal_photoz,gal_absg,gal_absr,gal_absi,gal_appg,gal_appr,gal_appi])

		# Absolute Magnitude Cut, abs_rmag < -20, -19.5, -19, and gal_z != 0
		cut = np.where((np.abs(gal_z)>1e-4)&(gal_absr < -19))
		gal_data = gal_data.T[cut].T
		gal_ra,gal_dec,gal_z,gal_photoz,gal_absg,gal_absr,gal_absi,gal_appg,gal_appr,gal_appi = gal_data

		# Project
		ang_d,lum_d = C.zdistance(clus_z,70)
		angles = C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
		rdata = angles * ang_d
		vdata = c * (gal_z - clus_z) / (1 + clus_z)

		# Select Member Gals and get Avg RA, DEC, Z, and new rdata, vdata
		select = np.where((rdata < clus_rvir) & (np.abs(vdata) < 1500))[0]
		avg_ra = astats.biweight_location(gal_ra[select])
		avg_dec = astats.biweight_location(gal_dec[select])
		avg_z = astats.biweight_location(gal_z[select])
		ang_d,lum_d = C.zdistance(avg_z,70)
		angles = C.findangle(gal_ra,gal_dec,avg_ra,avg_dec)
		rdata2 = angles * ang_d
		vdata2 = c * (gal_z - avg_z) / (1 + avg_z)

		AVG_RA.append(avg_ra)
		AVG_DEC.append(avg_dec)
		AVG_Z.append(avg_z)

	AVG_RA = np.array(AVG_RA)
	AVG_DEC = np.array(AVG_DEC)
	AVG_Z = np.array(AVG_Z)

	# Write File
	write = True
	if write == True:
		f = open('Avg_Halo_Centers.tab','w')
		f.write('# Average (robust median) Halo Centers in RA, DEC and Z using Rdata < 2.0 and abs(Vdata) < 2000\n')
		f.write('# HaloID, avg_ra, avg_dec, avg_z\n')
		for i in range(len(c4['halo_id'])):
			f.write(str(c4['halo_id'][i])+'\t'+str(AVG_RA[i])+'\t'+str(AVG_DEC[i])+'\t'+str(AVG_Z[i])+'\n')
		f.close()





