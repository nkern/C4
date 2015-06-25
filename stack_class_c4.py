## stack_class_c4.py
"""
This file contains classes and functions necessary in stacking C4 Data
halos and running the caustic technique over them, utilizing both 'causticpy' and 'caustic_stack'
programs.
"""

from caustic_stack import *

class CFOUR(object):
	"""
	A class with functions designed to work with C4 data
	"""

	def __init__(self,varib):
		self.__dict__.update(varib)
		self.U = Universal(varib)


	def load_halos(self,sort=True,avg_centers=True):
		'''This function loads halo data and makes cosmology corrections'''

		# Halo Data
		c4 = fits.open(self.root+'/C4/'+self.data_set+'C4_Halos.fits')[1].data
		RA = c4['ra_bcg']
		DEC = c4['dec_bcg']
		Z = c4['z_mean']
		RVIR = c4['rvir']
		SINGLE = c4['single']
		SUB = c4['sub']
		NC4 = c4['nc4']
	
		# Get Averaged Halo Centers
		if avg_centers == True:
			avg_HaloID,RA,DEC,Z = np.loadtxt(self.root+'/C4/'+'Avg_Halo_Centers.tab',delimiter='\t',unpack=True)
			avg_HaloID = np.array(avg_HaloID,int)

		# Richnesses based on Photo and Spectro Z
		clus_id,N200,HVD,Nspec = np.loadtxt(self.root+'/C4/'+'C4_Richnesses.tab',delimiter='\t',unpack=True)
		index = []
		for i in range(len(clus_id)):
			if clus_id[i] in avg_HaloID: index.append(np.where(avg_HaloID==clus_id[i])[0][0])
		index = np.array(index)

		RA,DEC,Z,RVIR,SINGLE,SUB,NC4 = RA[index],DEC[index],Z[index],RVIR[index],SINGLE[index],SUB[index],NC4[index]
		HaloData = np.vstack([RA,DEC,Z,RVIR,HVD,N200,SINGLE,SUB,NC4])
		HaloID = avg_HaloID[index]

		if sort == True:
			HaloID, HaloData = self.sort_halos(HaloID,N200,HaloData)
			RA,DEC,Z,RVIR,HVD,N200,SINGLE,SUB,NC4 = HaloData

		# Load M200 and HVD's from individual runs if available
		try:
			ID,M200,HVD2,R200 = np.loadtxt(self.root+'/C4/individual/Iteration1_ClusData.tab',delimiter='\t',unpack=True,usecols=(0,1,4,3))
			if len(M200) != len(HaloID):
				match = np.array(map(lambda x: np.where(ID == x)[0],HaloID))
				match = match[~np.isnan(match)]
				M200,HVD = M200[match],HVD2[match]
		except:
			M200 = np.zeros(len(HaloID))
			# Fix Certain Variables
			RVIR[np.where(RVIR > 2.5)] = 2.5
			HVD[np.where(HVD > 1000)] = 1000

		HaloData = np.vstack([RA,DEC,Z,M200,RVIR,HVD,N200,SINGLE,SUB,NC4,Nspec])

		return HaloID, HaloData


	def sort_halos(self,HaloID,Observable,HaloData):
		''' Sort Halo Data by some Criteria '''
		# Sort Arrays by descending Observable
		sort = np.argsort(Observable)[::-1]	
		HaloID = HaloID[sort]
		HaloData = HaloData.T[sort].T
	
		# Return packed array
		return HaloID, HaloData 


	def load_galaxies(self,haloid):
		''' Loads haloid galaxies from a local directory '''
		galdata = fits.open(self.root+'/C4/'+self.data_set+'Halo_'+str(haloid)+'_SpectroData.fits')[1].data
		gal_ra = galdata['ra']
		gal_dec = galdata['dec']
		gal_z = galdata['z']
		gal_gmags = galdata['absMagG']
		gal_rmags = galdata['absMagR']
		gal_imags = galdata['absMagI']

		return np.vstack([gal_ra,gal_dec,gal_z,gal_gmags,gal_rmags,gal_imags])


	def load_chris_gals(self,haloid):
		''' Loads haloid galaxies from Chris' data '''
		galra,galdec,gal_umag,gal_gmag, gal_rmag, gal_imag,gal_zmag, gal_z = np.loadtxt(self.chris_data_root+'/proj_data/C4_dr12_'+str(haloid)+'_data.csv',dtype='float',delimiter=',',skiprows=1,usecols=(1,2,3,4,5,6,7,13),unpack=True)

		return np.vstack([galra,galdec,gal_z,gal_gmag,gal_rmag,gal_imag])


	def load_project_append(self,haloid,clus_ra,clus_dec,clus_z,clus_m200,clus_r200,clus_hvd,PS):
		"""
		This function loads galaxy data then projects it, discarding the galaxy data after using it to save memory
		"""

		#Load Galaxies
		if self.chris_data == False:
			galdata = self.load_galaxies(haloid)
			gal_ra,gal_dec,gal_z,gal_gmags,gal_rmags,gal_imags = galdata
		else:
			galdata = self.load_chris_gals(haloid)
                        gal_ra,gal_dec,gal_z,gal_gmags,gal_rmags,gal_imags = galdata

		# Get angles and phase spaces
		ang_d,lum_d = self.U.C.zdistance(clus_z,self.H0)
		angles = self.U.C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
		rdata = angles * ang_d
		vdata = self.c * (gal_z - clus_z) / (1 + clus_z)

		# Append to PS
		PS.append( {'Rdata':rdata,'Vdata':vdata,'G_Mags':gal_gmags,'R_Mags':gal_rmags,'I_Mags':gal_imags,'HaloID':haloid,'R200':clus_r200,'HVD':clus_hvd,'M200':clus_m200,'Z':clus_z} )


        def load_project_append_bootstrap(self,HaloID,M200,R200,HVD,Z,Halo_P,Halo_V,PS,weight):
                """
                This function loads galaxy data then projects it, discarding the galaxy data after using it to save memory
                """

                # Load Galaxies
		if self.lightcone == True:
			Gal_RA,Gal_DEC,Gal_Z,G_Mags,R_Mags,I_Mags,Gal_P,Gal_V = self.configure_galaxies(HaloID,np.array([M200,R200,HVD,Z,Halo_P[0],Halo_P[1],Halo_P[2],Halo_V[0],Halo_V[1],Halo_V[2]]))
		else:
                	Gal_P,Gal_V,G_Mags,R_Mags,I_Mags = self.configure_galaxies(HaloID,np.array([M200,R200,HVD,Z,Halo_P[0],Halo_P[1],Halo_P[2],Halo_V[0],Halo_V[1],Halo_V[2]]))

                # Do Projection
                r, v, pro_pos = self.U.line_of_sight(Gal_P,Gal_V,Halo_P,Halo_V,project=False)

                r = np.array(r)
                v = np.array(v)
                pro_pos = np.array(pro_pos)
                G_Mags = np.array(G_Mags)
                R_Mags = np.array(R_Mags)
                I_Mags = np.array(I_Mags)

		for m in range(weight):
                	# Append to PS
                	PS.append( {'Rdata':r,'Vdata':v,'pro_pos':np.array(pro_pos),'G_Mags':G_Mags,'R_Mags':R_Mags,'I_Mags':I_Mags,'HaloID':HaloID,'M200':M200,'R200':R200,'HVD':HVD} )


#	Outdated Mass Mixing Routine
#	def mass_mixing(self,HaloID,HaloData,mass_scat):
#		'''
#		This function performs a mass mixing procedure with a given fractional scatter in assumed mass
#		'''
#		# Unpack Array
#		M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData
#		# Create lognormal distribution about 1 with width mass_scat, length HaloID.size
#		mass_mix = npr.lognormal(0,mass_scat,len(HaloID))
#		# Apply Mass Scatter
#		M_crit200 *= mass_mix
#		# Create M200_match array
#		M_crit200_match = np.copy(M_crit200)
#		# Sort by Descending Mass
#		sort = np.argsort(M_crit200)[::-1]
#		M_crit200 = M_crit200[sort]
#		R_crit200 = R_crit200[sort]
#		Z = Z[sort]
#		HVD = HVD[sort]
#		HPX = HPX[sort]
#		HPY = HPY[sort]
#		HPZ = HPZ[sort]
#		HVX = HVX[sort]
#		HVY = HVY[sort]
#		HVZ = HVZ[sort]
#		HaloID = HaloID[sort]
#		# Re-pack
#		HaloData = np.array([M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ])
#		return HaloID,HaloData,sort


#	Outdated Richness Estimator, see c4_richness.py for new one
#	def richness_est(self,rdata,vdata,gmags,rmags,imags,shiftgap=True):
#		''' Richness estimator, magnitudes should be absolute mags'''
#
#		# Rough Cut at < 7 Mpc and +/- 5000 km/s
#		rmag_limit = -19.0
#		cut = np.where((rdata < 7) & (np.abs(vdata) < 5000) & (rmags < rmag_limit))[0]
#		rdata = rdata[cut]
#		vdata = vdata[cut]
#		gmags = gmags[cut]
#		rmags = rmags[cut]
#		imags = imags[cut]
#
#		# Shiftgapper for Interlopers
#		if shiftgap == True:
#			clus_data = np.vstack([rdata,vdata,gmags,rmags,imags])
#			clus_data = self.U.C.shiftgapper(clus_data.T).T
#			rdata,vdata,gmags,rmags,imags = clus_data
#
#		# Take rough virial radius measurement
#		r_vir = np.exp(-1.86)*len(np.where((rmags < -19.55) & (rdata < 1.0) & (np.abs(vdata) < 3500))[0])**0.51
#
#		# Measure Velocity Dispersion of all galaxies within 1.25 * r_vir
#		vel_disp = astats.biweight_midvariance(vdata[np.where(rdata < (r_vir*1.25))])
#
#		# Find color of Red Sequence, measured as SDSS_g-SDSS_r vs. SDSS_r absolute magnitude
#		# This is almost always around an abs_gmags-abs_rmags ~ 0.8
#		color_data = gmags-rmags
#		color_cut = np.where((color_data<1.0)&(color_data>0.65))[0]
#		size = len(color_cut)
#		RS_color = astats.biweight_location(color_data[color_cut]) 
#		RS_sigma = astats.biweight_midvariance(color_data[color_cut])
#	
#		# Measure Richness and Background
#		signal = len(np.where((np.abs(vdata) < vel_disp*2)&(rdata <= r_vir))[0])
#		background = len(np.where((np.abs(vdata) < vel_disp*2)&(rdata <= r_vir*6)&(rdata >= r_vir*5)&(color_data<(RS_color+RS_sigma))&(color_data>(RS_color-RS_sigma)))[0])
#
#		richness = signal - background
#		return richness
















