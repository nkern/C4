'''
This file contains run-dependent parameters used by c4_stack.py file

file location : /n/Christoq1/nkern/C4/
'''
import time

### FLAGS ###
self_stack	= False				# Run self-stack or bin-stack
scale_data	= True				# Scale to-be-stacked phase space radial data by r200 if True
lightcone	= True				# If True, working on Henriques lightcone, if False, working on Guo data cube
write_data 	= False				# Write Data to Result directories if True
init_clean	= False				# Do an extra shiftgapper on ensemble before the lines of sight get stacked.
small_set	= False				# 100 Halo Set or 2000 Halo Set
mass_mix	= False				# Incorporate Mass Mixing Models?
bootstrap	= False				# Perform a bootstrapping technique to estimate error in mass estimation?
new_halo_cent	= True				# Use Updated Halo Centers instead of BCG Values
mirror		= True				# Mirror Phase Space in Caustic Surface Estimation?
true_mems	= False				# Run with only gals within r200?
run_los		= False				# Run line of sight mass estimation or not
cent_offset	= None                          # Either 'r', 'v', 'full', or None.

### CONSTANTS ###
# Run Dependent
ens_num         = 1				# Number of Ensembles to build and solve for IN THIS RUN
gal_num         = 25				# Number of galaxies taken per line of sight
line_num        = 20				# Number of lines of sight to stack over
method_num      = 0                             # Ensemble Build Method Number
cell_num        = 0                             # Cell Number ID corresponding to given gal_num & line_num geometry in a Run Table
table_num       = 1                             # Table Re-Run Version  
data_loc        = None				# Alternative data_loc, either None or String
write_loc       = None				# Alternative write_loc, either None or String

# Other Techniques
mm_est		= 'richness'			# "Mass Mix Estimator" - What scaling technique to get mass scatter? 'richness', 'vel_disp', 'luminosity', 
edge_perc	= 0.15				# Fractional percent of Top galaxies used in edge detection technique
center_scat     = None                          # If guessing halo center, fractional induced scatter into known center
avg_meth        = 'median'			# If bin stacking, by which method do you average bin properties? (ex. median, mean)
bootstrap_num   = None				# Highest directory marker for bootstrap data, ex. bootstrap1
bootstrap_rep   = None				# Bootstrap repetition directory marker, ex. bootstrap1/rep1

# Caustic Technique Dependent
q		= 10.0				# Scale of Gaussian Kernel Density Estimator
beta		= 0.2				# Velocity Anisotropy Beta parameter, if constant profile
fbeta		= 0.65				# fbeta value, see 'Diaferio 1999'
r_limit 	= 1.5				# Phase space radius Cut Scaled by R200
v_limit		= 4000.0			# Phase space velocity Cut in km/s

# Data Set
data_set	= 'DataSet1/'			# Data set to draw data from
photo_data	= 
halo_num	= 1000				# Total number of halos to work with

# Physical
c               = 2.99792e5                     # speed of light in km/s
h               = 0.7                           # Hubble Constant, unitless
H0              = h*100.0                       # Hubble Constant, km s-1 Mpc-1

# Other
run_time        = time.asctime()              	# Time when program was started
root            = '/nfs/christoq_ls/nkern'	# Root for OSDC

