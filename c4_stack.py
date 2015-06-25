## c4_stack.py
"""
This program uses the Caustic Technique to estimate the masses of galaxy clusters,
after having applied a stacking technique to create ensemble clusters.
This script does strictly bin stacking.

"""
## Import Modules
print '...importing modules'
from caustic_stack import *
from stack_class_c4 import *
from fits_table import fits_table

sys.path.insert(0,os.getcwd())
__import__('caustic_params')
from caustic_params import *
print "Loaded caustic_params from",sys.modules['caustic_params']

## Fed Positional Parameter ##
run_num	= int(sys.argv[1])		# run_num th iteration of the whole job array in PBS script

## Bin Stacking
stack_range = np.arange(run_num*ens_num*line_num,run_num*ens_num*line_num+ens_num*line_num)
if data_loc == None:
	data_loc = 'binstack/bs_run_table'+str(table_num)
if write_loc == None:
	write_loc = 'bs_m'+str(method_num)+'_run'+str(cell_num)			
if bootstrap == True:
	write_loc = 'bo_m'+str(method_num)+'_run'+str(cell_num)
	data_loc = 'binstack/bootstrap'+str(bootstrap_num)+'/rep'+str(bootstrap_rep)

## Make dictionary for loaded constants, doesn't matter if some don't exist
keys = ['c','h','H0','q','beta','fbeta','r_limit','v_limit','photo_data','spec_data','halo_num','gal_num','line_num','method_num','write_loc','data_loc','root','self_stack','scale_data','write_data','run_time','init_shiftgap','shiftgapper','edge_int_remove','small_set','run_los','bootstrap','run_num','ens_num','cell_num','stack_range','mass_mix','mass_scat','bootstrap_num','bootstrap_rep','avg_meth','cent_offset','center_scat','new_halo_cent','true_mems','mirror','edge_perc','lightcone','mm_est','data_set']

varib = ez.create(keys,locals())
varib.update({'_name_':'varib'})

## INITIALIZATION ##
S = Stack(varib)
U = Universal(varib)
C4 = CFOUR(varib)

U.print_separation('## Running c4_stack.py')
names = ['run_time','','run_num','gal_num','line_num','cell_num','ens_num','halo_num','method_num','avg_meth','','write_data','scale_data','run_los','mirror','new_halo_cent','init_shiftgap','shiftgapper','edge_int_remove','bootstrap','bootstrap_num','bootstrap_rep','','data_set','data_loc','write_loc','root']
U.print_varibs(names,varib)

## Load Halo Data
U.print_separation('# ...Loading Halos',type=1)
HaloID,HaloData = C4.load_halos(sort=True,avg_centers=True)
RA,DEC,Z,M200,RVIR,HVD,N200,SINGLE,SUB,NC4,Nspec = HaloData

# Create Dummy M200
M200 = np.zeros(len(RA))

# Load in Halo Data from previously written file, or write file if doesn't exist
try:
	try:
		f = open(root+'/C4/'+data_loc+'/'+write_loc+'/halo_arrays.pkl','rb')
		input = pkl.Unpickler(f)	
		run_dict = input.load()
		globals().update(run_dict)
		f.close()
	except:
                f = open(root+'/C4/'+data_loc+'/halo_arrays.pkl','rb')
                input = pkl.Unpickler(f)        
                run_dict = input.load()
                globals().update(run_dict)
                f.close()

except IOError:
	try:
		f = open(root+'/C4/'+data_loc+'/halo_arrays.pkl','wb')	
		output = pkl.Pickler(f)
		run_dict = {'HaloID':HaloID,'RA':RA,'DEC':DEC,'Z':Z,'HVD':HVD,'N200':N200,'RVIR':RVIR,'SINGLE':SINGLE,'SUB':SUB,'NC4':NC4,'M200':M200,'Nspec':Nspec}
		output.dump(run_dict)
		f.close()
	except IOError:
		pass

## Selection Function
# 2650 clusters for TRUTH_CASUTIC
select = False
if select == True:
	cut = np.where((N200>10)&(Z>0.05)&(Z<.13)&(Nspec>5))[0]
	for i in run_dict.keys():
		run_dict[i] = run_dict[i][cut]
	globals().update(run_dict)


## Load in Chris' C4 Data
chris_data = True
C4.chris_data = chris_data
C4.chris_data_root = '/nfs/christoq_ls/MILLENNIUM/Henriques/TRUTH_CAUSTIC'
#C4.chris_data_root = '/nfs/christoq_ls/C4/sdssdr12'
load_chris_halos = True
if load_chris_halos == True:
	halos = fits.open(root+'/C4/'+data_loc+'/halos.fits')[1].data
	HaloID = halos['orig_order']
	RA = halos['ra_bcg']
	DEC = halos['dec_bcg']
	Z = halos['z_biwt']
	Nspec = halos['Nspec']
	N200 = halos['N200']
	HVD = halos['HVD']
	RVIR = halos['RVIR']
	SINGLE = halos['single']
	SUB = halos['sub']
	NC4 = halos['nc4']
	M200 = np.zeros(len(N200))	

	#c4_caustics = fits.open(C4.chris_data_root+'/caustics_geom/c4_caustics.fits')[1].data[halos['orig_order']]
	#M200 = c4_caustics['m200']

## Bin on N200
bin = False
if bin == True:
	sort = np.argsort(N200)[::-1]
	HaloID = HaloID[sort]
	RA = RA[sort]
	DEC = DEC[sort]
	Z = Z[sort]
	HVD = HVD[sort]
	N200 = N200[sort]
	RVIR = RVIR[sort]
	SINGLE = SINGLE[sort]
	SUB = SUB[sort]
	NC4 = NC4[sort]
	M200 = M200[sort]
	Nspec = Nspec[sort]

## Cut down to X Nspec
cut_nspec = False
if cut_nspec == True:
	cut_num = 5
	cut = np.where((Nspec > cut_num)&(N200!=0)&(M200!=0)&(SUB<4))[0]
	HaloID = HaloID[cut][:halo_num]
        RA = RA[cut][:halo_num]
        DEC = DEC[cut][:halo_num]
        Z = Z[cut][:halo_num]
        HVD = HVD[cut][:halo_num]
        N200 = N200[cut][:halo_num]
        RVIR = RVIR[cut][:halo_num]
        SINGLE = SINGLE[cut][:halo_num]
        SUB = SUB[cut][:halo_num]
        NC4 = NC4[cut][:halo_num]
        M200 = M200[cut][:halo_num]
        Nspec = Nspec[cut][:halo_num]

	X = np.log10(N200)
	Y = np.log10(M200)


# Initialize Multi-Ensemble Array to hold resultant data
STACK_DATA = []

U.print_separation('# ...Starting Ensemble Loop',type=1)
#  j: index for each ensemble, w.r.t. the FOR loop
#  k: index for each final ensemble, w.r.t. total number of final ensembles
#  l: index for line of sight, w.r.t. total lines of sight of run
for j in range(ens_num):

	# Update S.j
	S.j = j

	# Total Halo Index
	k = run_num*ens_num + j

	# Define Container to Hold Phase Space Data
	PS = Data()

	# Iterate through lines of sight
	U.print_separation('## Working on Ensemble '+str(j),type=1)
	for l in range(line_num):
		print '...Loading galaxy data and projecting line of sight #'+str(l)

		C4.load_project_append(HaloID[stack_range][j*line_num:(j+1)*line_num][l],RA[stack_range][j*line_num:(j+1)*line_num][l],DEC[stack_range][j*line_num:(j+1)*line_num][l],Z[stack_range][j*line_num:(j+1)*line_num][l],M200[stack_range][j*line_num:(j+1)*line_num][l],RVIR[stack_range][j*line_num:(j+1)*line_num][l],HVD[stack_range][j*line_num:(j+1)*line_num][l],PS)

	PS.to_array(['Rdata','Vdata','HaloID','R200','HVD','Z','M200','G_mags','R_Mags','I_Mags'])

	# Check for empty-core clusters, skip those b/c shiftgapper can't handle them
	if np.where((PS.Rdata<0.5)&(np.abs(PS.Vdata)<3000))[0].size == 0:
		STACK_DATA.append( {} )
		print ''
		print '-'*40
		print 'EMPTY CORE CLUSTER'
		print '-'*40
		print ''
		continue

	# Build Ensemble and Run Caustic Technique
	stack_data = S.caustic_stack(PS.Rdata,PS.Vdata,PS.HaloID,np.vstack([PS.M200,PS.R200,PS.HVD,PS.Z]),line_num,feed_mags=True,G_Mags=PS.G_Mags,R_Mags=PS.R_Mags,I_Mags=PS.I_Mags,ens_shiftgap=shiftgapper,edge_int_remove=edge_int_remove)

	# Append other Arrays
	extra = {'_name_':'stack_data'}
	stack_data.update(extra)

	# Append to STACK_DATA
	STACK_DATA.append(stack_data)

# Finished Loop
U.print_separation('#...Finished Ensemble Loop',type=2)


# Make Arrays
make_arrays = False
if make_arrays == True:
	EDGEMASS = np.zeros(ens_num)
	CAUMASS = np.zeros(ens_num)
	EDGEMASS_EST = np.zeros(ens_num)
	CAUMASS_EST = np.zeros(ens_num)
	R200_EST = np.zeros(ens_num)
	ENS_HVD = np.zeros(ens_num)
	EDGESURF = []
	CAUSURF = []
	ENS_R,ENS_V = [], []
	for i in range(ens_num): 
		EDGEMASS[i] = STACK_DATA[i]['ens_edgemass'][0]
		CAUMASS[i] = STACK_DATA[i]['ens_caumass'][0]
		EDGEMASS_EST[i] = STACK_DATA[i]['ens_edgemass_est'][0]
		CAUMASS_EST[i] = STACK_DATA[i]['ens_caumass_est'][0]
		R200_EST[i] = STACK_DATA[i]['ens_r200_est'][0]
		ENS_HVD[i] = STACK_DATA[i]['ens_hvd'][0]
		EDGESURF.append(STACK_DATA[i]['ens_edgesurf'])
		CAUSURF.append(STACK_DATA[i]['ens_causurf'])
		ENS_R.append(STACK_DATA[i]['ens_r'])
		ENS_V.append(STACK_DATA[i]['ens_v'])

	EDGESURF,CAUSURF = np.array(EDGESURF),np.array(CAUSURF)
	ENS_R,ENS_V = np.array(ENS_R),np.array(ENS_V)
	BINN200 = np.array(map(np.median,N200[:halo_num].reshape(ens_num,halo_num/ens_num)))
        BINHVD = np.array(map(np.median,HVD[:halo_num].reshape(ens_num,halo_num/ens_num)))
	BINNSPEC = np.array(map(np.median,Nspec[:halo_num].reshape(ens_num,halo_num/ens_num)))

	d = {'M200':M200,'N200':N200,'RVIR':RVIR,'HVD':HVD,'Nspec':Nspec,'BINN200':BINN200,'BINNSPEC':BINNSPEC,'EDGEMASS':EDGEMASS,'CAUMASS':CAUMASS,'EDGEMASS_EST':EDGEMASS_EST,'CAUMASS_EST':CAUMASS_EST,'R200_EST':R200_EST,'EDGESURF':EDGESURF,'CAUSURF':CAUSURF,'x_range':S.C.x_range,'ENS_R':ENS_R,'ENS_V':ENS_V,'ENS_HVD':ENS_HVD}

	f = open('C4Stack_Ngal'+str(gal_num)+'_Nclus'+str(line_num)+'_Nens'+str(ens_num)+'.pkl','wb')
	output = pkl.Pickler(f)
	output.dump(d)
	f.close()


### Save Data into .pkl Files ###
if write_data == True:
	U.print_separation('##...Starting Data Write',type=2)
	for m in range(ens_num):
		n = run_num*ens_num + m
		U.print_separation("Writing Data for Ensemble #"+str(n),type=2)
		print 'Writing File: '+root+'/C4/'+data_loc+'/'+write_loc+'/Ensemble_'+str(n)+'_Data.pkl'
		pkl_file = open(root+'/C4/'+data_loc+'/'+write_loc+'/Ensemble_'+str(n)+'_Data.pkl','wb')
		output = pkl.Pickler(pkl_file)
		output.dump(STACK_DATA[m])
		output.dump(varib)
		pkl_file.close()

	U.print_separation('#...Finished Data Write',type=2)

## Wait until at least 60 seconds is up
duration = (float(time.asctime()[11:13])*3600+float(time.asctime()[14:16])*60+float(time.asctime()[17:19])) - (float(run_time[11:13])*3600+float(run_time[14:16])*60+float(run_time[17:19]))
if duration < 60:
	time.sleep(60-duration)
	

U.print_separation('## Finished c4_stack.py'+'\n'+'Start:'+'\t'+run_time+'\n'+'End:'+'\t'+time.asctime(),type=1)


