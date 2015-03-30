## c4_stack.py
"""
This program uses the Caustic Technique to estimate the masses of galaxy clusters,
after having applied a stacking technique to create ensemble clusters.
This script does strictly self stacking or bin stacking.

"""
## Import Modules
print '...importing modules'
from caustic_stack import *
from stack_class import *

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
keys = ['c','h','H0','q','beta','fbeta','r_limit','v_limit','photo_data','spec_data','halo_num','gal_num','line_num','method_num','write_loc','data_loc','root','self_stack','scale_data','write_data','run_time','init_clean','small_set','run_los','bootstrap','run_num','ens_num','cell_num','stack_range','mass_mix','mass_scat','bootstrap_num','bootstrap_rep','avg_meth','cent_offset','center_scat','new_halo_cent','true_mems','mirror','edge_perc','lightcone','mm_est','data_set']

varib = ez.create(keys,locals())
varib.update({'_name_':'varib'})

## INITIALIZATION ##
S = Stack(varib)
U = Universal(varib)
C4 = CFOUR(varib)

U.print_separation('## Running c4_stack.py')
names = ['run_time','','run_num','gal_num','line_num','cell_num','ens_num','halo_num','method_num','avg_meth','','write_data','scale_data','run_los','mirror','new_halo_cent','init_clean','bootstrap','','bootstrap_num','bootstrap_rep','','data_set','data_loc','write_loc','root']
U.print_varibs(names,varib)

## Load Halo Data
U.print_separation('# ...Loading Halos',type=1)
HaloID,HaloData = C4.load_halos(sort=True,avg_centers=True)
RA,DEC,Z,RVIR,N200,SINGLE,SUB,NC4 = HaloData

# Load in Halo Data from previously written file, or write file if doesn't exist
try:
	f = open(root+'/C4/'+data_loc+'/halo_arrays.pkl','rb')
	input = pkl.Unpickler(f)	
	run_dict = input.load()
	globals().update(run_dict)
	f.close()
except IOError:
	try:
		f = open(root+'/C4/'+data_loc+'/halo_arrays.pkl','wb')	
		output = pkl.Pickler(f)
		run_dict = {'HaloID':HaloID,'RA':RA,'DEC':DEC,'Z':Z,'N200':N200,'RVIR':RVIR,'SINGLE':SINGLE,'SUB':SUB,'NC4':NC4}
		output.dump(run_dict)
		f.close()
	except IOError:
		pass

# Unpack HaloData array into local namespace

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

		C4.load_project_append(HaloID[stack_range][j*line_num:(j+1)*line_num][l],RA[stack_range][j*line_num:(j+1)*line_num][l],DEC[stack_range][j*line_num:(j+1)*line_num][l],Z[stack_range][j*line_num:(j+1)*line_num][l],RVIR[stack_range][j*line_num:(j+1)*line_num][l],PS)

	PS.to_array(['Rdata','Vdata','HaloID','R200','G_mags','R_Mags','I_Mags'])

	# Build Ensemble and Run Caustic Technique
	stack_data = S.caustic_stack(PS.Rdata,PS.Vdata,PS.HaloID,np.vstack([PS.M200,PS.R200,PS.HVD]),line_num,feed_mags=True,G_Mags=PS.G_Mags,R_Mags=PS.R_Mags,I_Mags=PS.I_Mags)

	# Append other Arrays
	extra = {'pro_pos':PS.pro_pos,'_name_':'stack_data'}
	stack_data.update(extra)

	# Append to STACK_DATA
	STACK_DATA.append(stack_data)

# Finished Loop
U.print_separation('#...Finished Ensemble Loop',type=2)

### Save Data into .pkl Files ###
if write_data == True:
	U.print_separation('##...Starting Data Write',type=2)
	for m in range(ens_num):
		n = run_num*ens_num + m
		U.print_separation("Writing Data for Ensemble #"+str(n),type=2)
		print 'Writing File: '+root+'/OSDCStacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(n)+'_Data.pkl'
		pkl_file = open(root+'/OSDCStacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(n)+'_Data.pkl','wb')
		output = pkl.Pickler(pkl_file)
		output.dump(STACK_DATA[m])
		output.dump(varib)
		pkl_file.close()

	U.print_separation('#...Finished Data Write',type=2)

## Wait until at least 60 seconds is up
duration = (float(time.asctime()[11:13])*3600+float(time.asctime()[14:16])*60+float(time.asctime()[17:19])) - (float(run_time[11:13])*3600+float(run_time[14:16])*60+float(run_time[17:19]))
if duration < 60:
	time.sleep(60-duration)
	

U.print_separation('## Finished millennium_stack.py'+'\n'+'Start:'+'\t'+run_time+'\n'+'End:'+'\t'+time.asctime(),type=1)


