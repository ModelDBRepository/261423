// genesis - constants.g

/**************************************************************************

These are the definitions and default values of the global variables that
are used in the par-ACnet23-6 series of scripts.

Note the many flags to enable or disable various options.

Strings are defined to specify the paths to synapses on cells. e.g.
'Ex_inh_synpath' is the path to the compartment and the inhibitory synapse
on the excitatory (pyramidal) cell.

**************************************************************************/

float tmax = 10.2         // max simulation run time (sec)

/***** time steps for clocks *****/
float dt = 20e-6      // clock 0 - simulation time step
float out_dt = 0.0001     // clock 1 - for most output and graphics
float netview_dt = 0.0002 // clock 2 - for network Vm and Ik display
float facdep_dt = 0.001   // clock 3 - for fac/dep weight update, or spike counting
float weights_out_dt = 0.05  // clock 4 - output and plot the weights
float pulse_interval = 1.0   // clock 5 - temporary value used for oddball out
float wtscale_dt = 0.010     // clock 6 - for weight normalization and freq monitor

// Booleans indicating the type of calculations or output
int debug_level = 0  // display additional information during setup
                     // Use a higher number for more verbosity
int pdebug = 0       // print parallel debugging messages
int conn_debug = 0
int batch = 1        // if (batch) run the default simulation without graphics
int graphics = 0     // display control panel, graphs, optionally net view
int netview_output = 0  // Record network output (soma Vm) to a file
int binary_file = 0 // if 0, use asc_file to produce ascii output
            // else use disk_out to produce binary FMT1 file
int write_asc_header = 1 // write header information to ascii file
int EPSC_output = 0  // output Ex_ex EPS currents to file
int calc_EPSCsum = 0 // calculate summed excitatory post-synaptic currents
int pyrVm_output = 0 // For debugging, output Vm of pyr_23[792]
int calc_LFP = 1     // calculate summed local field potentials from efield
int use_syn_currents = 1 // calculate the efield using only synchan currents
// use leak and all channel currents, without capacitive currents

// these  options aren't available in this version
int use_all_currents = 0 
int calc_firing_freqs = 0 // calculate binned average firing rates
int weights_out = 0 // record and plot weight averages
int enhanced_output = 0 // if 1, output synaptic weights and firing rates
int use_weight_decay = 0 // Use exponential decay of weights with distance

int use_prob_decay = 1 // Use connection probablility exp(-r*r/(sigma*sigma))
int use_probdecay_function = 1 // Use function planarconnect_probdecay
int use_connect_inputs2 = 0 // if 1, use raddmsg instead of rvolumeconnect
int connect_network = 1  // Set to 0 for testing with unconnected cells

/***** Flags for synaptic plasiticity *****/

int normalize_weights = 0 // perform synaptic scaling to average weight
int change_weights = 1 // Allow weight changes. if = 0, facdep_update will be
                       // called, but return.  If < 0, facdep_applier not created.
int use_Ex_ex_FACDEP = 1 // set up FACDEP for excitatory synapses on Ex_cells
int use_Inh_ex_FACDEP = 1 // set up FACDEP for excitatory synapses on Inh_cells

// This has not yet been implemented:
// int use_Ex_inh_FACDEP = 0 // set up FACDEP for inhibitory synapses on Ex_cells

/* Definitions for the parallel version */

// typically n_nodes = n_slices + 2, so with 24 slices, use
// pgenesis -nox -nodes 26 par-ACnet23-6-2.g
// unless output = 0, then use one less node.

int n_slices            // each slice will do a horizontal slice of network
int n_nodes             // total number of nodes
str workers             // will be a comma-separated lis of workers
int worker0             // The first worker node
int control_node        // mode for console i/o and injection circuitry
int graphics_node	// node for XODUS display
int output_node         // which node is handling output to file
int i_am_control_node, i_am_worker_node // booleans indicating the function
int i_am_output_node, i_am_spare_node   //   assigned to this node
int i_am_graphics_node

control_node = 0 // This should be the terminal window
graphics_node = 0  // If there is a GUI, put it in this node also

/* Dimensions and grid spacing of the excitatory and inhibitory networks.

   The default size is a square patch of cortex with 48 x 48
   excitatory cells and 24 x 24 inhibitory, spanning two octaves.  For
   a longer piece, spanning six octaves (about 1.92 x 5.76 mm), use
   Ex_NY = 144 and Inh_NY = 72.  

   For now, assume the same grid for all layers. This may change.
*/

int Ex_NX = 48; int Ex_NY = 48
int Inh_NX = 24; int Inh_NY = 24

/* Neurons will be placed on a two dimensional NX by NY grid, with points
   SEP_X and SEP_Y apart in the x and y directions.

   Cortical networks typically have pyramidal cell separations on the order
   of 10 micrometers, and can have local pyramidal cell axonal projections
   of up to a millimeter or more.  For small network models, one sometimes
   uses a larger separation, so that the model represents a larger cortical
   area.  In this case, the neurons in the model are a sparse sample of the
   those in the actual network, and each one of them represents many in the
   biological network.  To compensate for this, the conductance of each
   synapse may be scaled by a synaptic weight factor, to represent the
   increased number of neurons that would be providing input in the actual
   network.  Here, we use a separation of 40 um that is larger than the
   spacing between cells in A1. With no axonal delays, actual spacing is
   irrelevant.

   The names of the grids for layer23 cells used here and in the
   included scripts are Ex_layer23 and  Inh_layer23. These can be
   changed to use "layer4" or any other layer name by changing the
   global variable "default_layer".
*/

// 40 micrometer spacing between cells
float Ex_SEP_X = 40e-6
float Ex_SEP_Y = 40e-6 
float Inh_SEP_X = 2*Ex_SEP_X  // There are 1/4 as many inihibitory neurons
float Inh_SEP_Y = 2*Ex_SEP_Y

/* "SEP_Z" should be set to the actual layer thickness, in order to allow
   possible random displacements of cells from the 2-D lattice.  Here, it
   needs to be large enough that any connections to distal dendrites will
   be within the range -SEP_Z to SEP_Z.
*/

float Ex_SEP_Z = 1.0; float Inh_SEP_Z = Ex_SEP_Z

/* For the parallel version, the network will be divided into horizontal
   slices of dimension NX x NY/n_slices.

   Be sure that NY is a multiple of n_slices !!!
*/

int Ex_cells_per_slice //  = Ex_NX*Ex_NY / n_slices
int Inh_cells_per_slice // = Inh_NX*Inh_NY / n_slices


/* Specification of stimulation input patterns (distribution of inputs)
   and type of input (e. g. a steady spike train or a randomspike source,
   gated with a pulsegen, or a more realistic thalamic input model "MGBv")
*/

// str input_type = "pulsed_spiketrain"  // pulsed spike train input

/*  pulsed  Poisson-distributed input with a specifed average frequency.
    This is implemented with GENESIS randomspike elements driving spikegens,
    gated by a pulsegen.
*/

str input_pattern = "row"     // input goes to an entire row of cells
str input_type = "pulsed_randomspike"

str input_paradigm = "" // default is one frequency input
str input_paradigm = "oddball" // two inputs

/* Settings for weighted random selections between two MGBv inputs of
   different frequency.
*/

/* the  "oddball" input uses the MGBv model or simple-inputs, but switches
   between two input rows with "frequent" input 1 having a larger probability
   than the "oddball" or "infrequent" input 2.
*/
// str input1_prob = 0.82 // input 1 has 0.82 probablitiy and input 2 has 0.18
str input1_prob = 1.0 // For this version, use 100% probability for input 1

/* The flag 'use_stim_order' determines the type of sequence to use with
   input_paradigm = "oddball".

   if use_stim_order = 0, a random sequence of frequent and infrequent tones
   will be used, with the probability for input tone 1 to be given by
   'input1_prob', typically 0.82. For a repeated single tone only,
   use input1_prob = 1.0.

   use_stim_order =  1 is a variation used with input_paradigm = "oddball"
   that uses non-random sequences for the two tones, taken from a file
   {stim_order_file}.

   use_stim_order = 2 produces an initial frequent tone, and then
   trains of four frequent tones, followed by one infrequent tone.

   use_stim_order = 3 produces trains of four frequent tones, followed by
   two infrequent tones. These are useful for short test runs.
*/
   
int use_stim_order = 0 // 1: stimulus order (1st, 2nd, 3rd, .. oddball from file)
     // 2: programmed sequence 1 1 1 1 1 2 1 1 1 1 2 1 1 1 1 2 ....

int hflag = 1    // use hsolve if hflag = 1
int hsolve_chanmode = 4  // chanmode to use if hflag != 0
int use_sprng = 1 // Use SPRNG random number generator, rather than default RNG

// default values of input pulse parameters that may be modified below
float pulse_width =  0.015    // width of pulse
float pulse_delay = 0.1        // delay before start of pulse
float pulse_interval = 1.0 // time from start of pulse to next (period)

/* Customize these strings and parameters to modify this simulation for
   other excitatory or inhibitory cells.
*/

// This global parameter is used to define layer names, such as
// "/Ex_"@{default_layer} == "Ex_layer4"
str default_layer = "layer23" // could be replaced by "layer4"

str Ex_cell_name = "pyr_23"   // name of the excitatory cell
str Inh_cell_name = "bask_4" // name of the inhibitory (basket) cell

str Ex_cellfile = "pyr_23_asym.p"  // name of the excitatory cell parameter file
str Inh_cellfile = "bask.p"  // name of the inhibitory cell parameter file

/***************************************************************************
 Paths to synapses on cells: cell_synapse = compartment-name/synchan-name
 e.g. Ex_inh_synpath is path to inhibitory synapse on excitatory cell

****************************************************************************/

str Ex_ex_synpath = "oblique2b/AMPA_pyr" // lower oblique
str Ex_inh_synpath = "apical0/GABA_pyr" // pyr apical trunk GABA
str Ex_drive_synpath = "basal2b/AMPA_pyr"
str Inh_ex_synpath = "dend/AMPA_bask"  // bask dend AMPA
str Inh_inh_synpath = "soma/GABA_bask" // bask soma GABA
str Ex_bg_synpath = "apical1/AMPA_pyr"  // not used when freq = 0

// Excitatory drive inputs - path to synapse on Ex and Inh cells to apply drive
str Inh_drive_synpath = "dend/AMPA_bask_drive"  // drive -> bask dendrite

// wildcard list of cells to be used for LFP calculations will be
// cellpath @ "[]"
str cellpath = "/Ex_" @ {default_layer} @ "/" @ {Ex_cell_name}
str solvepath = "solver" // relative path from cell to the hsolve element

/* Synaptic parameters */

float Ex_ex_gmax = 15e-9   // Ex_cell ex synapse
float Ex_inh_gmax = 4.0e-9  // Ex_cell inh synapse
float Inh_ex_gmax = 0.6e-9  // Inh_cell ex synapse

float Inh_inh_gmax = 0.0e-9 // Inh_cell inh synapse 

float Ex_bg_gmax = 80e-9  // Ex_cell background excitation

// Thalamic drive is to Ex and Inh cells only
float Ex_drive_gmax = 30e-9 // Ex_cell thalamic input
float Inh_drive_gmax = 0.8e-9 // Inh_cell thalamic input

// time constants for dual exponential synaptic conductance

//These will be used for all excitatory (AMPA) channels
float tau1_ex = 0.001  // rise time
float tau2_ex =  0.003 // decay time

// GABA inhibition from Inh (basket) cells
float tau1_inh = 0.005    // rise time for inhibitory synapses
float tau2_inh =  0.008   // decay time for inhibitory synapses

// make a special case for Inh cell excitatory channels
float tau1_Inh_ex = 0.003
float tau2_Inh_ex = 0.003

// Poisson distributed random excitation frequency of Ex_cells
// NOTE: For use with hsolve, the synchan frequency must be set to a non-zero
// value before the solver setup.  Then it may be set to any value, including zero.

float frequency = 5.0


/****************************************************************************
  flags and parameters for the Varela phenomenological model for synpaptic
  facilitation and/or depression. This uses the new (June 2016) GENESIS 2.4
  object 'facdep_rules'.
*****************************************************************************/

/***** flags and globals for FACDEP *****/

float avg_weight = 1.0 // Desired average synaptic weight per target synchan

int rand_delay = 0  // Assign delays from 0 - 2*delay
// assignment of random weights is only implemented for use_Ex_ex_FACDEP
int rand_weight = 0  // Assign weights  min_weight*weight - max_weight*weight

// Parameters used for FACDEP update

/* Facilitation/Depression parameters for RS, FS, and LTS  cells, taken from
best fits for excitatory synapses in layer 4 of mouse auditory cortex:
*/
float Ex_dD1 = 0.8; float Ex_tau_D1 = 0.2   // fast, strong depression
float Ex_dD2 = 0.9; float Ex_tau_D2 = 3.0  // slow, weak depression

float Inh_dD1 = 0.55; float Inh_tau_D1 = 0.2   // fast, strong depression
float Inh_dD2 = 0.9; float Inh_tau_D2 = 3.0 // slow, weak depression
         
// Set random number seed. If seed is 0, set randseed from clock,
// else from given seed
int seed = 0  // Simulation will give different random numbers each time

/****** the seed is set here to reproduce the results in Beeman (2013) *****/
int seed = 1369497795

if (use_sprng)
    setrand -sprng
end

if (seed)
    randseed {seed}
else
    seed = {randseed}
end

// Label to appear on the graph
str graphlabel = "Vm of row center cell"
str net_efile = "Ex_netview"  // filename prefix for Ex_netview data
str net_ifile = "Inh_netview" // filename prefix for Inh_netview data
str net_EPSC_file = "EPSC_netview" // filename prefix for Ex_ex_synpath Ik (EPSCs)
str EPSC_sum_file = "EPSC_sum" // filename prefix for summed Ex_ex_synpath Ik
str sum_file = "run_summary"    // text file prefix for summary of run params
str pyrVm_file = "pyrVm"
str avg_weight_file = "average_weight"
str stim_order_file = "Stimulus_order_binary.txt" // .txt

// These definitions depend on MGBv_input.g or simple_inputs.g
// int Ninputs = Ex_NY - 16 // Number of auditory input channels from the thalamus (MGB)
// In this case, there will be two inputs MGBv[1] and MGBv[2]
int Ninputs = 2

float drive_weight = 1.0 // Default weight of all input drive connections
float octave_distance = 0.96e-3 // approx 1 mm/octave - integer rows/octave = 24
// octave_distance = Ex_SEP_Y // just one row

int rows_per_octave = {round {octave_distance/Ex_SEP_Y}}

/* input_spread is the number of rows below and above the "target row"
   getting thalamic input.  Some typical values can be found in Miller, et
   al.  (2001) Neuron 32:151-160. The implementation in simple_inputs.g
   function 'connect_inputs' creates connections with an exponentially
   decaying probablility to adjacent rows +/- input_spread.
*/

 int input_spread = {round {rows_per_octave/6.0}} // 1/6 octave
// int input_spread = 0 // no spread of thalamic connections to other rows

/* input_delay and input_jitter provide a decorrelation between the inputs
   to a row by introducing a random delay to each cell uniformly
   distributed between input_delay*(1 - input_jitter) and input_delay*(1 +
   input_jitter).  A reasonable amount of decorrelation can be provided
   with input_delay = 0.002 seconds, and input_jitter = 0.4.
*/

float input_delay, input_jitter // used in simple-inputs.g or MGBv_input.g
// float input_delay = 0.002 // 2 msec is adequate for the original pyr_4 cell
float input_delay = 0.010 // provides sufficient de-correlation for pyr_23
float input_jitter = 0.4

// float spike_jitter = 0.0005 // 0.5 msec jitter in thalamic inputs
/* Richardson, et al. (2009) J. Neurosci. 2009;29 6406-6417
   Mean EPSP jitter for thalamic inputs in the auditory thalamocortical
   system was 0.41 +/- 0.13 msec.
*/

float spike_jitter = 0.0  //default is no jitter in arrival time

/* parameters for synaptic connections */

float syn_weight = 1.0 // synaptic weight, effectively multiplies gmax

/* 
   prop_delay is the delay per meter, or 1/cond_vel.  The value often used
   corresponds to Shlosberg et al. 2008 rat somatosenory cortex axonal cond
   velocities of RS and Martinotti cell axons ranging from 0.2 to 0.3
   m/sec.  With a value of 0.25 m/sec, cells 1 mm apart have a 4 msec
   conduction delay

   However this value is for longer distance interlaminar axons between
   layer 5 and layer 1 in rat somatosensory cortex.  Estimates for shorter
   distance (< a few mm) intralaminar unmeyelinated axons are in the range
   of 60-90 mm/s. (Salin and Price, 1996).

   The velocity used previously with the pyr_4 PC model was 0.08 m/sec,
   giving a delay of 12.5 msec per mm of separation between somatic spike
   generator and synapse. The current layer 23 model uses a velocity of
   0.125 m/sec with a propagation delay of 8.0  sec/m or msec/mm.
   
   This can have a signigicant effect in the delay of the onset
   of inhibition or excitation, as connected cells 400 um apart
   can have a delay of about 3 msec.
*/
float prop_delay = 8.0 //  delay per meter, or 1/cond_vel

