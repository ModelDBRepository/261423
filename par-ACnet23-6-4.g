//genesis

/**********************************************************************
** This simulation script and the files included in this package
** are Copyright (C) 2019 by David Beeman (dbeeman@colorado.edu)
** and are made available under the terms of the
** GNU Lesser General Public License version 2.1
*********************************************************************/

/*======================================================================

  ACnet23 a 'simple yet realistic' model of primary auditory layer 2/3
  based on:

  Beeman D (2013) A modeling study of cortical waves in primary auditory
  cortex. BMC Neuroscience, 14(Suppl 1):P23 doi:10.1186/1471-2202-14-S1-P23
  (http://www.biomedcentral.com/1471-2202/14/S1/P23)

  This simulation script reproduces the results shown in Fig. 1 of:

  David Beeman, Alfred Yu, Joshua Crone (2019) Studying evoked potentials
  in large cortical networks with PGENESIS 2.4.
  BMC Neuroscience 2019, 20(Suppl 1):P46
  https://bmcneurosci.biomedcentral.com/articles/10.1186/s12868-019-0538-0#Sec281

  Further description of the ACnet2 series of models and the installation/use
  of GENESIS and PGENESIS is given at the GENESIS website:
  http://genesis-sim.org/GENESIS/

  Note that this simulation requires the May 2019 official release version
  of PGENESIS 2.4

   This variation of the 2013 model uses a new compiled object
   facdep_rules to provide short term facilitation and depression
   of excitatory synapses of the two populations of neurons:

   Ex_cell_name = "pyr_23"   // name of the excitatory pyramidal cell
   Inh_cell_name = "bask_4" // name of the inhibitory (basket) cell
 
 The main script (script_name) includes:

   synapseinfo.g - optionally included for synapse debugging information
   par-ACnet23_netfuncs.g - defines functions for connecting the network
   par-ACnet23_funcs.g - defines most other functions used by the main script
   protodefs.g - defines the prototype cell elements and cells. Includes:
       pyrchans3.g - channel definitions for pyramidal cells
       FSchans.g - channel definitions for FS (basket) Inh_cell
   par-simple_inputs.g - provides for regular or random pulsed spike trains
   MGB-protodefs.g - prototypes for the input Medial Geniculate Body cells
       uses MGBcell.p and MGBchans.g for channel definitions. These are
       not used in this simulation.
   par-LFPfuncs.g - functions for measuring field potentials
  
 ======================================================================*/

str script_name = "par-ACnet23-6-4.g" // version with default parameters
str RUNID = "6404"        // default ID string used in output file names

// Include the definitions and default values for all global variables
// constants.g defines many flags for simulation options and provides the
// default values of most cell and network model parameters.

include constants3.g

/* for debugging and exploring - see comments in file for details
   Usage:   synapse_info path_to_synchan
   Example: synapse_info /Ex_layer4/Ex_cell[5]/dend/Ex_channel
*/

if (debug_level)
    include synapseinfo.g
end


/* ============ Customizations for this run =========================
   Change these global variables from their default values before
   setting up the simulation.
   ==================================================================
*/

// These need to be defined before the includes

n_slices = 24
n_nodes = n_slices + 2
Ex_cells_per_slice = Ex_NX*Ex_NY / n_slices
Inh_cells_per_slice = Inh_NX*Inh_NY / n_slices

/* Fill the workers list (a string of comma-separated worker numbers) */

output_node = n_slices + 1 // Use a separate node for output to files
// output_node = 0 // share a node with control_node

int worker0 = 1     // Used to identify the first worker node
// int worker0 = 0     // Used to identify the first worker node
// str workers = "0"   // "0" is first in list of nodes
str workers = "1"   // "1" is first in list of nodes
int i
for (i = 2; i <= n_slices; i = i + 1)
// for (i = 1; i < n_slices; i = i + 1)
        workers = workers @ "," @ {i} // add a comma and the next worker
end

if (pdebug)
   echo "Created list of workers = " {workers}
end

/* Optional changes to override definitions in contstants.g */

// RUNID = 6410
// Ex_inh_gmax = 10.0e-9
tmax = 11.1 // original run was for 10.2 sec

// =================================
//   Parallel function definitions
// =================================
/* The function quit_sim is designed to be issued from just one node,
   and causes all nodes, including itself, with "quit@all".  One can quit at
   the command prompt with "quit_sim", or with the QUIT button in the GUI,
   which invokes quit_sim.
*/
function quit_sim
    echo@{control_node} Quitting simulation
    quit@all
end

// ===========-==================
//   Include Function definitions
// ==============================

include  par-ACnet23_funcs-6.g
include  par-netfuncs-5.g

/* For running simulations in batch mode without using the GUI when
   "batch = 1". The function do_batch may need to be edited for the specific
   simulation. Note that it is called last, after all setup of the network
   has been performed. Some global options and parameters used in setup
   must be set at the beginning of this main script.
*/ 
function do_batch  // performed by all nodes
    str path
    if (i_am_output_node)
       change_RUNID {RUNID}
    end
    barrier // all wait to change_RUNID
    if(pdebug)
      echo@{control_node} "Node: " {mynode} "  Memory used: " {getstat -memory}	
    end
    if (i_am_control_node)
        echo "RUNID: " {RUNID}
        if (conn_debug)
            print_avg_syn_number
            print_avg_input_number
        end
        do_reset
        do_reset
       step_tmax
    end
    barrier // wait until the run is finished
    flush_buffers
    if (pdebug)
      echo@{control_node} "Node: " {mynode} "  Memory used: " {getstat -memory}	
    end
    barrier // wait before quitting
    paroff
    quit
end // function do_batch

//===============================
//    Main simulation section
//===============================

/***** set the clocks *****/
setclock  0  {dt}               // set the simulation clock
setclock  1  {out_dt}       // for most output and graphics
setclock  2  {netview_dt}   // for network Vm and Ik display
setclock 3 {facdep_dt}        // for FACDEP weight update
setclock 4 {weights_out_dt} // output and plot the weights
setclock 5 {pulse_interval} // used for oddball_freqs output
if (normalize_weights)  // rarely used option
    setclock 6 {wtscale_dt}
end

//===================================================================
//    Start up parallel nodes - console output goes to control_node
//===================================================================

paron -parallel -silent 1 -nodes {n_nodes}


/* Parse optional arguments with alternative parameter values. For example:
       pgenesis -nox -nodes 26 par-ACnet23-6B.g 6206B 6.0e-6
   The statements below extract a new RUNID and Ex_inh_gmax value from the
   extra arguments. At some point, this will be generalized.
*/
int count = {argc}
if (count > 1)
    RUNID =  {argv 1}
    Ex_inh_gmax = {argv 2}
    echo "Changing parameters - RUNID: "{RUNID} "  gmax = " {Ex_inh_gmax}
else
    echo "No arguments were given. Default values will be used."
end

/* When all nodes but the control_node are waiting at the final barrier,
   only report this every 60 seconds instead of the default 3 seconds.
   This is described in BoG Sec. Sec. 21.9.6.  A longer pvm_hang_time gives
   time to enter interactive commands between the printed dots.
*/
setfield /post pvm_hang_time 60   // give time for connection set up
setfield /post msg_hang_time 3600 // give time for for console input

i_am_control_node = {mynode} == {control_node}
i_am_graphics_node =  (graphics) && ({mynode} == {graphics_node})
i_am_worker_node = ({mynode} >= worker0) && ({mynode} <  {n_slices + worker0})
i_am_output_node =  {mynode} == {output_node}

randseed {seed + {mynode}*347} // this avoids artifacts of each node using
                               // the same random number sequence
if(pdebug)
    echo@{control_node} {nnodes} " nodes in " {nzones} " zones, " {npvmcpu} "CPUs"
    echo@{control_node} I am node {mynode}
    echo@{control_node} Completed startup at {getdate}
end
barrier // wait for everyone to catch up

if ((i_am_control_node) && (pdebug > 0))
  echo "START of network setup: " {getdate}
end
barrier

/* Set up the network on the worker nodes */
if (i_am_worker_node)

/* Including the protodefs file creates prototypes of the channels,
   and other cellular components under the neutral element '/library'.
   Calling the function 'make_prototypes', defined earlier in this
   script, uses these and the cell reader to add the cells.
*/
    include protodefs3-0.g
    // Now /library contains prototype channels, compartments, spikegen

    make_prototypes // This adds the prototype cells to /library

    make_network_slice {default_layer} // Copy cells into network layers

    // make_network should do some of this, but set all synchan gmax values
    set_all_gmax

    /* synchan tau values should not be changed after hsolve SETUP */
    // Change the synchan tau1 and tau2 from the values used in protodefs
    set_all_taus

    // set the random background excitation frequency
    set_frequency {frequency}

    make_hsolve {default_layer}
    if (pdebug)
      echo@{control_node} {mynode} " waiting at barrier 1 after creating slice"
    end
    barrier 1 // wait for every slice to be set up

    // Now connect them
    if (connect_network)
        //if (debug_level)
        if (pdebug > 0)
          echo@{control_node} {mynode} " Starting connection set up: " {getdate}
        end
        // defined in par-netfuncs.g
        connect_cells // connect up the cells in the network layers
        // if (debug_level)
        if (pdebug > 0)
          echo@{control_node} {mynode} " Finished connection set up: " {getdate}
        end
    end
    if (pdebug)
//      echo@{control_node} {mynode} " waiting at barrier 2 after connecting"
      echo "node " {mynode} " waiting at barrier 2 after connecting"
    end
    barrier 2 // wait for every slice to be connected

  // This version uses no random background excitation
  frequency = 0.0
  set_frequency {frequency}

  // set weights and delays
  set_weights {syn_weight}

  set_delays  {prop_delay}

    if (pdebug)
      echo@{control_node} {mynode} " Set wights/delays; waiting at barrier 3"
    end
    barrier 3 // wait for every slice to be connected

  /* set up FACDEP -- setup functions defined in ACnet23_funcs.g */

  if (use_Ex_ex_FACDEP)
    setup_depdep_Ex_ex
    if(debug_level)
        echo "Using facdep_rules2 with step "  {getclock 3}
    end
    useclock /depdep_Ex_ex 3
  end

  if (use_Inh_ex_FACDEP)
    setup_depdep_Inh_ex
    if(debug_level)
        echo "Using facdep_rules2 with step "  {getclock 3}
    end
    useclock /depdep_Inh_ex 3
  end
    if (pdebug)
      echo@{control_node} {mynode} " setup facdep: waiting at barrier 4"
    end
    barrier 4 // wait for facdep setup

else // non-worker nodes
  barrier 1
  barrier 2
  barrier 3
  barrier 4
end // if/else i_am_worker

if ((i_am_control_node) && (pdebug > 0))
  echo "END of network setup (create, connect, facdep): " {getdate}
end
barrier

/* Set up the inputs to the network.  Depending on the type of input
   to be used, include the appropriate file for defining the functions
   make_inputs and connect_inputs.

   par-simple_inputs.g recognizes the input_types "pulsed_spiketrain"
   and "pulsed_randomspike", and the patterns "row" (the default),
   "line" and "box".

   The "ACnet2-6" series of scripts have been specialized to use the two
   types "MGBv" (with two tones and the MGBv input model) and the variant
   "oddball", which alternately applies two tones with probabilities
   input1_prob and 1 - input1_prob. 
*/

/* However, this version uses par-simple_inputs3.g for a pulsed randomspike.
   Unlike earlier versions, this provides thalamic input to Inh_cells also.
   The function connect_inputs2 in par-simple_inputs3.g uses a loop with raddmsg
   to connect the spike generators to cells in a target row and and adjacent
   rows. The connection probability decreaases with distance from the target
   row.

   The older function connect_inputs used rvolumeconnect to connect
   to probabalistically connect to one of the rows, but connecting
   to all the cells in the row with 100% proability if that row is picked.

   To reproduce the published CNS 2019 results, set use_connect_inputs2 = 0
   in constants.g.
*/

if (i_am_control_node)
  include par-simple_inputs3.g
  // Now increase row_sep from 7 to 12
  row_sep = 12

  /* Modify the pulse parameters to reflect a version of Eliades/Boatman expt */

  float pulse_width =  0.015    // width of pulse
  float pulse_delay = 0.1        // delay before start of pulse
  float pulse_interval = 1.0 // time from start of pulse to next (period)
  float spikefreq =  793.7   // frequency "f0" of first input row

  make_inputs {spikefreq} // Create array of network inputs starting at spikefreq

  if(input_paradigm == "oddball") // from ACnet2-9-7_funcs.g
    setup_oddball_inputs 1 2 {input1_prob}
    apply_oddball_inputs // start it off before the first step
    setclock 5 {pulse_interval}
  end
  if(use_connect_inputs2)
    connect_inputs2
  else
    connect_inputs // Connect the inputs to the network
  end
setall_driveweights {drive_weight} // Initialize the weights of input drive
  setfield /MGBv[1]/spikepulse level1 1.0
  if (pdebug)
      echo@{control_node} {mynode} "waiting at barrier 5 after connecting inputs"
  end
  if (conn_debug)
      echo "Messages from /MGBv[1]/spikepulse/spike: "
      rshowmsg  /MGBv[1]/spikepulse/spike
  end
  barrier 5
else
  barrier 5
end // if (i_am_control_node)

if (pdebug)
    echo@{control_node} {mynode} "waiting at barrier 6 after setting up inputs"
end
barrier 6


/* For this version, efield electrodes will be placed in the
   middle of a row, with an offset to avoid being on top of a
   neuron. z = 2 mm above the somata of the network layer
   2/3. Electrodes are placed above row 16. make_electrode_cells applies
   only to workers and creates efield_cells elements in each slice. The fields
   are summed in the output node calculator.
*/
if (calc_LFP)
    include par-LFPfuncs.g
    use_syn_currents = 1
    make_electrode_cells "Im_field" {cellpath}[] \
      {(Ex_NX-1)*Ex_SEP_X/2.0} {39.5*Ex_SEP_Y} 0.002 // use efield_cells
end // if (calc_LFP)
barrier 7
make_electrode_cells_messages "Im_field" // only applies to workers
barrier 8 // wait for everyone to catch up

if (pdebug)
    echo@{control_node} "waiting at barrier 8 after setting up efield_cells"
end

if (i_am_output_node)
    if(input_paradigm == "oddball")
        setup_oddball_output
    end
    // Create disk_out elements /output/{net_efile}, {net_ifile}, {net_EPSC_file}
    make_network_out

    // set up the calculator and file for summed EPSCs
    if (calc_EPSCsum)
        do_EPSCsum_out
    end
    if (calc_LFP)
    /* include these for fields arising from only synchans or thalamic drive
      make_electrode "Electrode_synchans" {(Ex_NX-1)*Ex_SEP_X/2.0} \
        {39.5*Ex_SEP_Y} 0.002 // uses efield2
      make_electrode_synpath "Ex_drive_synchans" \
          {(Ex_NX-1)*Ex_SEP_X/2.0} {39.5*Ex_SEP_Y} 0.002
      make_electrode_synpath "Ex_inh_synchans" \
          {(Ex_NX-1)*Ex_SEP_X/2.0} {39.5*Ex_SEP_Y} 0.002
      do_electrodes_out
   */
      do_electrodes_cells_out
    end
    if (pyrVm_output)
        make_pyrVm_out 792 // middle of row 16
    end
end // if (i_am_output_node)

if (pdebug)
  echo@{control_node} "node " {mynode} " waiting at barrier 10 to setup outputs"
end
barrier 10 // everyone waits

if (i_am_worker_node)
    make_network_out_messages // function checks flags
    if (calc_EPSCsum)
         make_EPSCsummer_messages
    end
/* Don't calculate or output these
    if (calc_LFP)
        make_electrode_messages  "Electrode_synchans" {cellpath}
        make_electrode_synpath_messages  "Ex_drive_synchans" \
          {cellpath}  {Ex_drive_synpath}
        make_electrode_synpath_messages  "Ex_inh_synchans" \
          {cellpath}  {Ex_inh_synpath}
    end
*/
    if (pyrVm_output)
        make_pyrVm_out_messages 792 // these return if cell_numm not on mynode
    end
end // if (i_am_worker_node)

if (pdebug)
  echo@{control_node} "node " {mynode} " waiting at barrier 11 to setup output messages"
end
barrier 11 // everyone waits


/* if running in batch mode, the function do_batch may need to be modified */
if (batch)
  do_batch
end

/* Otherwise, run the simulation and stop at a genesis prompt on control_node */
if(i_am_control_node)
    echo "Network of "{Ex_NX}" by "{Ex_NY}" excitatory cells with separations" \
        {Ex_SEP_X}" by "{Ex_SEP_Y}
    echo "and "{Inh_NX}" by "{Inh_NY}" inhibitory cells with separations" \
         {Inh_SEP_X}" by "{Inh_SEP_Y}
    echo "Random number generator seed initialized to: " {seed}

    print_avg_syn_number
    barrier 12
    echo "resetting all nodes"
    reset@all
    echo "stepping all nodes"
    step_tmax
    barrier 13
else
  barrier 12
  barrier 13
end // if(i_am_control_node)

/* for all nodes: these will return if the cell is not on that node */
if (debug_level>1)
    print_Ex_synapse_info 792
end

// All the other nodes will get stuck at this barrier
// and the genesis prompt will be available in Node 0
if (!i_am_control_node)
   barrier 15
end
