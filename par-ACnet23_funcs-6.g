//genesis

/**********************************************************************
** This simulation script and the files included in this package
** are Copyright (C) 2013 by David Beeman (dbeeman@colorado.edu)
** and are made available under the terms of the
** GNU Lesser General Public License version 2.1
** See the file copying.txt for the full notice.
*********************************************************************/

/*======================================================================

  Summary of function definitions for the "ACnet2-9" series of simulations

  * functions to set synchan gmax and tau1, tau2 parameters. e.g.
    set_Ex_ex_gmax, set_Inh_ex_gmax, set_Ex_bg_gmax, etc.
    set_all_gmax -- set all gmax to Ex_ex_gmax, ..., Ex_bg_gmax
    set_all_taus -- set all synchan tau1,tau2 to values 'tau1_ex', etc.

  * functions for setting FACDEP parameters for the Ex and Inh
    facdep_rules elements, e.g.
    set_facdep_dt, set_Ex_dD1, set_Ex_D1_tau, ..., set_Inh_D1_tau

  * functions to output results
    make_output, make_header, do_disk_out, do_network_out, ...

  * functions for applying 'oddball' stimuli
    make_trial_params, apply_oddball_inputs, setup_oddball_output,
    setup_oddball_inputs(input1, input2, prob1)

  * functions to calculate and output the summed Ex_ex synaptic currents
    make_EPSCsummer, do_EPSCsum_out

    Note that functions for calculating local field potentials are in
    LFPfuncs.g. Network creation functions are in ACnet2-9_netfuncs.g.

  * functions for post-run output
    do_run_summary

  * functions that are usually invoked by the GUI
    change_RUNID, do_reset, step_tmax


  * functions to set up FACDEP 
    setup_depdep_Ex_ex, setup_depdep_Inh_ex

 ======================================================================*/

// =============================
//   Function definitions
// =============================

/**** set synchan parameters ****/

function set_Ex_ex_gmax(value)  // excitatory synchan gmax in Ex_cell
   float value                  // value in nA
    str Ex_layer = "Ex_" @ {default_layer}
   Ex_ex_gmax = {value}*1e-9	// use this for driver input also
   setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_ex_synpath} gmax {Ex_ex_gmax}
end

function set_Ex_drive_gmax(value)  // thalamic drive synchan gmax in Ex_cell
   float value                  // value in nA
    str Ex_layer = "Ex_" @ {default_layer}
   Ex_drive_gmax = {value}*1e-9
   setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_drive_synpath} gmax {Ex_drive_gmax}
end

function set_Ex_bg_gmax(value)  // background ex  synchan gmax in Ex_cell
   float value                  // value in nA
    str Ex_layer = "Ex_" @ {default_layer}
   Ex_bg_gmax = {value}*1e-9
   setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_bg_synpath} gmax {Ex_bg_gmax}
end

function set_Inh_ex_gmax(value)   // excitatory synchan gmax in Inh_cell
   float value			  // value in nA
   str Inh_layer = "Inh_" @ {default_layer}
   Inh_ex_gmax = {value}*1e-9	  // use this for driver input also
   setfield /{Inh_layer}/{Inh_cell_name}[]/{Inh_ex_synpath} gmax {Inh_ex_gmax}
end

function set_Inh_drive_gmax(value) // thalamic drive synchan gmax in Inh_cell
   float value			  // value in nA
   str Inh_layer = "Inh_" @ {default_layer}
   Inh_drive_gmax = {value}*1e-9
   setfield /{Inh_layer}/{Inh_cell_name}[]/{Inh_drive_synpath} gmax {Inh_drive_gmax}
end

function set_Ex_inh_gmax(value)  // inhibitory synchan gmax in Ex_cell
   float value			  // value in nA
   str Ex_layer = "Ex_" @ {default_layer}
   Ex_inh_gmax = {value}*1e-9
   setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_inh_synpath} gmax {Ex_inh_gmax}
end

function set_Inh_inh_gmax(value)  // inhibitory synchan gmax in Inh_cell
   float value			  // value in nA
   str Inh_layer = "Inh_" @ {default_layer}
   Inh_inh_gmax = {value}*1e-9
   setfield /{Inh_layer}/{Inh_cell_name}[]/{Inh_inh_synpath} gmax {Inh_inh_gmax}
end

function set_all_gmax // set all gmax to Ex_ex_gmax, ..., Ex_bg_gmax
    str Ex_layer = "Ex_" @ {default_layer}
    str Inh_layer = "Inh_" @ {default_layer}

    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_ex_synpath} gmax \
        {Ex_ex_gmax}
    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_inh_synpath} gmax \
        {Ex_inh_gmax}
    setfield /{Inh_layer}/{Inh_cell_name}[]/{Inh_ex_synpath} gmax \
        {Inh_ex_gmax}
    setfield /{Inh_layer}/{Inh_cell_name}[]/{Inh_inh_synpath}  gmax \
        {Inh_inh_gmax}
    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_drive_synpath} gmax \
        {Ex_drive_gmax}
    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_bg_synpath} gmax \
        {Ex_bg_gmax}
    setfield /{Inh_layer}/{Inh_cell_name}[]/{Inh_drive_synpath}  gmax \
        {Inh_drive_gmax}
end

/* NOTE: Functions to set synchan tau1 and tau2 values should be called
   only prior to hsolve SETUP.  Unlike the case with gmax, which may be
   changed if followed by a reset, changing the taus of hsolved synchans
   causes erroneous results.
*/
function set_all_taus  // assume all AMPA synchans have same taus
    str Ex_layer = "Ex_" @ {default_layer}
    str Inh_layer = "Inh_" @ {default_layer}

    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_ex_synpath} tau1 {tau1_ex} \
        tau2 {tau2_ex}
    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_drive_synpath} tau1 {tau1_ex} \
        tau2 {tau2_ex}
    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_bg_synpath} tau1 {tau1_ex} \
        tau2 {tau2_ex}
    setfield /{Inh_layer}/{Inh_cell_name}[]/{Inh_ex_synpath} tau1 {tau1_Inh_ex} \
        tau2 {tau2_Inh_ex}
    // GABA inhibition from Inh (basket) cells 
    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_inh_synpath} tau1 {tau1_inh} \
        tau2 {tau2_inh}
    setfield /{Inh_layer}/{Inh_cell_name}[]/{Inh_inh_synpath} tau1 {tau1_inh} \
        tau2 {tau2_inh}

end

/***** functions for setting FACDEP parameters *****/

extern setup_depdep_Ex_ex // defined later belowx
extern setup_depdep_Inh_ex // defined later below

/* Be sure that the strings for the default facdep applier element
   and synchan path has been set before using these, which are
   invoked by the GUI.
*/

function set_facdep_dt(value)
    float value
    facdep_dt = value
    setclock 3 {facdep_dt}
end


/* set facdep  (dep+dep) parameters for the Ex (pyramidal) cells */

function set_Ex_dD1(value)
    float value
    Ex_dD1 = value
    if ({exists /depdep_Ex_ex})
        setup_depdep_Ex_ex
    end
end

function set_Ex_tau_D1(value)
    float value
    Ex_tau_D1 = value
    if ({exists /depdep_Ex_ex})
        setup_depdep_Ex_ex
    end
end

function set_Ex_dD2(value)
    float value
    Ex_dD2 = value
    if ({exists /depdep_Ex_ex})
        setup_depdep_Ex_ex
    end
end

function set_Ex_tau_D2(value)
    float value
    Ex_tau_D2 = value
    if ({exists /depdep_Ex_ex})
        setup_depdep_Ex_ex
    end
end


/* set facdep parameters for the Inh (basket) cells */
function set_Inh_dD1(value)
    float value
    Inh_dD1 = value
    if ({exists /depdep_Inh_ex})
        setup_depdep_Inh_ex
    end
end

function set_Inh_tau_D1(value)
    float value
    Inh_tau_D1 = value
    if ({exists /depdep_Inh_ex})
        setup_depdep_Inh_ex
    end
end

function set_Inh_dD2(value)
    float value
    Inh_dD2 = value
    if ({exists /depdep_Inh_ex})
        setup_depdep_Inh_ex
    end
end

function set_Inh_tau_D2(value)
    float value
    Inh_tau_D2 = value
    if ({exists /depdep_Inh_ex})
        setup_depdep_Inh_ex
    end
end

function set_frequency(value) // set Ex_ex average random firing freq
    float value
    str Ex_layer = "Ex_" @ {default_layer}
    frequency = value
    setfield /{Ex_layer}/{Ex_cell_name}[]/{Ex_bg_synpath} frequency {frequency}
end

/**** functions to output results ****/

function make_output(rootname) // asc_file to {rootname}_{RUNID}.txt
    str rootname, filename
    if ({exists {rootname}})
        call {rootname} RESET // this closes and reopens the file
        delete {rootname}
    end
    filename = {rootname} @ "_" @ {RUNID} @ ".txt"
    create asc_file {rootname}
    setfield ^    flush 0  leave_open 1 filename {filename} float_format %0.8g
    setclock 1 {out_dt}
    useclock {rootname} 1
end

/* This function returns a string to be used for a one-line header at the
   beginning of an ascii file that contains values of Vm or Ik for each
   cell in the network, at time steps netview_dt. When binary file output
   is used with the disk_out object, the file contains information on the
   network dimensions that are needed for the xview object display.  When
   asc_file is used to generate the network data, this information is
   provided in a header of the form:

   #optional_RUNID_string Ntimes start_time dt NX NY SEP_X SEP_Y x0 y0 z0

cell in the network, at time steps netview_dt.

   This header may be read by a data analysis script, such as netview.py.

   The line must start with "#" and can optionally be followed immediately
   by any string.  Typically this is some identification string generated
   by the simulation run.  The following parameters, separated by blanks or
   any whitespace, are:

   * Ntimes - the number of lines in the file, exclusive of the header

   * start_time - the simulation time for the first data line (default 0)

   * dt - the time step used for output (netview_dt)

   * NX, NY - the integer dimensions of the network

   *  SEP_X, SEP_Y - the x,y distances between cells (optional)

   * x0, y0, z0 - the location of the compartment (data source) relative to
     the cell origin (often ignored)

   A typical header generated by this simulation is:

   #B0003	5000 0.0 0.0002	 48 48 4e-05  4e-05 0.0	 0.0	0.0
*/

function make_header(diskpath)
    str diskpath
    int Ntimes = {round {tmax/netview_dt}} // outputs t = 0 thru tmax - netview_dt
    str header_str = "#" @ {RUNID} @ "  " @ {Ntimes} @ " 0.0 " @ {netview_dt} @ "  "
    if(diskpath == {net_efile})
        header_str = {header_str} @ {Ex_NX} @ " " @ {Ex_NY} @ " " @ {Ex_SEP_X} @ \
        "  " @ {Ex_SEP_Y} @ " 0.0  0.0  0.0"
    elif(diskpath == {net_ifile})
        header_str = {header_str} @ {Inh_NX} @ " " @ {Inh_NY} @ " " @ {Inh_SEP_X} @ \
        "  " @ {Inh_SEP_Y} @ " 0.0  0.0  0.0"
    elif(diskpath == {net_EPSC_file})
	header_str = {header_str} @ {Ex_NX} @ " " @ {Ex_NY} @ " " @ {Ex_SEP_X} @ \
	"  " @ {Ex_SEP_Y} @ " 0.0  0.0	0.0"
    else
        echo "Wrong file name root!"
        header_str = ""
    end    
    return {header_str}
end

/* Create disk_out element to write netview data to a binary or ascii  file */
function make_disk_out(diskpath)
    str diskpath, filename
    if ({exists /output/{diskpath}})
            delete /output/{diskpath}
    end
    if(binary_file==0) // use asc_file
	filename = {diskpath}  @ "_" @ {RUNID} @ ".txt"
        create par_asc_file /output/{diskpath}
        setfield /output/{diskpath} leave_open 1 flush 0 filename {filename}
        setfield /output/{diskpath} float_format %.3g notime 1
        if(write_asc_header==1)
            setfield /output/{diskpath} append 1
            call /output/{diskpath} OUT_OPEN
            call /output/{diskpath} OUT_WRITE {make_header {diskpath}}
            call /output/{diskpath} RESET
        end
    else  // use disk_out to make FMT1 binary file
        filename = {diskpath}  @ "_" @ {RUNID} @ ".dat"
        create par_disk_out /output/{diskpath}
	setfield /output/{diskpath} leave_open 1 flush 0 filename {filename}
    end //if(binary_file==0)
end

function make_disk_out_messages(diskpath,srcpath, srcelement, field)
    str name, diskpath, srcpath, srcelement, field
// to be invoked by workers
    if({hflag} && {hsolve_chanmode > 1})
      foreach name ({getelementlist {srcpath}})
        raddmsg {name}/solver /output/{diskpath}@{output_node} \
          SAVE {getfield {name}/soma io_index} \
          {findsolvefield {name}/solver {name}/{srcelement} {field}}
      end
    else
      foreach name ({getelementlist {srcpath}})
         raddmsg {name}/{srcelement} /output/{diskpath}@{output_node} \
           SAVE {getfield {name}/soma io_index} {field}
        end
    end
end

function make_network_out
   str Ex_layer = "Ex_" @ {default_layer}
   if(netview_output)
      setclock 2 {netview_dt}
      make_disk_out {net_efile}
      useclock /output/{net_efile} 2
      make_disk_out {net_ifile}
      useclock /output/{net_ifile} 2
   end
   if (EPSC_output)
      make_disk_out {net_EPSC_file}
      setclock 2 {netview_dt}
      useclock /output/{net_EPSC_file} 2
   end
end

function make_network_out_messages
   str Ex_layer = "Ex_" @ {default_layer}
   str Inh_layer = "Inh_" @ {default_layer}
   if(netview_output)
      make_disk_out_messages {net_efile} /{Ex_layer}/{Ex_cell_name}[] soma Vm
      make_disk_out_messages {net_ifile} /{Inh_layer}/{Inh_cell_name}[] soma Vm
   end
   if (EPSC_output)
      make_disk_out_messages {net_EPSC_file} \
         /{Ex_layer}/{Ex_cell_name}[] {Ex_ex_synpath} Ik
   end
end

/***** functions to provide 'oddball' stimuli *****/

/* These stimulus pulse parameters are set after calling function
   make_inputs in MGBv_input2-7.g or simple_inputs2-7.g
*/
float pulse_delay, pulse_width, pulse_interval

function apply_oddball_inputs
    int ignore_steps = 0 // intials steps to ignore
    int input1 = {getfield  /oddball_inputs input1}
    int input2 = {getfield  /oddball_inputs input2}
    float prob1 = {getfield  /oddball_inputs prob1}
    int freqnum, stim
    float curr_time = {getstat -time}
    int facdep_step = {round {curr_time/pulse_interval}}
    if (use_stim_order >0)
      ignore_steps = 5
    end
    if (facdep_step < ignore_steps) // use input1 and set freqnum = 0
        setfield /MGBv[{input1}]/spikepulse level1 1.0
        setfield /MGBv[{input2}]/spikepulse level1 0.0 
      setfield /oddball_inputs freqnum 0
    else // process following steps
      if (use_stim_order == 1) // read stim from a file
          // 0 i= "oddball" input2; 1 = "frequent" input1
          stim = {readfile {stim_order_file}}
      elif (use_stim_order == 2)
        if ({facdep_step % 5} == 0)   // every 5th trial is infrequent
          stim = 0
        else
          stim = 1
        end
      elif (use_stim_order == 3)
        if ({facdep_step % 6} < 2)   // every 5th and 6th trial is infrequent
          stim = 0
        else
          stim = 1
        end
      else // use probabilities
        if ({rand 0 1} <= {prob1}) // frequent input
           stim = 1
        else
            stim = 0
        end
      end // if (use_stim_order)
      if (stim == 1) // frequent
        setfield /MGBv[{input1}]/spikepulse level1 1.0 
        setfield /MGBv[{input2}]/spikepulse level1 0.0 
        setfield /oddball_inputs freqnum 1
      else // oddball input
        setfield /MGBv[{input1}]/spikepulse level1 0.0 
        setfield /MGBv[{input2}]/spikepulse level1 1.0 
        setfield /oddball_inputs freqnum 2
      end
    end // process following steps
    if (debug_level > 1)
        echo "*** trial " {facdep_step} " stim " {stim} " prob1 " {prob1}
    end
end

function setup_oddball_inputs(input1, input2, prob1)
    int input1, input2
    float prob1 // probability of input1 - passed from global 'input1_prob'
    create script_out /oddball_inputs
    addfield  /oddball_inputs input1
    addfield  /oddball_inputs input2
    addfield  /oddball_inputs prob1
    setfield /oddball_inputs command apply_oddball_inputs
    setfield /oddball_inputs input1 {input1}
    setfield /oddball_inputs input2 {input2}
    setfield /oddball_inputs prob1 {prob1}
    // A place to store the current step input number (1 or 2, or 0 for none)
    addfield  /oddball_inputs freqnum 
    setfield  /oddball_inputs freqnum 0
    setclock 5 {pulse_interval}
    useclock /oddball_inputs 5
	  if (use_stim_order == 1)
		  openfile {stim_order_file} r
	  end
end

function make_oddball_header
    /* create a one line header with the format:

      #RUNID pulse_delay pulse_width pulse_interval sample_duration
       The sample duration and input pulse parameters are used by the
       PSD analysis tools for averaging over trial repetions of the stimulus.
    */
    float sample_duration = 0.5 // default length of sample to be taken
    str  header_str
    header_str = "#" @ {RUNID} @ "  " @ {pulse_delay}  @ "  " \
      @ {pulse_width}  @ "  " @ {pulse_interval}  @ "  " @ {sample_duration}
    return {header_str}
end

function setup_oddball_output
    /* Now create an asc_file to output the information on when the inputs
       input1 and input2 occurred.
    */
    str rootname
    rootname = "oddball_freqs"
    make_output {rootname}
    useclock {rootname} 5
    setfield  {rootname} append 1 flush 1
    call  {rootname} OUT_OPEN
    call  {rootname} OUT_WRITE {make_oddball_header}
    call  {rootname} RESET
    raddmsg@{control_node} /oddball_inputs {rootname}@{output_node} SAVE freqnum
end

/* ------------------------------------------------------------------
 Create a calculator object to hold summed Ex_ex currents, and then
 send all Ex_ex synaptic currents to it for summation

  NOTE:  Previous versions of make_EPSCsummer also summed the excitatory
  currents from the thalamic input drive. Here, the effect of any input
  drive would be indirect through propagation of its effect from cell to
  cell via synaptic connections.
------------------------------------------------------------------- */

function do_EPSCsum_out // called from output_node
    str data_source = "/EPSCsummer" // data_source sums Ex_cell ex currents
    create calculator {data_source}
    useclock {data_source} 1 // the clock for out_dt
    make_output {EPSC_sum_file}
    // should be OK to add the message here
        addmsg {data_source} {EPSC_sum_file} SAVE output
end

function make_EPSCsummer_messages // data_source sums Ex_cell ex currents
    int i
    str data_source = "/EPSCsummer" // data_source sums Ex_cell ex currents
    str Ex_layer = "Ex_" @ {default_layer}
    for (i=0; i < Ex_cells_per_slice; i = i + 1)
        if({hflag} && {hsolve_chanmode > 1})
            raddmsg /{Ex_layer}/{Ex_cell_name}[{i}]/solver {data_source}@{output_node} \
              SUM {findsolvefield /{Ex_layer}/{Ex_cell_name}[{i}]/solver \
              {Ex_ex_synpath} Ik}
        else
            raddmsg /{Ex_layer}/{Ex_cell_name}[{i}]/{Ex_ex_synpath} \
                {data_source}@{output_node} SUM Ik
        end
    end
end // make_EPSCsummer_messages

/*****  Functions for post-run output *****/

function do_run_summary
    str filename = {sum_file} @ "_" @ {RUNID} @ ".txt"
    openfile {filename} w
    writefile {filename} "Script:" {script_name} "  RUNID:" {RUNID} "  seed:" {seed} \
        "  date:" {getdate}
    writefile {filename} "tmax:" {tmax} " dt:" {dt} " out_dt:" {out_dt} \
        " netview_dt:" {netview_dt} 
    writefile {filename} "EPSC_output:" {EPSC_output} "  netview_output:" {netview_output}
    writefile {filename} "Ex_NX:" {Ex_NX} " Ex_NY:" {Ex_NY} " Inh_NX:" \
        {Inh_NX} "  Inh_NY:" {Inh_NY}
    writefile {filename} "Ex_SEP_X:" {Ex_SEP_X} " Ex_SEP_Y:" {Ex_SEP_Y} \
        "  Inh_SEP_X:" {Inh_SEP_X}  "  Inh_SEP_Y:" {Inh_SEP_Y}
    writefile {filename} "Ninputs:" {Ninputs} " bg Ex freq:" {frequency} \
        " bg Ex gmax:" {Ex_bg_gmax}
    writefile {filename} "================== Network parameters ================="
    writefile {filename} "Ex_ex_gmax:" {Ex_ex_gmax} " Ex_inh_gmax:" {Ex_inh_gmax} \
        "  Inh_ex_gmax:" {Inh_ex_gmax} "  Inh_inh_gmax:" {Inh_inh_gmax}
    writefile {filename} "tau1_ex:" {tau1_ex} "  tau2_ex:" {tau2_ex}  \
        "  tau1_inh:" {tau1_inh} "  tau2_inh:" {tau2_inh}
    writefile {filename} "tau1_Inh_ex:" {tau1_Inh_ex} "  tau2_Inh_ex:" {tau2_Inh_ex}
    writefile {filename} "syn_weight: " {syn_weight} \
        " prop_delay:" {prop_delay} " default drive weight:" {drive_weight}
    writefile {filename} "Ex_drive_gmax: " {Ex_drive_gmax} \
        "  Inh_drive_gmax: " {Inh_drive_gmax}
    /* The code below is specific to the input model -- see MGBv_input.g */
    writefile {filename} "MGBv input delay: " {input_delay} \
	"  MGBv input jitter :" {input_jitter}
    writefile {filename} "input_spread : " {input_spread} \
        "  connect_network :" {connect_network} "use_weight_decay: " \
        {use_weight_decay}
    writefile {filename} "change_weights: " {change_weights} \
         " use_Ex_ex_FACDEP: " {use_Ex_ex_FACDEP} \
         " use_Inh_ex_FACDEP: " {use_Inh_ex_FACDEP}
    writefile {filename} "input paradigm: "  {input_paradigm} "input1_prob: " {input1_prob}

    writefile {filename} "================== Thalamic inputs ================="
    writefile {filename} "Input" " Row " "Frequency" "    State" "   Weight" \
        "    Delay" "    Width" " Period"
    int input_num, input_row
    str  input_source = "/MGBv" // Name of the array of input elements
    float input_freq, delay, width, interval, drive_weight
    str pulse_src, spike_out
    str toggle_state
    floatformat %10.3f
    for (input_num = 1; input_num <= {Ninputs}; input_num= input_num +1)
        pulse_src = {input_source} @ "[" @ {input_num} @ "]" @ "/spikepulse"
        spike_out = {input_source} @ "[" @ {input_num} @ "]" @ "/soma/spike"
        input_row  = \
          {getfield {{input_source} @ "[" @ {input_num} @ "]"} input_row}
        input_freq  = \
          {getfield {{input_source} @ "[" @ {input_num} @ "]"} input_freq}
        drive_weight = \
          {getfield {{input_source} @ "[" @ {input_num} @ "]"} output_weight}
        delay = {getfield {pulse_src} delay1 }
        width = {getfield {pulse_src} width1}
        interval = {getfield {pulse_src} delay1} + {getfield {pulse_src} delay2}
        // get the spiketoggle state
        toggle_state = "OFF"
        if ({getfield {pulse_src} level1} > 0.5)
            toggle_state = "ON"
        end
        writefile {filename} {input_num} {input_row} -n -format "%5s"
        writefile {filename} {input_freq} {toggle_state} {drive_weight} \
            {delay} {width} {interval}  -format %10s
    end
    writefile {filename} "----------------------------------------------------------"
    writefile {filename} "Notes:"
    closefile {filename}
    floatformat %0.10g
end

function change_RUNID(value)
    str value, path, new_name, rootname
    RUNID =  value
    // Set up new file names for output
    if (i_am_output_node)
        foreach path ({el /##[OBJECT=asc_file],/##[OBJECT=par_asc_file]})
            rootname = {getfield {path} name} 
            new_name = {rootname} @ "_" @ {RUNID} @ ".txt"
            if (debug_level)
                echo@{control_node} "New filename: " {new_name} "path: " {path}
            end
            setfield {path} append 0 leave_open 0
            call {path} RESET
            setfield {path} filename {new_name}
            setfield {path} append 1 leave_open 1
            if ((rootname == {net_efile}) || (rootname == {net_ifile}) \
              || (rootname == {net_EPSC_file}))
                call {path} OUT_OPEN
                call {path} OUT_WRITE {make_header {rootname}}
            elif (rootname == "oddball_freqs")
                call {path} OUT_OPEN
                call {path} OUT_WRITE {make_oddball_header}
		call {path} RESET
            end
        end
    end // (i_am_output_node)
end

extern reset_freqmons // defined below
function do_reset // from control_node
    if(calc_firing_freqs) // not used here
        reset_freqmons
    end
    reset@all
end

function step_tmax
    echo "dt = "{getclock 0}"   tmax = "{tmax}
    echo "RUNID: " {RUNID}
    echo "START: " {getdate}
    step@all {tmax} -time
    // step@all // do one more step be sure tmax is included
    echo "END  : " {getdate}
    do_run_summary
end

/*** functions to set up FACDEP *****/

function setup_depdep_Ex_ex
    str facdep_applier = "depdep_Ex_ex"
    str cellpath = "/Ex_" @ {default_layer} @ "/" @ {Ex_cell_name} @ "[]"
    str relsynpath = {Ex_ex_synpath}
    if (change_weights < 0) // Don't even make the object
        return
    end
    if (!{exists {facdep_applier}})
            create facdep_rules2 {facdep_applier}
    end
    setfield {facdep_applier} cellpath {cellpath}
    setfield {facdep_applier} synpath {relsynpath}
    setfield {facdep_applier} change_weights {change_weights}
    setfield {facdep_applier} use_depdep 1 // two depression factors
    setfield {facdep_applier} debug_level 0 // {debug_level}
    setfield {facdep_applier} dD1 {Ex_dD1}
    setfield {facdep_applier} dD2 {Ex_dD2}
    setfield {facdep_applier} tau_D1 {Ex_tau_D1}
    setfield {facdep_applier} tau_D2 {Ex_tau_D2}
    setclock 3 {facdep_dt}
    useclock  {facdep_applier} 3
end // function setup_depdep_Ex_ex

function setup_depdep_Inh_ex
    str facdep_applier = "depdep_Inh_ex"
    str cellpath = "Inh_" @ {default_layer} @ "/" @ {Inh_cell_name} @ "[]"
    str relsynpath = {Inh_ex_synpath}
    if (change_weights < 0) // Don't even make the object
        return
    end
    if (!{exists {facdep_applier}})
            create facdep_rules2 {facdep_applier}
    end
    setfield {facdep_applier} cellpath {cellpath}
    setfield {facdep_applier} synpath {relsynpath}
    setfield {facdep_applier} change_weights {change_weights}
    setfield {facdep_applier} use_depdep 1 // two depression factors
    setfield {facdep_applier} debug_level 0 // {debug_level}
    setfield {facdep_applier} dD1 {Inh_dD1}
    setfield {facdep_applier} dD2 {Inh_dD2}
    setfield {facdep_applier} tau_D1 {Inh_tau_D1}
    setfield {facdep_applier} tau_D2 {Inh_tau_D2}
    setclock 3 {facdep_dt}
    useclock  {facdep_applier} 3
end // function setup_depdep_Inh_ex

function make_pyrVm_out(cell_num) // execute from output_node
    int cell_num
    str rootname = {pyrVm_file} @ {cell_num} // string concat, not a node
    // {pyrVm_file} is defined in main script as "pyrVm"
    make_output {rootname}
end

function make_pyrVm_out_messages(cell_num) // execute from node with cell_num 
    int cell_num, index, node
    str Ex_layer = "Ex_" @ {default_layer}
    str rootname = {pyrVm_file} @ {cell_num}
    // {pyrVm_file} is defined in main script as "pyrVm"
    index = {cell_index {cell_num} {Ex_cells_per_slice}}
    node = {cell_node {cell_num} {Ex_cells_per_slice}}
    if ({mynode} != node)
        return
    end
    if({hflag} && {hsolve_chanmode > 1})
        raddmsg /{Ex_layer}/{Ex_cell_name}[{index}]/solver \
          {rootname}@{output_node} SAVE \
          {findsolvefield /{Ex_layer}/{Ex_cell_name}[{index}]/solver \
          soma Vm}
    else
        raddmsg /{Ex_layer}/{Ex_cell_name}[{index}]/soma \
	  {rootname}@{output_node} SAVE  Vm
    end
end // function make_pyrVm_out_messages

function print_Ex_synapse_info(cell_num) // execute from node with cell_num
    int cell_num, index, node
    str Ex_layer = "Ex_" @ {default_layer}
    index = {cell_index {cell_num} {Ex_cells_per_slice}}
    node = {cell_node {cell_num} {Ex_cells_per_slice}}
    if ({mynode} == node)
        synapse_info /{Ex_layer}/{Ex_cell_name}[{index}]/{Ex_ex_synpath}
        synapse_info /{Ex_layer}/{Ex_cell_name}[{index}]/{Ex_inh_synpath}
    end
end

function print_Inh_synapse_info(cell_num) // execute from node with cell_num
    int cell_num, index, node
    str Inh_layer = "Inh_" @ {default_layer}
    index = {cell_index {cell_num} {Inh_cells_per_slice}}
    node = {cell_node {cell_num} {Inh_cells_per_slice}}
    if ({mynode} == node)
        synapse_info /{Inh_layer}/{Inh_cell_name}[{index}]/{Inh_ex_synpath}
    end
end

function reseed_RNG(seed)
    int seed, newseed
    if ({mynode} == 0) // control node (and worker0 if = 0) gets orig seed
        newseed = seed
     else // worker 1 and higher get multiples of 347 added to seed
	newseed = seed + ({mynode} - worker0)*347
     end
     echo@{control_node} "New RNG seed = " {newseed}
     randseed {newseed}
end

// Flush the buffers of all text output elements
function flush_buffers
    str path
    if (i_am_output_node)
        foreach path ({el /##[OBJECT=asc_file],/##[OBJECT=par_asc_file]})
            if(pdebug)
                echo@{control_node} "Final step and flush for: " {path}
            end
	    setfield {path} flush 1
        end
        step@all {netview_dt} -time // add extra steps to flush buffers
    end //  (i_am_output_node)
end
