 //genesis - LFPfuncs.g

/* Functions for calculating local field potentials, EEG, and ECoG
   potentials. Uses an efield2 object and can be used with hsolve.

  Will be included in the main script if calc_LFP = 1.

  The global string variable 'solvepath' is defined in the main script
  to give the name that is being used for the hsolve element, e.g. "solver".

  Output clock will be clock 1 (out_dt)

  function make_electrode(name, x, y, z, cellpath) sums over all pyr cells,
  is most often used to sum over all pyrmidal cells.
*/

/* When recording field potentials from many cells, it may be
   useful to record only synaptic currents rather than include
   transmembrane currents generated by action potential, leak,
   and capacitive currents.

   The global variable 'use_syn_currents' is used to allow this option.

   set in main script:
   int use_syn_currents = 1 // calculate the efield using only synchan currents

   Typical usage:
   make_electrode "Electrode_1" 0.0 0.0 0.001 {cellpath} // Lower Left: cell [0]
   make_electrode "Electrode_2" 0.00096 0.00064 0.001 {cellpath} // over pyr_4[792]

   To include all  currents, including capacitive currents, one can use
   make_electrode without hsolve and use_syn_currents = 0. This will use
   the axial currents Im, which are  equal to the total membrane and
   capacitive Cm*(dCvm/dt) currents. The hsolve does not include the
   capacitive currents in the calculation of Im.

   Alternatively, use hsolve and the function make_electrode_cells.
   This creates an efield_cells object that correctly calculates Im
   and the resulting field. It also places the Im value in the disabled
   object that was taken over by the solver, where it is accessible for use
   by efield or efield2 as if it were not hsolved.

*/

// This function uses efield2 to correctly set the distance to hsolved elements
// It is to be executed by output_node
function make_electrode(name, x, y, z) // sums over all cells in path
    str name        // electrode name
    float x, y, z   // electrode position
    if (calc_LFP == 0)
         return
    end
    if (!{exists /electrodes})
        create neutral /electrodes
    end
    pushe /electrodes
    if ({exists {name}})
       delete {name}
    end
    float scale // scale factor for currents
    if (use_syn_currents)
       scale = -0.32 // invert sign of trans-membrane currents
    else // Use Im to give axial currents
        scale = 0.32 // positive - axial currents are positive-in
    end
    create efield2 {name}
    setfield {name} x {x} y {y} z {z}
    setfield {name} scale {scale}
    pope
    useclock /electrodes/{name} 1 // clock at out_dt
    
// In order to output the results, use function do_electrodes_out
end // function make_electrode

function make_electrode_messages(name, cellpath) // executed by worker nodes
    str name        // electrode name
    str cellpath    // path to cells producing the field
    float x, y, z   // electrode position
    x = {getfield@{output_node} /electrodes/{name} x}
    y = {getfield@{output_node} /electrodes/{name} y}
    z = {getfield@{output_node} /electrodes/{name} z}
    float cx, cy, cz   // channel position
    float distance
    str compt, chan, src
    int i
    for (i = 0; i < Ex_cells_per_slice; i = i + 1)
        src = {cellpath} @ "[" @ {i} @ "]"
        if({hflag} && {hsolve_chanmode > 1})
            foreach compt ({el {src}/#[OBJECT=compartment]})
                if (use_syn_currents)
                  foreach chan  ({el {compt}/#[OBJECT=synchan]})
                    cx = {getfield {chan} x}
                    cy = {getfield {chan} y}
                    cz= {getfield {chan} z}
                    distance = {sqrt { {pow {{cx} - {x}} 2.0} + \
                        {pow {{cy} - {y}} 2.0} + \
                        {pow {{cz} - {z}} 2.0} }}
                      if(debug_level > 1)
                        echo {getfield {compt} name}  {getfield {chan} name} \
                        "distance = " {distance}
                      end
                      raddmsg {src}/{solvepath} /electrodes/{name}@{output_node} \
                        CURRENT {findsolvefield {src}/{{solvepath}} {chan} Ik} \
                        {distance}
                  end // foreach
                else // Use Im to give axial currents
                  /* Assume that efield_cells has updated hsolved compt Im */
                  raddmsg {compt} /electrodes/{name}@{output_node} CURRENT Im 0.0
                end
            end // loop over compts
        else  // not hsolve - Ileak is not calculated
            foreach compt ({el {src}/#[OBJECT=compartment]})
                if(use_syn_currents)
                   foreach chan  ({el {compt}/#[OBJECT=synchan]})
                       raddmsg {src}/{synpath}  /electrodes/{name}@{output_node} Ik 0.0
                   end // foreach chan
                else
                   raddmsg {compt}  /electrodes/{name}@{output_node} CURRENT Im 0.0
                end
            end // loop over compts
        end // if({hflag} && {hsolve_chanmode > 1})
      end // loop over cells
end

/* A variation that calculates the Efield for only a specified synpath
   e.g. "apical3/AMPA_pyr", or {Ex_ex_synpath}
*/

function make_electrode_synpath(name, x, y, z)
    str name        // electrode name
    str synpath     // path to compartment/synchan
    float x, y, z   // electrode position
    str cellpath    // path to cells
    str chan        // full path to the particular synchan
    float cx, cy, cz   // channel position
    float distance
    if (calc_LFP == 0)
         return
    end
    if (!{exists /electrodes})
        create neutral /electrodes
    end
    pushe /electrodes
    if ({exists {name}})
       delete {name}
    end
    create efield2 {name}
    setfield {name} x {x} y {y} z {z}
    setfield {name} scale -0.32 // membrane currents should be 'negative-in'
    pope

    useclock /electrodes/{name} 1 // clock at out_dt
    call /electrodes/{name} RECALC // calculate the distances to the compartments
    // In order to output the results, use function do_electrodes_out
end // function make_electrode_synpath

function make_electrode_synpath_messages(name, cellpath, synpath)
    str name        // electrode name
    str synpath     // path to compartment/synchan
    float x, y, z   // electrode position
    str cellpath    // path to cells
    str chan        // full path to the particular synchan
    x = {getfield@{output_node} /electrodes/{name} x}
    y = {getfield@{output_node} /electrodes/{name} y}
    z = {getfield@{output_node} /electrodes/{name} z}
    float cx, cy, cz   // channel position
    float distance
    str src
    int i
    for (i = 0; i < {Ex_cells_per_slice}; i = i + 1)
        src = {cellpath} @ "[" @ {i} @ "]"
        if({hflag} && {hsolve_chanmode > 1})
        // need to provide the distance to the electrode
            chan = {src} @ "/" @ {synpath}
            cx = {getfield {chan} x}
            cy = {getfield {chan} y}
            cz= {getfield {chan} z}
            distance = {sqrt { {pow {{cx} - {x}} 2.0} + \
                        {pow {{cy} - {y}} 2.0} + \
                        {pow {{cz} - {z}} 2.0} }}
            raddmsg {src}/{solvepath}  /electrodes/{name}@{output_node} CURRENT \
              {findsolvefield {src}/{solvepath} {synpath} Ik} {distance}
        else  // not hsolve
            raddmsg {src}/{synpath}  /electrodes/{name}@{output_node} CURRENT Ik 0.0
        end // if({hflag} && {hsolve_chanmode > 1})
    end // loop over cells
end

function make_electrode_cells(name, cellpath, x, y, z)
    str name        // electrode name
    str  cellpath // here cellpath is {cellpath}[]
    str worker_name // name used for worker electrodes
    /* Example:
       make_electrode_cells "Electrode_cells_Im" {cellpath}[]  \
       {(Ex_NX-1)*Ex_SEP_X/2.0} {39.5*Ex_SEP_Y} 0.001 // use efield_cells
    */
    float x, y, z   // electrode position
    if (calc_LFP == 0)
         return
    end
    if (i_am_worker_node)
        if (!{exists /electrodes})
            create neutral /electrodes
        end
        worker_name = {name} @ "_" @ {mynode}
        pushe /electrodes
            if ({exists {worker_name}})
                delete {worker_name}
            end
        create efield_cells {worker_name}
        setfield {worker_name} x {x} y {y} z {z}
        setfield {worker_name} cellpath {cellpath}
        setfield {worker_name} solvepath {solvepath} // global variable
        setfield {worker_name} scale 0.32
        setfield {worker_name} debug_level 0
        pope
        useclock /electrodes/{worker_name} 1 // clock at out_dt
    end // if (i_am_worker_node)
    if (i_am_output_node)
        if (!{exists /electrodes})
            create neutral /electrodes
        end
        create calculator /electrodes/{name}
        useclock /electrodes/{name} 1 // clock at out_dt
    end
end // function make_electrode_cells

function make_electrode_cells_messages(name)
    str name, worker_name
    if (i_am_worker_node)
        worker_name = {name} @ "_" @ {mynode}
        raddmsg /electrodes/{worker_name} /electrodes/{name}@{output_node} SUM field
    end
end

function do_electrodes_out
    str name, path
    foreach path ({el /electrodes/#[OBJECT=efield],/electrodes/#[OBJECT=efield2]})

/* create an asc_file to output the LFP to a file
           {name} @ "_" @ {RUNID} @ ".txt". Here, the asc_file element
            will be  /{name} and the electrode will be /electrodes/{name}
        */
        name = {getpath {path} -tail}
        make_output {name}
        addmsg /electrodes/{name}  {name} SAVE field
    end
end

function do_electrodes_cells_out
    str name, path
    foreach path ({el /electrodes/#[OBJECT=calculator]})
        name = {getpath {path} -tail}
        make_output {name}
        addmsg /electrodes/{name}  {name} SAVE output
    end
end


/********** graphics functions *********/

extern makegraphscale
function make_LFPplot
    str graph_form = "/electrode_graphs"
    str name, path, title_str
    float x, y, z
    floatformat %0.4g // limit precision of displayed x, y, z
    create xform {graph_form} [280,0,400,800]
    pushe {graph_form}
      foreach path ({el /electrodes/#[OBJECT=efield_cells],/electrodes/#[OBJECT=efield],/electrodes/#[OBJECT=efield2]})
       name = {getpath {path} -tail}
       x = {getfield {path} x}; y = {getfield {path} y}; z = {getfield {path} z}
       title_str = {name} @ " at " @ {x} @ " " @ {y} @ " " @ {z}
       // make the graph
       create xgraph {name} -hgeom 33% -title  {title_str}
       setfield {name} XUnits sec YUnits V
       setfield {name} xmax {tmax}   ymin -0.00030  ymax 0.000005
    end // making electrode graphs

    floatformat %0.10g // restore default precision
    if (calc_EPSCsum)
        create xgraph EPSCsum -hgeom 25% -title "Summed Ex_ex Currents"
        setfield EPSCsum XUnits sec YUnits Amp
        setfield EPSCsum xmax {tmax}   ymin 0 ymax 3e-6
    end
    pope
    // e.g. makegraphscale {graph_form}/Electrode_3
    foreach path ({el /electrodes/#[OBJECT=efield_cells],/electrodes/#[OBJECT=efield],/electrodes/#[OBJECT=efield2]})
      
       name = {getpath {path} -tail}
       makegraphscale {graph_form}/{name}
    end
    if (calc_EPSCsum)
        makegraphscale {graph_form}/EPSCsum
    end
    useclock  {graph_form} 2
    xshow {graph_form}
end

extern colors

function make_LFPplot_messages 
    str name, path, col_name
    int col_num = 0
    str graph_form = "/electrode_graphs"
    foreach path ({el /electrodes/#[OBJECT=efield_cells],/electrodes/#[OBJECT=efield],/electrodes/#[OBJECT=efield2]})
       name = {getpath {path} -tail}
       col_name = {colors {col_num}}
       addmsg /electrodes/{name} {graph_form}/{name} \
         PLOT field *{name} *{col_name}
       col_num = col_num + 1
    end // loop over electrodes
    if (calc_EPSCsum)
      addmsg {data_source} {graph_form}/EPSCsum \
        PLOT output *EPSC_sum *magenta
    end
end
