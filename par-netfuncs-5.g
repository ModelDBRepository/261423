//genesis

// /**********************************************************************
// ** This simulation script and the files included in this package
// ** are Copyright (C) 2013 by David Beeman (dbeeman@colorado.edu)
// ** and are made available under the terms of the
// ** GNU Lesser General Public License version 2.1
// ** See the file copying.txt for the full notice.
// *********************************************************************/

// /*============================================================================

//  Network setup function definitions for the "par_ACnet2" series of simulations

//  ============================================================================*/
/* These functions assume that network has been divided horizontal
   slices that are on nodes 0 through  n_slices-1.  Each slice has
   cells_per_slice = NX*NY / n_slices cells.
*/
function cell_node(cell_num, cells_per_slice) // return the node that cell_num is in
    int cell_num, cells_per_slice
    return {trunc {cell_num/cells_per_slice}} + {worker0}
end
// return the cell index within that node
function cell_index(cell_num, cells_per_slice)
    int cell_num, node, cells_per_slice
    node = {cell_node {cell_num} {cells_per_slice}}
    return {cell_num - (node - worker0)*cells_per_slice}
end

// /***** functions to connect network *****/ 
// see function connect_cells below

/***** functions to set weights and delays *****/

function set_weights(weight)
    float weight
    str Ex_layer = "Ex_" @ {default_layer}
    str Inh_layer = "Inh_" @ {default_layer}
    syn_weight = weight  // set the global variable to the new weight
    // In these simulations, use_weight_decay = 0
    if (!use_weight_decay) // use fixed weights
        // use fixed weights - the default
        rvolumeweight /{Ex_layer}/{Ex_cell_name}[]/soma/spike \
            /{Ex_layer}/{Ex_cell_name}[]/{Ex_ex_synpath} \
              -fixed {syn_weight}
        rvolumeweight /{Ex_layer}/{Ex_cell_name}[]/soma/spike \
            /{Inh_layer}/{Inh_cell_name}[]/{Inh_ex_synpath} \
              -fixed {syn_weight}
        rvolumeweight /{Inh_layer}/{Inh_cell_name}[]/soma/spike  \
            /{Ex_layer}/{Ex_cell_name}[]/{Ex_inh_synpath} \
              -fixed {syn_weight}
        rvolumeweight /{Inh_layer}/{Inh_cell_name}[]/soma/spike  \
            /{Inh_layer}/{Inh_cell_name}[]/{Inh_inh_synpath} \
              -fixed {syn_weight}

        echo "All maxium synaptic weights set to "{weight}
        /* Newer versions of this script use the test:
             if(use_Ex_ex_FACDEP && rand_weight)
        with rand_weight = 0, so that the following randomization of
        Ex_ex weights is not performed. The rand_weight option was intended
        for use with stdp_rules for modeling LTP and LTD. Here the weights
        are randomized, to reproduce the published simulation results. Thie
        has a small, but observable effect on the recovery from the N! peak.
        */
        // Make an exception for Ex_cell excitatory (Ex_ex) weights
        if(use_Ex_ex_FACDEP)
          int n, i
          /* reseed the RNG to give the same weights
             This is necessary only with the option "rand_weight = 1"
             to randomly initialize the Ex_ex weights
          */
          randseed {seed} // Better to use reseed_RNG {seed}
          for (n=0; n < Ex_cells_per_slice; n=n+1)
            str chanpath =  "/" @ {Ex_layer} @"/" \
             @ {Ex_cell_name} @  "[" @ {n} @ "]/" @ {Ex_ex_synpath}
            for (i=0; i <  {getfield {chanpath} nsynapses}; i=i+1)
              if (rand_weight == 1)
                  weight = {rand 0 {2*avg_weight}}
              end
              setfield {chanpath} synapse[{i}].weight {weight}
            end // for (i=0
          end // for (n=0
          if (debug_level)
         echo@{control_node} "All Ex_ex synaptic weights set to average of "{avg_weight}
         echo@{control_node} "All other synaptic weights set to "{syn_weight}
           end
         end // if(use_Ex_ex_FACDEP && rand_weight)
    end //   if (!use_weight_decay)
end


/* If the delay is zero, set it as the fixed delay, else use the conduction
   velocity to calculate delays based on radial distance to the target.
   For a fixed delay, either planardelay or volumedelay can be used.
   As axonal conduction velocities are the same in both the vertical
   direction, and in the horizontal plane, volumedelay should be
   used to account for the vertical distance traveled to the dendrites.
*/

function set_delays(delay)
    float delay
    str Ex_layer = "Ex_" @ {default_layer}
    str Inh_layer = "Inh_" @ {default_layer}
    prop_delay = delay
    if (delay == 0.0)
        rvolumerdelay /{Ex_layer}/{Ex_cell_name}[]/soma/spike -fixed {delay}
        rvolumedelay /{Inh_layer}/{Inh_cell_name}[]/soma/spike -fixed {delay}
    else
        rvolumedelay /{Ex_layer}/{Ex_cell_name}[]/soma/spike -radial {1/delay}
        rvolumedelay /{Inh_layer}/{Inh_cell_name}[]/soma/spike -radial {1/delay}
    end
    if (debug_level)
        echo@{control_node}  "All propagation delays set to "{delay}" sec/m"
    end
end

//=============================================================
//    Functions to set up the network
//=============================================================

function make_prototypes
  /* Step 1: Assemble the components to build the prototype cell under the
     neutral element /library.  This is done by including "prododefs.g"
     before using the function make_prototypes.
  */
  
  /* Step 2: Create the prototype cell specified in 'cellfile', using readcell.
     This should set up the apropriate synchans in specified compartments,
     with a spikegen element "spike" attached to the soma.  This will be
     done in /library, where it will be available to be copied into a network

     In this case there are two types of cells "Ex_cell", "Inh_cell",

  */

  readcell {Ex_cellfile} /library/{Ex_cell_name}
  readcell {Inh_cellfile} /library/{Inh_cell_name}
  addfield /library/{Ex_cell_name}/soma io_index // for ordering i/o messages
  addfield /library/{Inh_cell_name}/soma io_index // for ordering i/o messages

  // In this case, use different values from the defaults in 'cellfile'
  setfield /library/{Ex_cell_name}/{Ex_ex_synpath} gmax {Ex_ex_gmax}
  setfield /library/{Ex_cell_name}/{Ex_inh_synpath} gmax {Ex_inh_gmax}
  setfield /library/{Ex_cell_name}/{Ex_drive_synpath} gmax {Ex_ex_gmax}
  setfield /library/{Ex_cell_name}/soma/spike thresh 0 abs_refract 0.002 \
	output_amp 1
  // Note: the Ex-cell has wider and slower spikes, so use larger abs_refract

  setfield /library/{Inh_cell_name}/{Inh_ex_synpath} gmax {Inh_ex_gmax}
  setfield /library/{Inh_cell_name}/{Inh_inh_synpath} gmax {Inh_inh_gmax}
  setfield /library/{Inh_cell_name}/{Inh_drive_synpath} gmax {Inh_ex_gmax}
  setfield /library/{Inh_cell_name}/soma/spike thresh 0  abs_refract 0.001 \
	 output_amp 1
end

/***** functions to set up the network *****/

/* Setting up hsolve for a network requires setting up a solver for
   one cell of each type in the network and then duplicating the
   solvers.  The procedure is described in the advanced tutorial
   'Simulations with GENESIS using hsolve by Hugo Cornelis' from
   genesis-sim.org/GENESIS/UGTD/Tutorials/advanced-tutorials
*/
function make_hsolve(layer)
    if(!hflag)
        return
    end
    str layer
    // do for Ex_layer
    str netlayer = "/Ex_" @ {layer}
    str cellname = {netlayer} @ "/" @ {Ex_cell_name}
    pushe {cellname}[0]
    create hsolve {solvepath}
    setmethod solver 11 // Use Crank-Nicholson
    setfield {solvepath} chanmode {hsolve_chanmode} \
        path "../[][TYPE=compartment]"
    setfield {solvepath} computeIm 1
    call {solvepath} SETUP
    int i
    for (i = 1 ; i < Ex_cells_per_slice ; i = i + 1)
        call {solvepath} DUPLICATE \
          {cellname}[{i}]/{solvepath}  ../##[][TYPE=compartment]
        setfield {cellname}[{i}]/{solvepath} \
            x {getfield {cellname}[{i}] x} \
            y {getfield {cellname}[{i}] y} \
            z {getfield {cellname}[{i}] z}
    end
    pope
    str netlayer = "/Inh_" @ {layer}
    str cellname = {netlayer} @ "/" @ {Inh_cell_name}
    pushe {cellname}[0]
    create hsolve {solvepath}
    setmethod {solvepath} 11 // Crank-Nicholson
    setfield {solvepath} chanmode  {hsolve_chanmode} path "../[][TYPE=compartment]"
    call {solvepath} SETUP
    int i
    for (i = 1 ; i < Inh_cells_per_slice ; i = i + 1)
        call {solvepath} DUPLICATE \
            {cellname}[{i}]/{solvepath}  ../##[][TYPE=compartment]
        setfield {cellname}[{i}]/{solvepath} \
            x {getfield {cellname}[{i}] x} \
            y {getfield {cellname}[{i}] y} \
            z {getfield {cellname}[{i}] z}
    end
    pope
    if(pdebug)
        echo@{control_node} "Set up hsolve on node " {mynode}
    end
end

function make_network_slice(layer)
    /*
       For the unsliced network, there will be NX cells along the
       x-direction, separated by SEP_X, and NY cells along the y-direction,
       separated by SEP_Y.  The default origin is (0, 0).  This will be the
       coordinates of cell[0].  The last cell, cell[{NX*NY-1}], will be at
       (NX*SEP_X -1, NY*SEP_Y-1).

       For the parallel version, each worker node makes a slice of the net
       with the origin set at the lower left corner of the slice, and a
       soma io_index that represents the cell number on the unsliced net.
       Each slice node will have a set of cells cell[0] -
       cell[{cells_per_slice - 1}], but the soma io_index will be in the
       appropriate part of the range from 0 - NX*NX-1, in order to
       identify its location on the full net.

       Here, each network layer is divided into horizontal slices of dimension
       NX x NY/n_slices.
    */
    int i, slice
    slice = {mynode} - {worker0}  // The slice numbers start with 0
    str layer
    str Ex_layer = "Ex_" @ {layer}
    str Inh_layer = "Inh_" @ {layer}

    createmap /library/{Ex_cell_name} /{Ex_layer} {Ex_NX} {Ex_NY/n_slices} \
      -delta {Ex_SEP_X} {Ex_SEP_Y} \
      -origin 0 {slice * Ex_SEP_Y * Ex_NY / n_slices}
        // Label the cell soma with the io_index. For horizontal slices,
        // this is equal to cell_num on the unsliced network.
        for (i = 0; i < Ex_cells_per_slice; i=i+1)
            setfield  /{Ex_layer}/{Ex_cell_name}[{i}]/soma \
                     io_index {slice * Ex_cells_per_slice + i}
        end


    // Displace the /Inh_layer4 origin to be in between Ex_cells
    createmap /library/{Inh_cell_name} /{Inh_layer} {Inh_NX} {Inh_NY/n_slices} \
      -delta {Inh_SEP_X} {Inh_SEP_Y} \
      -origin  {Ex_SEP_X/2}  {Ex_SEP_Y/2 + slice * Inh_SEP_Y * Inh_NY / n_slices}
        // Label the cell soma with the io_index. For horizontal slices,
        // this is equal to cell_num on the unsliced network.
        for (i = 0; i < Inh_cells_per_slice; i=i+1)
            setfield  /{Inh_layer}/{Inh_cell_name}[{i}]/soma \
                     io_index {slice * Inh_cells_per_slice + i}
        end
end // function make_network_slice

/***** functions to connect network *****/
/* Use function planarconnect_probdecay when volumeconnect3 is not available */

/* planarconnect_probdecay is a specialized version of planarconnect where
   the connection probability is weighted with a probability
   prob0*exp(-r*r/(sigma*sigma)) where r is the radial distance
   between source and destination in the x-y plane.

   Unlike the more general version, it is assumed that all cells in the
   source network are used as sources.  The destination region is defined
   as a ring between a radius rmin and rmax in the destination network.
   The sources are all assumed to be spikegen objects, and the destinations
   to be synchans, or variants such as facsynchans.
   
   An example of the ith cell in the source network would be:

   {src_cell_path}[{i}]{src_spikepath}

   with
   src_cell_path = /Ex_layer/{Ex_cell_name}
   src_spikepath = "soma/spike"

   An example destination would use:

   dst_cell_path = /Inh_layer/{Inh_cell_name}
   dst_synpath = {Inh_ex_synpath}

   For example:

   planarconnect_probdecay /Ex_layer4/{Ex_cell_name} "soma/spike" {Ex_NX*Ex_NY} \
      /Inh_layer4/{Inh_cell_name} {Inh_ex_synpath} {Inh_NX*Inh_NY} \
      0.0 {pyr_range} {Ex_Inh_prob} {Ex_sigma}

  For the parallel version, n_src is the number of cells in a slice,
  Ex_cells_per_slice or Inh_cells_per_slice. But, n_dstis the total number
  of destination cells across all slices, Ex_NX*Ex_NY or Inh_NX*Inh_NY.
*/

function planarconnect_probdecay(src_cell_path, src_spikepath, n_src, \
    dst_cell_path, dst_synpath, n_dst, rmin, rmax, prob0, sigma)

    str src, src_cell_path, src_spikepath, dst, dst_cell_path, dst_synpath
    int n_src, n_dst, i, j, dst_node, dst_index
    float rmin, rmax, prob0, sigma, x, y, dst_x, dst_y, prob, rsqr

    float decay = 1.0/(sigma*sigma)
    float rminsqr = rmin*rmin
    float rmaxsqr = rmax*rmax
    int cells_per_slice = n_dst/n_slices
    /* loop over all source cells - this will be all cells in layer */
    for (i = 0 ; i < {n_src} ; i = i + 1)

       src = {src_cell_path} @ "[" @ {i} @ "]/" @ {src_spikepath}
       x = {getfield {src} x}
       y = {getfield {src} y}

        /* loop over all possible destination cells - all in dest layer */
        for (j = 0 ; j < {n_dst} ; j = j + 1)
	    if (n_dst > Inh_NX*Inh_NY) // The destination is an Ex_layer
                dst_x = (j%Ex_NX)*Ex_SEP_X
                dst_y = ({trunc {j/Ex_NX}})*Ex_SEP_Y
            else // It is an Inh_layer
                dst_x = (j%Inh_NX)*Inh_SEP_X + Ex_SEP_X/2
                dst_y = ({trunc {j/Inh_NX}})*Inh_SEP_Y + Ex_SEP_Y/2
            end
            rsqr = (dst_x - x)*(dst_x - x) + (dst_y - y)*(dst_y - y)
            if( ({rsqr} >= {rminsqr}) && ({rsqr} <= {rmaxsqr}) )
                prob = {prob0}*{exp {-1.0*{rsqr*decay}}}
                if ({rand 0 1} <= prob)
                  // dst_node = {trunc {j/cells_per_slice}} + {worker0}
                  dst_node =  {cell_node {j} {cells_per_slice}}
                  // dst_index = j  - (dst_node - worker0)*cells_per_slice
                  dst_index = {cell_index {j} {cells_per_slice}}
                  dst = {dst_cell_path} @ "[" @ {dst_index} @ "]/" @ {dst_synpath}
                  raddmsg {src} {dst}@{dst_node} SPIKE
                end
            end
        end // dst loop
    end // src loop
end // function planarconnect_probdecay



/* Functions to Connect each source spike generator in src_layer to
   target synchans in dst_layer within the specified range. The
   connection probability is weighted with a probability
   exp(-d*d/(sigma*sigma)) where r is the radial distance between
   source and destination in the x-y plane.

   With the flag use_probdecay_function = 1, the function planarconnect_probdecay
   will be used to make the connections between cells.

   When use_probdecay_function = 0, the built-in command volumeconnect3
   is used with the '-planar' option. This provides faster setup, but a
   rvolumeconnect3 command is not available in PGENESIS.


  Default values used in this simulation:
  float Ex_Ex_prob = 0.15
  float Ex_Inh_prob = 0.45
  float Inh_Ex_prob = 0.6
  float Inh_Inh_prob = 0.6  // just a guess

 Typical usage:

  connect_Ex_Inh /Ex_layer4 /Inh_layer4
  connect_Ex_Inh "Ex_"@{default_layer} "Inh_"@{default_layer}  

*/


function connect_Ex_Ex(src_layer,dst_layer)
  str src_layer, dst_layer
  float Ex_Ex_prob = 0.15
  float Ex_Inh_prob = 0.45
  float Inh_Ex_prob = 0.6
  float Inh_Inh_prob = 0.6  // just a guess
  // distance at which number of targets cells becomes very small
  float pyr_range = 25*Ex_SEP_X 
  float bask_range = 25*Ex_SEP_X

  // Ex_layer cells connect to excitatory synchans
  float Ex_sigma = 10.0*Ex_SEP_X // values from K Yaun SfN 2008 poster fit
  float Inh_sigma = 10.0*Ex_SEP_X

  planarconnect_probdecay /{src_layer}/{Ex_cell_name} "soma/spike" \
      {Ex_cells_per_slice} \
      /{dst_layer}/{Ex_cell_name}  {Ex_ex_synpath} {Ex_NX*Ex_NY} \
      0.0 {pyr_range} {Ex_Ex_prob} {Ex_sigma}
end


function connect_Ex_Inh(src_layer,dst_layer)
  str src_layer,dst_layer
  float Ex_Ex_prob = 0.15
  float Ex_Inh_prob = 0.45
  float Inh_Ex_prob = 0.6
  float Inh_Inh_prob = 0.6  // just a guess
  // distance at which number of targets cells becomes very small
  float pyr_range = 25*Ex_SEP_X 
  float bask_range = 25*Ex_SEP_X

  // Ex_layer cells connect to excitatory synchans
  float Ex_sigma = 10.0*Ex_SEP_X // values from K Yaun SfN 2008 poster fit
  float Inh_sigma = 10.0*Ex_SEP_X

  planarconnect_probdecay /{src_layer}/{Ex_cell_name} "soma/spike" \
      {Ex_cells_per_slice} \
      /{dst_layer}/{Inh_cell_name}  {Inh_ex_synpath} {Inh_NX*Inh_NY} \
      0.0 {pyr_range} {Ex_Inh_prob} {Ex_sigma}
end

function connect_Inh_Ex(src_layer,dst_layer)
  str src_layer,dst_layer
  float Ex_Ex_prob = 0.15
  float Ex_Inh_prob = 0.45
  float Inh_Ex_prob = 0.6
  float Inh_Inh_prob = 0.6  // just a guess
  // distance at which number of targets cells becomes very small
  float pyr_range = 25*Ex_SEP_X 
  float bask_range = 25*Ex_SEP_X

  // Ex_layer cells connect to excitatory synchans
  float Ex_sigma = 10.0*Ex_SEP_X // values from K Yaun SfN 2008 poster fit
  float Inh_sigma = 10.0*Ex_SEP_X

  planarconnect_probdecay /{src_layer}/{Inh_cell_name} "soma/spike" \
      {Inh_cells_per_slice} \
      /{dst_layer}/{Ex_cell_name}  {Ex_inh_synpath} {Ex_NX*Ex_NY} \
      0.0 {bask_range} {Inh_Ex_prob} {Inh_sigma}
end

function connect_Inh_Inh(src_layer,dst_layer)
  str src_layer,dst_layer
  float Ex_Ex_prob = 0.15
  float Ex_Inh_prob = 0.45
  float Inh_Ex_prob = 0.6
  float Inh_Inh_prob = 0.6  // just a guess
  // distance at which number of targets cells becomes very small
  float pyr_range = 25*Ex_SEP_X 
  float bask_range = 25*Ex_SEP_X

  // Ex_layer cells connect to excitatory synchans
  float Ex_sigma = 10.0*Ex_SEP_X // values from K Yaun SfN 2008 poster fit
  float Inh_sigma = 10.0*Ex_SEP_X

  planarconnect_probdecay /{src_layer}/{Inh_cell_name} "soma/spike" \
      {Inh_cells_per_slice} \
      /{dst_layer}/{Inh_cell_name}  {Inh_inh_synpath} {Inh_NX*Inh_NY} \
      0.0 {bask_range} {Inh_Inh_prob} {Inh_sigma}
end

/* function connect_cells is performed by each worker node to make
   connections from its "slice" of the network
*/
function connect_cells
  /* Default values used in this simulation */
  float Ex_Ex_prob = 0.15
  float Ex_Inh_prob = 0.45
  float Inh_Ex_prob = 0.6
  float Inh_Inh_prob = 0.6  // just a guess

  // distance at which number of targets cells becomes very small
  float pyr_range = 25*Ex_SEP_X 
  float bask_range = 25*Ex_SEP_X

  // Ex_layer cells connect to excitatory synchans

    float Ex_sigma = 10.0*Ex_SEP_X // values from K Yaun SfN 2008 poster fit
    float Inh_sigma = 10.0*Ex_SEP_X

    if(use_probdecay_function)
        str src_layer = "Ex_" @ {default_layer}
        str dst_layer = "Ex_" @ {default_layer}
        connect_Ex_Ex {src_layer} {dst_layer}

        str dst_layer = "Inh_" @ {default_layer}
        connect_Ex_Inh  {src_layer} {dst_layer}

        str src_layer = "Inh_" @ {default_layer}
        str dst_layer = "Ex_" @ {default_layer}
        connect_Inh_Ex  {src_layer} {dst_layer}

        str dst_layer = "Inh_" @ {default_layer}
        connect_Inh_Inh  {src_layer} {dst_layer}
    else // only for serial versions with use_probdecay_function = 0
      str src_layer = "Ex_" @ {default_layer}
      str dst_layer = "Ex_" @ {default_layer}
      volumeconnect3 /{src_layer}/{Ex_cell_name}[]/soma/spike \
       /{dst_layer}/{Ex_cell_name}[]/{Ex_ex_synpath} \
		-relative \  // Destination coordinates are measured relative to source
       -sourcemask box -1 -1 -1 1 1 1 \   // Larger than source area ==> all cells
       -destmask ellipsoid 0 0 0 {pyr_range} {pyr_range} {Ex_SEP_Z*0.5}  \
       -desthole box {-Ex_SEP_X*0.5} {-Ex_SEP_Y*0.5} {-Ex_SEP_Z*0.5} \
         {Ex_SEP_X*0.5} {Ex_SEP_Y*0.5} {Ex_SEP_Z*0.5} \
        -probability  {Ex_Ex_prob} -decay {1.0/Ex_sigma} {Ex_Ex_prob} 0.0 2.0 \
       -planar // restrict distance to x-y plane
    str dst_layer = "Inh_" @ {default_layer}
    // Ex-Inh connections don't need an intial desthole
    volumeconnect3 /{src_layer}/{Ex_cell_name}[]/soma/spike \
    /{dst_layer}/{Inh_cell_name}[]/{Inh_ex_synpath} \
    -relative \	    // Destination coordinates are measured relative to source
    -sourcemask box -1 -1 -1  1 1  1 \  // Larger than source area ==> all cells
    -destmask ellipsoid 0 0 0 {pyr_range} {pyr_range}  {Ex_SEP_Z*0.5} \
   -probability {Ex_Inh_prob} -decay {1.0/Ex_sigma} {Ex_Inh_prob} 0.0 2 \
   -planar
   str src_layer = "Inh_" @ {default_layer}
   str dst_layer = "Ex_" @ {default_layer}
    // Inh-Ex connections don't need an intial desthole
    volumeconnect3 /{src_layer}/{Inh_cell_name}[]/soma/spike \
    /{dst_layer}/{Ex_cell_name}[]/{Ex_inh_synpath} \
    -relative \	    // Destination coordinates are measured relative to source
    -sourcemask box -1 -1 -1  1 1  1 \  // Larger than source area ==> all cells
    -destmask ellipsoid 0 0 0  {bask_range} {bask_range}  {Inh_SEP_Z*0.5} \
    -probability {Inh_Ex_prob}  -decay {1.0/Inh_sigma} {Inh_Ex_prob} 0.0 2 \
    -planar // restrict distance to x-y plane

    // This could be optional
    // Inh_layer4 cells connect to inhibitory synchans
    str src_layer = "Inh_" @ {default_layer}
    str dst_layer = "Inh_" @ {default_layer}
    volumeconnect3 {src_layer}/{Inh_cell_name}[]/soma/spike \
    /{dst_layer}/{Inh_cell_name}[]/{Inh_inh_synpath} \
    -relative \	    // Destination coordinates are measured relative to source
    -sourcemask box -1 -1 -1  1 1  1 \   // Larger than source area ==> all cells
    -destmask ellipsoid 0 0 0 {bask_range} {bask_range}  {Inh_SEP_Z*0.5}  \
    -desthole box {-Inh_SEP_X*0.5} {-Inh_SEP_Y*0.5} {-Inh_SEP_Z*0.5} \
       {Inh_SEP_X*0.5} {Inh_SEP_Y*0.5} {Inh_SEP_Z*0.5} \
    -probability {Inh_Inh_prob} -decay {1.0/Inh_sigma} {Inh_Inh_prob} 0.0 2 \
    -planar // restrict distance to x-y plane           
  end // if
end // function connect_cells

// Utility functions to calculate statistics
function print_avg_syn_number // to be executed by control_node
    int n, dst_node, dst_index
    float num_Ex_ex = 0
    float num_Ex_inh = 0
    str Ex_layer = "Ex_" @ {default_layer}
    str Inh_layer = "Inh_" @ {default_layer}
    for (n=0; n < Ex_NX*Ex_NY; n=n+1)
        dst_node = {cell_node {n} {Ex_cells_per_slice}}
        dst_index = {cell_index {n} {Ex_cells_per_slice}}
        num_Ex_ex = num_Ex_ex + {getsyncount@{dst_node} \
            /{Ex_layer}/{Ex_cell_name}[{dst_index}]/{Ex_ex_synpath}}
        num_Ex_inh = num_Ex_inh + {getsyncount@{dst_node} \
            /{Ex_layer}/{Ex_cell_name}[{dst_index}]/{Ex_inh_synpath}}
    end
    num_Ex_ex = num_Ex_ex/(Ex_NX*Ex_NY)
    num_Ex_inh = num_Ex_inh/(Ex_NX*Ex_NY)

    float num_Inh_ex = 0
    float num_Inh_inh = 0
    for (n=0; n < Inh_NX*Inh_NY; n=n+1)
        dst_node = {cell_node {n} {Inh_cells_per_slice}}
        dst_index = {cell_index {n} {Inh_cells_per_slice}}
        num_Inh_ex = num_Inh_ex + {getsyncount@{dst_node} \
            /{Inh_layer}/{Inh_cell_name}[{dst_index}]/{Inh_ex_synpath}}
        num_Inh_inh = num_Inh_inh + {getsyncount@{dst_node} \
            /{Inh_layer}/{Inh_cell_name}[{dst_index}]/{Inh_inh_synpath}}
    end
    num_Inh_ex = num_Inh_ex/(Inh_NX*Inh_NY)
    num_Inh_inh = num_Inh_inh/(Inh_NX*Inh_NY)

    echo "Average number of Ex_ex synapses per cell: " {num_Ex_ex}
    echo "Average number of Ex_inh synapses per cell: " {num_Ex_inh}
    echo "Average number of Inh_ex synapses per cell: " {num_Inh_ex}
    echo "Average number of Inh_inh synapses per cell: " {num_Inh_inh}

end // print_avg_syn_number

function print_avg_input_number
// to be executed by control_node
   int n, dst_node, dst_index, row, min, max
    float num_Ex_drive
    str Ex_layer = "Ex_" @ {default_layer}
    for (row=0; row < Ex_NY; row=row+1)
        num_Ex_drive = 0.0
        min = row*Ex_NX; max = min + Ex_NX
        for (n=min; n<max; n=n+1)
            dst_node = {cell_node {n} {Ex_cells_per_slice}}
            dst_index = {cell_index {n} {Ex_cells_per_slice}}
            num_Ex_drive = num_Ex_drive + {getsyncount@{dst_node} \
	    /{Ex_layer}/{Ex_cell_name}[{dst_index}]/{Ex_drive_synpath}}
	end
        // num_Ex_drive = 1.0*num_Ex_drive
        echo "Row "{row} "Average number of Ex_drive synapses per row: " {num_Ex_drive}
    end // loop over rows
end // function print_avg_input_number
