// MBG-protodefs.g - Definition of prototype elements for MGBv relay cell

/* NOTE: This assumes that a previously included protodefs file
   has created /library with the elements below:

include compartments

create spikegen spike
setfield spike  thresh 0.00  abs_refract 1.0e-3  output_amp 1
make_cylind_compartment

*/

// functions to create modified Traub and Miles tabchannels
include MGBchans.g

/* MGBchans.g assigns values to the global variables EREST_ACT, ENA, EK, and
   SOMA_A.  The first three will be superceded by any values defined below.
   The value of SOMA_A is not relevant, as the cell reader calculates the
   compartment area.

   Currently:

   EREST_ACT = -0.063
   ENA       =  0.050
   EK        = -0.090
*/
// Use default values

/* file for synaptic channels */
include synchans // from Neurokit/prototypes

/* synchans.g defines:
   EGlu = 0.045
   EGABA = -0.082
*/

EGlu = 0.0
EGABA = -0.08

// Make a "library element" to hold the prototypes which will be used
// by the cell reader to add compartments and channels to the cell.

if (!{exists /library})     // But, only if it doesn't already exist
    create neutral /library
end

// We don't want the library to try to calculate anything, so we disable it
disable /library

// To ensure that all subsequent elements are made in the library
pushe /library

/* Functions in MGBchans.g are used to create prototype channels
   Na_traub_mod and K_traub_mod   
*/
make_Na_traub_mod
make_K_traub_mod

// Make a prototype excitatory channel, "Ex_channel" - from synchans.g
make_Ex_channel     /* synchan with Ek = 0.0, tau1 = tau2 = 3 msec */

// Make a prototype inhibitory channel, "Inh_channel"
make_Inh_channel     /* synchan with Ek = -0.08, tau1 = tau2 = 20 msec */

pope // Return to the original place in the element tree
