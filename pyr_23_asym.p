// Cell parameter file for the reduced model of "cell 945" at 945 um below
// the surface, from Eyal et al. 2016 Fig 1c4.
// RM and Cm values are scaled by a ratio of the model segment area to
// that of same region of the reconstructed cell.
// These were obtained by parsing 2013_03_13_cell06_945_H42_05.CNG.swc
// to identify all terminal and non terminal branch intervals in regions
// at specified distances from the some

*relative
*cartesian
*asymmetric

// membrane constants (SI units) - scaled to give area corrections
*set_global        RM      3.10     // ohm*m^2
*set_global        RA      2.47      // ohm*m
*set_global        CM      0.0061    // farad/m^2
*set_global     EREST_ACT  -0.082     // volts (leakage potential)

// These include an additional scaling of RM and RA by 0.85 and CM by 1.18
// This gives another 15% increase in area, preserving  RM/RA ratio,
// resulting in better fit to 2 msec, 200pA injection pulse.

// Populate soma wih modified traub91 channels

soma    none  0  0 19 18.5   \
    Na_pyr             3000  \
    Kdr_pyr             1600  \
    Ca_hip_traub91      100  \
    Kahp_pyr             30  \
    Ca_conc            -7.769e12    \
    spike 0.0

*set_compt_param RA 1.24 // use average soma-apical0 Ra
apical0 soma       0   0  66    2.9  GABA_pyr 5.0 // apical trunk
*set_compt_param RA 2.47 // back to default
apical1 apical0    0   0  335   1.8  AMPA_pyr 5.0 // lower apical at 50 - 300 um

// For upper apical divide Rm and multiply CM by 2.08
// Apply additional spine area correction of 1.9
*set_compt_param RM 0.923
*set_compt_param CM 0.0206
apical2 apical1    0   0  267   2.3 GABA_pyr 5.0 // mid apical at 300 - 600 um
apical3 apical2    0   0  84    2.3 // apical tuft non-terminal > 600 um
apical4a apical3    -189   0  189   2.3 // terminal tuft btanches
apical4b apical3     189   0  189   2.3

// For oblique, scale by 3.7
// Apply additional spine area correction of 1.9*0.85
*set_compt_param RM 0.441
*set_compt_param CM 0.0432

oblique1 apical0 -34   0   0    4.9 //  AMPA_pyr 2.6526
oblique2a oblique1 -71   0  -71 2.3   GABA_pyr 2.6526
oblique2b oblique1 -71   0   71 2.3   AMPA_pyr 2.6526

// For basal0, scale factor is 2.78 - no spine correction, 0.85,1.18
*set_compt_param RM 1.116
*set_compt_param CM 0.0171
*set_compt_param RA 1.31 // use average soma-basal0 Ra
basal0  soma     0     0  -32   6.0  // GABA_pyr 2.6526
*set_compt_param RA 2.47// back to default

// Apply additional spine area correction of 1.9
// For basal1, scale factor is 5.0 (*0.85)
*set_compt_param RM 0.327
*set_compt_param CM 0.0583
basal1a  basal0 -22.0   0  -22.0 3.5
basal1b basal0   22.0   0 -22.0 3.5
// For basal2, scale factor is 3.22*1.9 (*0.85)
*set_compt_param RM 0.507
*set_compt_param CM 0.0375
basal2a basal1a -200    0  0    3.5
basal2b basal1a  0    0  -200   3.5   AMPA_pyr 5.0
basal2c basal1b  0    0  -200   3.5
basal2d basal1b  200    0  0   3.5
