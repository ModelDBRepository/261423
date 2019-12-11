#!/usr/bin/env python

# LFPtrials3 - a utility to plot a local field potential and optional
# Power Spectral Density from time series data that is organized into
# 'trials' with a specified 'stimulus type' applied during each trial. 

import sys, os, math
import numpy as np


def print_usage():
    usage_str = """
Usage: LFPtrials-out.py trial_param_file time_data_file ["optional comment string"]
e.g.: LFPtrials-out.py oddball_freqs_6009.txt Electrode_3_6009.txt

LFPtrials3 and PSDtrials2 are utilities to plot a trial averaged Local
Field Potential and Power Spectral Density from time series data that
is organized into 'trials' with an identical stimulus applied during
each trial.

LFPtrials-out.py outputs the LFP to be plotted for stimulus type 1
(the frequent stimulus) to a file 'LFP-out1_<time_data_file>'.
Likewise, it performes the calculation and outputs a file
'LFP-out2_<time_data_file>' if samples exist for stimulus 2.

for stimulus type 2.
These versions do not plot the LFP.

The 'trial_param_file' contains a header with information about the
stimulus parameters. The header should be in the format:

    #RUNID  stim_delay stim_length trial_interval sample_length

This is followed by lines, with each one containing the trial start time and
the stimulus type:

    0 = none, 1 = preferred stimulus, 2 = oddball (infrequent) stimulus

NOTE: Calculations are performed for all trials up to, but not including
the final entry in this file.

The 'time_data_file' typically contains some measure of Excitatory
Post-Synaptic Currents in the network over a run that has multiple
trials.  This typically contains Local Field Potentials recorded at a
distance from a network of cells. The average LFP(t) or EPSC(t) is
calculated by averaging over trials.  It is plotted at the top.

The PSDtrials version calculates the PSD from the trial averaged EPSC(t).
The PSD is shown in the lower plot.

The LFPtrials2 and PSDtrials3 versions calculate a PSD from the time
series for each trial, first subtracting the trial averaged EPSC(t) at
each time point. The average PSD is shown in the lower plot. This
isolates the components of the response that are not phase-locked to
the stimulus.  """
    print(usage_str)
    sys.exit


def plot_dataset(trial_param_file,time_data_file,comment_str,pref_stim):
    print(trial_param_file, time_data_file,pref_stim)

    # Get information from the trial_data_file header
    pfp = open(trial_param_file, 'r')
    lines = pfp.readlines()
    
    # check if first line is a header line starting with '#'
    # This should be in the format:
    #RUNID  stim_delay stim_length trial_interval sample_length
    header = lines[0].split()
    if header[0][0] == "#": 
        RUNID = header[0][1:]
        stim_delay = float(header[1])
        stim_length  = float(header[2])
        trial_interval = float(header[3])
        sample_length = float(header[4])
    else:
        print("No header information found!")
        sys.exit()
    
    # print RUNID, stim_delay, stim_length, trial_interval, sample_length

    ntrials = len(lines) - 1
    # count only the trials with this stimulus type
        # allocate array space for the stimulus parameter data
    trial_time =  np.zeros(ntrials)
    stim_type =   np.zeros(ntrials)
    
    # Then fill arrays from remaining lines
    first_pref_trial = -1
    for i in range(1, ntrials+1):
        data = lines[i].split()
        j = i -1 # trial number 0 is line 1
        trial_time[j] = float(data[0])
        stim_type[j] = int(data[1])
        if (first_pref_trial < 0) and (stim_type[j] == pref_stim):
            first_pref_trial = j
        # print trial_time[j], stim_type[j]

    pfp.close ()

    file_out = True

    print('Reading %s' % time_data_file)
    fp = open(time_data_file, 'r')
    lines = fp.readlines()
    count = len(lines)
    print(count, " lines")

    # determine size of the time step
    line_0 = lines[0].split(" ")
    line_1 = lines[1].split(" ")
    dt = float(line_1[0]) - float(line_0[0])
    print("time step = ", dt)

    # determine size of the sample data
    sample_items =  int(sample_length/dt + 0.5) + 1

    # Before the first trial, allocate space
    tn = np.zeros(sample_items)
    yn = np.zeros(sample_items)
    # print len(tn), sample_items

    # Now read the data
    i=0
    trial_num = 0
    nsamples = 0
    line_num = 0
    finished = False
    # last trial to count is the one before final enty in trial_param_file
    last_trial = ntrials - 2
    trial_end = trial_time[trial_num + 1]
    stype = stim_type[trial_num]
    if stype == pref_stim:
        nsamples +=1
    for line in lines:
        data = line.split(" ")
        t = dt*line_num
        line_num += 1
        # see if we have passed the end of the trial
        if t > trial_end:
            if trial_num < last_trial:
                # start over at the begining of the next sample
                trial_num += 1
                i = 0
                trial_end = trial_time[trial_num + 1]
                stype = stim_type[trial_num]
                if stype == pref_stim:
                     nsamples +=1
                
            else: # skip (likely incomplete) data for following trial
                finished = True

        sample_t = math.fmod(t, trial_interval)
        
        # Add to the array if the sample time is within the sample_length
        if not finished and (i < sample_items) and (stype == pref_stim):
            # Note that tn[i] is replaced, and yn[i] is added to
            tn[i] = sample_t; yn[i] = yn[i] + float(data[1])
            i += 1

    # end of loop over lines
    
    if (nsamples < 1):
        return
    
    t_final = trial_end*trial_interval # run duration to print

    # Average the data over number of samples having 'pref_stim' stim type
    avg_LFP = yn/nsamples

    # Option to subtract the base level LFP (averaged over t < stim_delay)
    # from LFP(t)

    normalize_baseline = True
    n_base_samples = 0
    LFPbase = 0.0

    imax = int(stim_delay/dt)

    if (normalize_baseline):
        for i in range (0,imax):
            LFPbase += yn[i]
            n_base_samples += 1

        LFPbase /= (n_base_samples*nsamples)
        # print ("baseline average ", LFPbase, " over pref_stim ", pref_stim)
        avg_LFP -= LFPbase

    if (file_out):
        outfile_name = 'LFP-out' + str(pref_stim) + '_' + time_data_file
        print('Writing ',  outfile_name)
        outfp = open(outfile_name, 'w')
        imax = int(sample_length/dt)
        for i in range (0,imax):
            outfp.write(str(tn[i]) +  " " + str(avg_LFP[i]) + "\n")

    # now do the calculation
    print('alculating average of ', nsamples, ' preferred stimulus samples from ', RUNID)
    print("number of trials = ", ntrials-1, " preferred stimulus type = ", pref_stim)
    print('run duration (sec): ', t_final, ' data aquisition interval: ', dt)

    # print run description to the "plot" ax3
    text_str = 'Run ' + runid + ' with trial and stimulus parameters from: ' + \
    trial_param_file +'\n' + 'summed Excitatory Local Field Potentials from: ' \
    + time_data_file + '\n' + str(nsamples) +  ' preferred stimulus samples from ' + \
    str(ntrials-1) + ' trials with preferred stiumulus type ' + str(pref_stim) + '\n' + \
    'sample length (sec): ' +  str(sample_length) + '; stimulus delay: ' + str(stim_delay) + \
     '; pulse width: ' + str(stim_length) + '\n' + \
     'trial interval: ' + str(trial_interval) + '; data interval: ' + str(dt) + \
     '; run time ' + str(t_final)  + '\n\n' + 'Comment: ' + comment_str

    if (pref_stim == 1):
        print(text_str)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print('Need a trial parameter file and a EPSC data file')
        print_usage()
        sys.exit()
    # Get the filename arguments
    trial_param_file = sys.argv[1]
    time_data_file = sys.argv[2]
    comment_str = ''
    if len(sys.argv) > 3:
        comment_str = sys.argv[3]
    # Generate a RUNID from a string like "EPSC_sum_M0004E.txt"
    fn1 = time_data_file
    fnbase,ext = os.path.splitext(fn1)
    # get string following final '_' and remove 1 char suffix
    runid = fnbase.split('_')[-1][:]

    print(trial_param_file, time_data_file, runid)


    plot_dataset(trial_param_file, time_data_file, comment_str, 1)
    plot_dataset(trial_param_file, time_data_file, comment_str, 2)

