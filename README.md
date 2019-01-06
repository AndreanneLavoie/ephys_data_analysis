# ephys_data_analysis_prep
Converts data generated by ephys recording, potogenetics, visual stimuli order etc. into organised data structure.

INSTALLATION:

Dependencies:
- matlab
- NPY-matlab

TO RUN:
1. Save all data into a single folder; Incluing:

from intan: (record of analog inputs for different stim; led, eye track, vs timing, behaviour etc.)

   - 'time.dat';
   - 'board-ADC-00.dat' (00-04, as required);
   - 'info.rhd'; 
   
from psychopy: 
        
   - '..._vs_date.csv' (order of stimuli)

from kilosort/phy: 

   - 'spike_times.npy';
   - 'spike_clusters.npy';
   - 'cluster_groups.tsv;

2. Run create_spikes.m
      Note: you will be asked to select the info.rhd file (which should be in your folder containing everything), therefore
            your path will be automatically saved;
            
      BUG: must modify the vs_stim_order file name, as there is a current bug that prevents it from automatically being found
      
