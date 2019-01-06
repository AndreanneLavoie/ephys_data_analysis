%abs_spiketimes = [1.1 2.2 3.3 4.4 5.5 6.6 7.7 8.8 9.9];
%trigger_times = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
%vs_stim_order = [30 60 90 120 150 180 210 240 270];

function [rel_spiketimes, trials, stim_cond] = make_rel_spiketimes(abs_spiketimes,trigger_times, vs_stim_order)
%Using matrix computation, outputs the relative spiketimes (to vs trigger
%time) of each spike, given absolute spike times and trigger times.

%Also returns the trial ID of each spikes and the type of stim condition
%based on relative spiketimes calculated.

%using matrix algebra to speed up computation 
Abs_mat = repmat(abs_spiketimes, length(trigger_times), 1);
Trig_mat = repmat(trigger_times', 1, length(abs_spiketimes));

Dif = Trig_mat - Abs_mat; %calculate the time difference between nearest trigger and spike 

Temp = Dif < 0; %gen matrix where 1 negative num and 0 positive; interface (1 0) is the vs identity 
Temp_der = diff(Temp); %take derivitive to find values that change from 1 to 0;

rel_spiketimes = Dif(find(Temp_der == -1)); %return the relative spiketimes corresponding to each spike

[r c] = find(Temp_der == -1);  %extract coordinates where -1 located in matrix 
triggerID = 1:length(trigger_times); %ID number each trigger_time starting at 1
trials = triggerID(r); %identify the vs stim trial num for each spike

stim_cond = vs_stim_order(trials); %return the type of stimuli assigned to each spike based on newly assigned trial ID
end

%rel_spiketimes
%trials
%stim_cond