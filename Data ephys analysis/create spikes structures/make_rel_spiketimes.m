
function spikes = make_rel_spiketimes(spikes)
%Using matrix computation, outputs the relative spiketimes (to vs trigger
%time) of each spike, given absolute spike times and trigger times.

%Also returns the trial ID of each spikes and the type of stim condition
%based on relative spiketimes calculated.

%First, to remove spikes that happen before visual stimulation start and
%after it ends,arbitrarily assign 0 as the number of trigger for spikes before the first
% trigger and after the last trial
expttime=[spikes.trigtiming(1,1) spikes.trigtiming(end,1)+sum(spikes.stimorder(1, 7:9))/1000.0]; %range of time that visual stimulation is presented
spikes.trials=ones(size(spikes.abs_spiketimes)); %the trial number for each individual spike; 
spikes.trials(spikes.abs_spiketimes < expttime(1) | spikes.abs_spiketimes > expttime(2))= 0; %0 = no trial assigned to spike

abs_spiketimes_durexp = spikes.abs_spiketimes(spikes.trials > 0); %extract only spikes that happen within exp

%Creating two matrices; #1: spikestimes x triggertimes; #2 triggertimes x spiketimes
%this will allow every possible difference between each spiketime and each
%triggertime to be calculated efficiently, and therefore we will be able to
%assign a trial to each spike.
Abs_mat = repmat(abs_spiketimes_durexp, 1, length(spikes.trigtiming));
Trig_mat = repmat(spikes.trigtiming(:,1)', length(abs_spiketimes_durexp),1); %use the first column of trigger times (c1 = start times; c2 = end times)

%calculate the time difference between nearest trigger and spike
tDif_trig_spike = Trig_mat - Abs_mat; 

%gen matrix were True if negative num (when spike is before trigger) and Flase when positive (when spike is after trigger); 
triggers_bef_spikes = tDif_trig_spike < 0; %determin all triggers before each spikes 
triggers_bef_spikes(:,end+1) = 0; %add last colum = 0 to make sure there is a 1 0 interface at the end;
trigger_ind_bef_spikes = diff(triggers_bef_spikes, [], 2); %find interface (1 change to 0) by taking derivitive to find 
[r c] = find(trigger_ind_bef_spikes < 0);  %extract coordinates where -1 located in matrix 
spikes.trials(spikes.trials > 0) = c; %assign all exp spiketimes a trial ID number corresponding to order of stimuli presented

%return the type of stimuli assigned to each spike based on newly assigned trial ID
spikes.stim_cond = nan(size(spikes.trials));
spikes.stim_cond(spikes.trials > 0) = spikes.stimorder(spikes.trials(spikes.trials > 0), 10);  %column 10 in stim_cond matrix is the column that contains the tuning parameter that is changing 

%calculate relative spiketimes
spikes.spiketimes = spikes.abs_spiketimes;
spikes.spiketimes(spikes.trials > 0) = spikes.abs_spiketimes(spikes.trials > 0) - spikes.vstiming(spikes.trials(spikes.trials > 0),1);
spikes.spiketimes(spikes.trials <= 0)= nan;
end
