
function spikes = make_spike_struct(spikes)
%Using matrix computation, outputs the relative spiketimes (to vs trigger
%time) of each spike, given absolute spike times and trigger times.

%Also returns the trial ID of each spikes and the type of stim condition
%based on relative spiketimes calculated.

%error check
size_vs_para = size(spikes.vs_params, 1);
size_vstiming = size(spikes.vstiming, 1);
size_trigtiming = size(spikes.trigtiming, 1);

if ~ ((size_vs_para== size_vstiming) && ( size_vs_para == size_trigtiming ))
    disp('error: size spikes.vs_params, spikes.vstiming and/or spikes.trigtiming do not match');
    
    min_size = min([size_vs_para size_vstiming size_trigtiming]);
    spikes.vs_params = spikes.vs_params(1:min_size, :);
    spikes.vstiming = spikes.vstiming(1:min_size, :);
    spikes.trigtiming = spikes.trigtiming(1:min_size, :);
end  

%First, to remove spikes that happen before visual stimulation start and
%after it ends,arbitrarily assign 0 as the number of trigger for spikes before the first
% trigger and after the last trial
expttime=[spikes.trigtiming(1,1) spikes.trigtiming(end,1)+sum(spikes.vs_params(1, 7:9))/1000.0]; %range of time that visual stimulation is presented
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
spikes.trials = spikes.trials';
spikes.trigger = spikes.trials; %required for filtspikes function

%return the type of stimuli assigned to each spike based on newly assigned trial ID
spikes.stimcond = nan(size(spikes.trials));
spikes.stimcond(spikes.trials > 0) = spikes.vs_params(spikes.trials(spikes.trials > 0), 10);  %column 10 in stim_cond matrix is the column that contains the tuning parameter that is changing 

%calculate relative spiketimes (realtive to vstiming)
spikes.spiketimes = spikes.abs_spiketimes; 
spikes.spiketimes(spikes.trials > 0) = spikes.abs_spiketimes(spikes.trials > 0) - spikes.vstiming(spikes.trials(spikes.trials > 0),1);
spikes.spiketimes(spikes.trials <= 0)= nan;
spikes.spiketimes = spikes.spiketimes';


%make sweeps field for filter function
%sweeps info for each trials
spikes.sweeps.trials = 1: size(spikes.vs_params, 1);
spikes.sweeps.fileInd = ones(size(spikes.sweeps.trials));
spikes.sweeps.trigger = spikes.sweeps.trials;
spikes.sweeps.time = spikes.trigtiming(:,1)';
spikes.sweeps.stimcond = spikes.vs_params(:,10)';
spikes.sweeps.led = zeros(size(spikes.sweeps.trials));  

%make another field to match psth_bl code; total time that psth plot will be over
spikes.info.detect.dur = spikes.sweeps.trials %proper size 
spikes.info.detect.dur(:) = ( spikes.vs_params(1,9) + spikes.vs_params(1,7) + 1000)/1000 %trigger time + t_dur + 1000ms than convert to sec.


%reshape spikes.abs_spiketimes for filtspikes function
spikes.abs_spiketimes = spikes.abs_spiketimes';
spikes.led = zeros(size(spikes.abs_spiketimes));
spikes.time = spikes.abs_spiketimes;


spikes.fileInd = ones(size(spikes.abs_spiketimes)); %if multiple files are generated;


end
