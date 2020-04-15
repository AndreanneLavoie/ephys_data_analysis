function [stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, assigns, tuning, plot_flag, led)
%this function will compute the tuning curve of given units and can plot it

if nargin<4 %plot 1 = true' 0 = false
    plot_flag=1;
end
if nargin<5 % only led off
    led=0;
end
%INPUT
%spike construct 
%assigns single vector containing any number of assigns
%tuning: 'sf', 'ori', 'tf'


%OUTPUT
%plot tuning curve
%fr for each stim cond

%compute firing rate for tunig parameter

%assing tuning parameter;
switch tuning
    case 'sf'
        tuning_index = 1;
    
    case 'tf'
        tuning_index = 2;
        
    case 'ori'
        tuning_index = 10;
end
%receptive field parameters;
stim_cond_list = unique(spikes.vs_params(:, tuning_index));
num_stim = length(stim_cond_list);

%time window parameters
pretrig = spikes.vs_params(1, 9)/1000;
time_bef_stim = spikes.vs_params(1, 6)/1000;
duration =  spikes.vs_params(1, 7)/1000;
time_aft_stim =  spikes.vs_params(1,8)/1000;
vs_start = pretrig; 
vs_stop = vs_start + duration; %in sec


%place holder
tuning_fr_total = zeros(size(stim_cond_list));
tuning_fr_spon = zeros(size(stim_cond_list));
%firing_rates_nosubstr=zeros(size(stim_cond_list));
%firing_rates_base= zeros(size(stim_cond_list));


for i = 1:num_stim
    
    if isempty(assigns)
        filtered_spikes = filtspikes(spikes, 0, 'stimcond', single(stim_cond_list(i)),'led', led);
        filtered_spikes.sweeps.trials = 1: size(spikes.vs_params, 1); %field required for computeFR function
    else
        filtered_spikes = filtspikes(spikes, 0, 'stimcond', stim_cond_list(i), 'assigns', assigns,'led',led);
        filtered_spikes.sweeps.trials = 1: size(spikes.vs_params, 1); %field required for computeFR function
    end
    %compute firing rate for each stim cond and substract baseline
    %firing rate
    tuning_fr_total(i) = computeFR(filtered_spikes, [vs_start vs_stop]);
    tuning_fr_spon(i) = computeFR(filtered_spikes, [0 vs_start]);
    %firing_rates_nosubstr(i)=computeFR(filtered_spikes, [vs_start vs_stop]);
    %firing_rates_base(i)=computeFR(filtered_spikes, [0 vs_start]);
end


%plot
if plot_flag
    tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;
    %figure;
    switch tuning
    case 'sf'
        plot(stim_cond_list, tuning_fr_evoked);

    case 'tf'
        plot(stim_cond_list, tuning_fr_evoked);

    case 'ori'
        if any(tuning_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
            tuning_fr_evoked = tuning_fr_evoked-min(tuning_fr_evoked);
            fprintf('tuning_fr_evoked contains negative value, force minimum to be zero\n');
        end
        theta = stim_cond_list/180*pi;
        theta = [theta; theta(1)];
        tuning_fr_evoked = [tuning_fr_evoked; tuning_fr_evoked(1)];
        %tuning_fr_total = [tuning_fr_total; tuning_fr_total(1)];
        polar(theta, tuning_fr_evoked);
        %polar(theta, tuning_fr_total);
    end
end
