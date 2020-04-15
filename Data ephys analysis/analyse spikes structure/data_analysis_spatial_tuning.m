function sf_tuning_fr = data_analysis_spatial_tuning(spikes, assigns, plot)
%INPUT
%spike construct 
%assigns single vector containing any number of assigns
%plot 1 = true' 0 = false

%OUTPUT
%plot tuning curve
%fr for each stim cond

%compute firing rate for orientation tuning:

%receptive field parameters;
stim_cond_list = unique(spikes.vs_params(:, 1));
num_stim = length(stim_cond_list);

%time window parameters
pretrig = spikes.vs_params(1, 9)/1000;
time_bef_stim = spikes.vs_params(1, 6)/1000;
duration =  spikes.vs_params(1, 7)/1000;
time_aft_stim =  spikes.vs_params(1,8)/1000;
vs_start = pretrig; 
vs_stop = vs_start + duration; %in sec


%place holder
sf_tuning_fr = zeros(size(stim_cond_list));
%firing_rates_nosubstr=zeros(size(stim_cond_list));
%firing_rates_base= zeros(size(stim_cond_list));


for i = 1:num_stim
    
    if isempty(assigns)
        filtered_spikes = filtspikes(spikes, 0, 'stimcond', stim_cond_list(i));
        filtered_spikes.sweeps.trials = 1: size(spikes.vs_params, 1); %field required for computeFR function
    else
        filtered_spikes = filtspikes(spikes, 0, 'stimcond', stim_cond_list(i), 'assigns', assigns);
        filtered_spikes.sweeps.trials = 1: size(spikes.vs_params, 1); %field required for computeFR function
    end
    %compute firing rate for each stim cond and substract baseline
    %firing rate
    temp_fr = computeFR(filtered_spikes, [vs_start vs_stop]) - computeFR(filtered_spikes, [0 vs_start]);
    sf_tuning_fr(i) = temp_fr;
    %firing_rates_nosubstr(i)=computeFR(filtered_spikes, [vs_start vs_stop]);
    %firing_rates_base(i)=computeFR(filtered_spikes, [0 vs_start]);
end


%plot
h1 = figure;

%subplot(3,1,1);
plot(stim_cond_list, sf_tuning_fr);
%ylim([-1 2]);

%set(h1, 'position', [600,1,500, 500]);
