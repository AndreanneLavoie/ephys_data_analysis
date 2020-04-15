function data_analysis_temporal_tuning(spikes, assigns)

%compute firing rate for orientation tuning:

%receptive field parameters;
stim_cond_list = unique(spikes.vs_params(:, 2));
num_stim = length(stim_cond_list);

%time window parameters
pretrig = spikes.vs_params(1, 9)/1000;
time_bef_stim = spikes.vs_params(1, 6)/1000;
duration =  spikes.vs_params(1, 7)/1000;
time_aft_stim =  spikes.vs_params(1,8)/1000;
vs_start = pretrig 
vs_stop = vs_start + duration + time_aft_stim; %in sec


%place holder
fr = zeros(size(stim_cond_list));
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
    fr(i) = temp_fr;
    %firing_rates_nosubstr(i)=computeFR(filtered_spikes, [vs_start vs_stop]);
    %firing_rates_base(i)=computeFR(filtered_spikes, [0 vs_start]);
end

%firing_rates(end) = firing_rates(1);
%firing_rates_nosubstr(end)=firing_rates_nosubstr(1);
%firing_rates_base(end)=firing_rates_base(1);
%plot
h1 = figure;

%subplot(3,1,1);
plot(stim_cond_list, fr);
%subplot(3,1,2);
%plot(stim_cond_list, firing_rates_nosubstr);
%subplot(3,1,3);
%plot(stim_cond_list, firing_rates_base);

set(h1, 'position', [600,1,500, 500]);
