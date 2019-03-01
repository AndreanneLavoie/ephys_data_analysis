function data_analysis_ori(spikes, assigns)

%compute firing rate for orientation tuning:

%receptive field parameters;
num_angles = length(unique(spikes.vs_params(:, 10)));

%time window parameters
twin_start = 0;
twin_stop = 4; %in sec
twin_baseline = -0.9;

%generate recpetive field index matrix (of size row * col); %on = white;
%black = off
stim_cond_list = linspace(0, 360, num_angles + 1);
theta = stim_cond_list/180*pi;

%place holder
firing_rates = zeros(size(stim_cond_list));

firing_rates_nosubstr=zeros(size(stim_cond_list));

firing_rates_base= zeros(size(stim_cond_list));
%first do white only; use filter function to extract all the stim cond from
%
for i = 1:num_angles
       
    filtered_spikes = filtspikes(spikes, 0, 'stimcond', stim_cond_list(i), 'assigns', assigns);

    %compute firing rate for each stim cond and substract baseline
    %firing rate
    temp_fr = computeFR(filtered_spikes, [twin_start twin_stop]) - computeFR(filtered_spikes, [twin_baseline twin_start]);
    firing_rates(i) = temp_fr;
    firing_rates_nosubstr(i)=computeFR(filtered_spikes, [twin_start twin_stop]);
    firing_rates_base(i)=computeFR(filtered_spikes, [twin_baseline twin_start]);
end

firing_rates(end) = firing_rates(1);
firing_rates_nosubstr(end)=firing_rates_nosubstr(1);
firing_rates_base(end)=firing_rates_base(1);
%plot
figure;
polar(theta, firing_rates);
figure;
polar(theta, firing_rates_nosubstr);
figure;
polar(theta, firing_rates_base);
