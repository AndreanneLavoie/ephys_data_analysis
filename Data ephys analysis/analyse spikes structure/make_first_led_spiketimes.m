function [length_led, led_spiketimes, led_hist] = make_first_led_spiketimes(spikes, assign, bin, bin_limit, first_spike_only, flag_plot)

%for given cluster (assign #), extract the first spike for each LED trial

%INPUT
%spike structure (requires: ) 
%assign (cluster) that want to extract spiketimes from [1x1]
%bin in ms; should be between 0.1ms and 5ms; 
if nargin<4 %bin limit in ms; e.g. max: [0, 100] or min: [0, 15] for plotting histogram 
    bin_limit=[0, 15];
end

if nargin<5%first_spike_only 1 = true; 0 = faslse
    first_spike_only=1;
end

if nargin<6%plot 1= plot; 0 = not plot
    flag_plot=1;
end

%OUTPUT
%first_led_spiketimes [1 x #led trials]

%%%
%filter spiketimes by assign
if ~isempty(assign)
    filt_spikes = filtspikes(spikes, 0, 'assigns', assign,'led', 1);
else
    filt_spikes = spikes;
end

%create extended matrix to allow comparision between all possible led
%times and all possible abs spiketimes
length_st = length(filt_spikes.abs_spiketimes);
length_led = length(filt_spikes.led_abstime(:,1));

abs_spikestimes_matrix = repmat(filt_spikes.abs_spiketimes, length_led, 1);%sec
abs_ledtimes_matrix = repmat(filt_spikes.led_abstime(:,1), 1, length_st); %sec

%calculate rel spiketime to ref led
diff_matrix = abs_spikestimes_matrix - abs_ledtimes_matrix;

%calculate LED period
led_period_temp = filt_spikes.led_abstime(2:(end),1) - filt_spikes.led_abstime(1:(end-1),1);
led_period_temp(led_period_temp > spikes.vs_params(1,13)/1000) = [];
led_period = mean(led_period_temp);

%remove spikes that are over ledperiod sec  %remove neg value %%ADD as option led period from vsparams to reduce 2.6sec 
diff_matrix(diff_matrix < 0  | diff_matrix >= led_period) = nan; 

if first_spike_only
    %find min of for each led trial (ie the first spike for the cluster)
    led_spiketimes = min(diff_matrix, [], 2);
    %remove all spikes that are larger than led period
    %led_spiketimes = led_spiketimes(led_spiketimes < led_period);
else
    %remove nan
    led_spiketimes = diff_matrix(~ isnan(diff_matrix) );
end 
%convert to ms

%plot
if flag_plot
    led_spiketimes_ms = led_spiketimes*1000;
    led_hist = histogram(led_spiketimes_ms, 'BinWidth', bin, 'BinLimits', bin_limit);
    
    %label axis
    xlabel('Latency (ms)');
    ylabel('# spikes');
else
    led_hist = [];
end

%check size of output var
% num_spiketimes = length(led_spiketimes_ms)
% num_non_nan_spiketimes = length(led_spiketimes_ms(~ isnan(led_spiketimes_ms)))


end