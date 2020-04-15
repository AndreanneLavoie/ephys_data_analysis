function [led_flag, real_mean_latency] = shuffle_led_spiketimes_10000(spikes, assign, analysis_window)
%INPUT
%analysis_window in ms

if nargin < 3
    analysis_window = 25; %ms
end

%parameters
plot_flag = 0;
shuffling_number = 10000; 
cut_off_shuffling = 99.9; %percent
%assign = 274;
bin = 1; %ms
bin_limit = [0 100]; %ms
plotting_iter = round(shuffling_number*0.02);
num_bins_window = round(analysis_window/bin);

%calculate led_spiketimes for assign
if plot_flag
    subplot(1,2,1)
    [length_led, led_spiketimes, led_hist] = make_first_led_spiketimes(spikes, assign, bin, bin_limit, 1, plot_flag);
    title('Real data');
    real_fr_bins = led_hist.Values; %already converted to ms in hist
    
else
    [length_led, led_spiketimes] = make_first_led_spiketimes(spikes, assign, bin, bin_limit, 1, plot_flag);
    led_spiketimes = led_spiketimes*1000; %convert to ms
    [real_fr_bins,edges_temp,bin_temp]=histcounts(led_spiketimes, 'BinWidth', bin, 'BinLimits', bin_limit);
end


%spike firing rate in each 1ms bin
real_mean3consBins = movmean(real_fr_bins, 3); %mean of three consecutive 1 ms bins
[max_real_mean3consBins, max_real_index] = max(real_mean3consBins(1:num_bins_window)); %calculate max three bins for real data within the analysis window
real_mean_latency = real_fr_bins(max_real_index)

num_spikes = length(led_spiketimes);

%% 10 000 shuffling
max_shuffled_mean_3consBins = nan(shuffling_number, 1);
tic

%generate random number
%formula will create array of random numbers drawns form a uniform distribution with the open interval [min_range, max_range] (line 20)
for i=1:shuffling_number
    rand_range = [0 100];
    random_num = (rand_range(2)-rand_range(1)).*rand(num_spikes, 1) + rand_range(1);  

    %add random numbers to led_spiketimes
    shuffled_led_spiketimes = led_spiketimes + random_num;

    %plot
    if plot_flag
        subplot(1,2,2)
        %shuffled_led_spiketimes_ms = shuffled_led_spiketimes*1000;
        shuffled_led_hist = histogram(shuffled_led_spiketimes, 'BinWidth', bin, 'BinLimits', bin_limit);
        title('Shuffled data');
        shuffled_fr_bins = shuffled_led_hist.Values; %save the firing rate of each bin
            %label axis
        xlabel('Latency (ms)');
        ylabel('# spikes');
        ylim(led_hist.BinLimits);
    else
        [shuffled_fr_bins, edges_temp, bin_temp]=histcounts(shuffled_led_spiketimes, 'BinWidth', bin, 'BinLimits', bin_limit);
    end
    
   %moving average
    shuffled_mean3consBins = movmean(shuffled_fr_bins, 3); %mean of three consecutive 1 ms bins
    max_shuffled_mean_3consBins(i) = max(shuffled_mean3consBins);
    
    %give current update
    if i==plotting_iter
       [assign i toc]
    end
end
%calculate 99.9th percentile cutoff value
%sort_shuffled_max = sort(max_shuffled_sum_3consBins);
% cutoff_99th_ind = round(shuffling_number*0.999);
% cutoff_99th_val = sort_shuffled_max(cutoff_99th_ind);
cutoff_99th_val = prctile(max_shuffled_mean_3consBins, cut_off_shuffling);

%determine if sum of top 3 real spiketime bins is larger than randomely
%shuffled ones
led_flag = max_real_mean3consBins > cutoff_99th_val;
end
