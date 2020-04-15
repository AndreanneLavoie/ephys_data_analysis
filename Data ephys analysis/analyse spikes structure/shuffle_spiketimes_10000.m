function led_flag = shuffle_led_spiketimes_10000(spikes, assign)

%parameters
shuffling_number = 10000; 
%assign = 274;
bin = 1; %ms
bin_limit = [0 100]; %ms
num_bins = bin_limit(2) - bin_limit(1);

%calculate led_spiketimes for assign
subplot(1,2,1)
[length_led, led_spiketimes, led_hist] = make_first_led_spiketimes(spikes, assign, bin, bin_limit, 1, 1);
title('Real data');

%spike firing rate in each 1ms bin
real_fr_bins = led_hist.Values;
%SHOULD THIS BE A MOVING MEAN??? or should i add all three values
%together???
real_mean3consBins = movmean(real_fr_bins, 3); %mean of three consecutive 1 ms bins
max_real_mean3consBins = max(real_mean3consBins);
real_sum_3consbins = movsum(real_fr_bins, 3);
max_real_sum_3consBins = max(real_sum_3consbins);

%remove nan vales (nan are assigned to spiketimes that are smaller than 0 and larger than led period)
led_spiketimes(isnan(led_spiketimes)) = [];
num_spikes = length(led_spiketimes);

%% 10 000 shuffling
max_shuffled_sum_3consBins = nan(shuffling_number, 1);

%generate random number
%formula will create array of random numbers drawns form a uniform distribution with the open interval [min_range, max_range] (line 20)
for i=1:shuffling_number
    rand_range = [0 100];
    random_num = (rand_range(2)-rand_range(1)).*rand(num_spikes, 1) + rand_range(1);  

    %add random numbers to led_spiketimes
    shuffled_led_spiketimes = led_spiketimes + random_num;

    %plot
    subplot(1,2,2)
    %shuffled_led_spiketimes_ms = shuffled_led_spiketimes*1000;
    shuffled_led_hist = histogram(shuffled_led_spiketimes, 'BinWidth', bin, 'BinLimits', bin_limit);
    title('Shuffled data');

    %save the firing rate of each bin
    shuffled_fr_bins = shuffled_led_hist.Values;
    shuffled_mean3consBins = movmean(shuffled_fr_bins, 3); %mean of three consecutive 1 ms bins
    max_shuffled_mean3consBins = max(shuffled_mean3consBins);
    shuffled_sum_3consbins = movsum(shuffled_fr_bins, 3);
    max_shuffled_sum_3consBins(i) = max(shuffled_sum_3consbins);

    %label axis
    xlabel('Latency (ms)');
    ylabel('# spikes');
    ylim([0 500]);

end
%calculate 99.9th percentile cutoff value
sort_shuffled_max = sort(max_shuffled_sum_3consBins);
cutoff_99th_ind = round(shuffling_number*0.999);
cutoff_99th_val = sort_shuffled_max(cutoff_99th_ind);

%determine if sum of top 3 real spiketime bins is larger than randomely
%shuffled ones
led_flag = max_real_sum_3consBins > cutoff_99th_val;
real_over_shuffle_ratio = max_real_sum_3consBins/cutoff_99th_val;

end