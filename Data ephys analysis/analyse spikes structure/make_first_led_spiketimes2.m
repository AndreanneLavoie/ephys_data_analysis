%%
assign = 103; 
bin = 0.5; %ms

%filter spiketimes by assign
filt_spikes = filtspikes(spikes, 0, 'assigns', assign );

%create extended matrix to allow comparision between all possible led
%times and all possible abs spiketimes
length_st = length(filt_spikes.abs_spiketimes);
length_led = length(filt_spikes.led_abstime(:,1));

abs_spikestimes_matrix = repmat(filt_spikes.abs_spiketimes, length_led, 1);%sec
abs_ledtimes_matrix = repmat(filt_spikes.led_abstime(:,1), 1, length_st); %sec

%calculate rel spiketime to ref led
diff_matrix = abs_spikestimes_matrix - abs_ledtimes_matrix;

%remove neg value
led_period_temp = filt_spikes.led_abstime(2:(end),1) - filt_spikes.led_abstime(1:(end-1),1);
led_period_temp(led_period_temp > 0.050005) = [];
led_period = mean(led_period_temp);
%remove spikes that are less tham 0 and over 0.05 sec
diff_matrix(diff_matrix < 0  | diff_matrix >= led_period) = nan; 

%led_ID = nan(size(diff_matrix));
led_pulse_num = spikes.vs_params(1,16);
led_ID_temp = 1:length_led;
led_ID = rem(led_ID_temp, led_pulse_num); %find remainder between number led events and led pulse number so that every vis stim trial has an repeating ordered number of led events (e.g. 1 to 40)
led_ID(led_ID == 0) = 40; %change remaninder 0 to id num of 40 for consistency;
led_ID(isnan(diff_matrix)) = nan;

%remove nan
led_spiketimes = diff_matrix(~ isnan(diff_matrix) );

%find min of for each led trial (ie the first spike for the cluster)
first_led_spiketimes = min(diff_matrix, [], 2);
first_led_spiketimes_logic = (diff_matrix ==first_led_spiketimes); %same dimensions as diff_matrix

%convert to ms
led_spiketimes_ms = led_spiketimes*1000;

%plot
figure
led_hist = histogram(led_spiketimes_ms, 'BinWidth', bin, 'BinLimits', [0 25]);

%check size of output var
num_spiketimes = length(led_spiketimes_ms)
num_non_nan_spiketimes = length(led_spiketimes_ms(~ isnan(led_spiketimes_ms)))

%label axis
xlabel('Latency (ms)');
ylabel('Frequency');