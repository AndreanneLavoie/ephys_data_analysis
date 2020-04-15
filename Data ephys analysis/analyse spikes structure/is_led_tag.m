%IS CLUSTER DRIVEN BY LED?

function cluster_led = is_led_tag(spikes, assign)
cluster_led.assigns = assign;
cluster_led.binsize = 0.5; %ms
% Convert binsize from s to ms
binsize = cluster_led.binsize/1000;

%num of led events
cluster_led.numtrials = length(spikes.led_abstime(:,1));

%calculate led period
led_period_temp = spikes.led_abstime(2:(end),1) - spikes.led_abstime(1:(end-1),1);
led_period_temp(led_period_temp > 0.050005) = [];
cluster_led.led_period = mean(led_period_temp);

cluster_led.range = [-0.02 cluster_led.led_period]; %sec

start_time = range(1);
end_time = range(2);
edges = start_time:binsize:end_time;

%initalize size of variable
cluster_led.spiketimes = nan(length(cluster_led.assigns), cluster_led.numtrials);
cluster_led.freq_count = nan(length(cluster_led.assigns), cluster_led.numtrials);
cluster_led.mean_baseline = nan(length(cluster_led.assigns), 1);
cluster_led.sd_baseline = nan(length(cluster_led.assigns), 1);
cluster_led.cut_off = nan(length(cluster_led.assigns), 1);
cluster_led.mean_led = nan(length(cluster_led.assigns), 1);
cluster_led.led_over_baseline = nan(length(cluster_led.assigns), 1);
cluster_led.led_fast_tag = nan(length(cluster_led.assigns), 1);

cluster_led.jitter = nan(length(cluster_led.assigns), 1);
cluster_led.reliability = nan(length(cluster_led.assigns), 1);

%%%%COULD MAKE A LOOP TO GO OVER EACH CLUSTER
for i=1:length(cluster_led.assigns)

    % Set spiketimes
    cluster_led.spiketimes(i,:) = make_first_led_spiketimes(spikes, cluster_led.assigns(i), cluster_led.binsize);

    % Get counts
    n = histc(cluster_led.spiketimes(i),edges);
    cluster_led.freq_count(i,:) = n/cluster_led.numtrials/binsize;

    %compare led spikes average to baseline
    cluster_led.mean_baseline(i) = mean(cluster_led.freq_count(cluster_led.freq_count(i,:) < 0));
    cluster_led.sd_baseline(i) = std(cluster_led.freq_count(cluster_led.freq_count(i,:) < 0));
    cluster_led.cut_off(i) = abs(2*cluster_led.sd_baseline(i));
    cluster_led.mean_led(i) = mean(cluster_led.freq_count(cluster_led.freq_count(i,:) > 0));
    cluster_led.jitter(i) = std(cluster_led.freq_count(cluster_led.freq_count(i,:) > 0));
    cluster_led.led_over_baseline(i) = abs(cluster_led.mean_led(i))/abs(cluster_led.mean_baseline(i));
    cluster_led.led_fast_tag(i) = cluster_led.mean_led(i)  > cluster_led.cut_off(i);
    cluster_led.reliability(i) = length(cluster_led.spiketimes(i,:))/cluster_led.numtrials;

end

plot3(cluster_led.reliability, cluster_led.jitter, cluster_led.led_over_baseline);
plot3(cluster_led.reliability, cluster_led.jitter, cluster_led.led_fast_tag);

end