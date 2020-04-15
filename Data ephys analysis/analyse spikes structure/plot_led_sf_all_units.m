%this script will plot all of the psth, led histogram and sf tuning for each unit in
%labels and then calculate mean and std for each photo-tagging parameter in
%PHOTO_TAGG_DATA structure 

%close all;
mkdir([spikes.kilosort_path 'good_unit']);
cd([spikes.kilosort_path 'good_unit']);
flag_plot = 1; %turn on or off plotting
labels = spikes.labels;
good_units = labels(labels(:,2)== 2, 1); %good = 2 in spikes.labels(:, 2)
photo_tagg_data.unit = good_units;

%initalize variables with correct size
num_labels = length(good_units);
photo_tagg_data.labels = good_units;
photo_tagg_data.jitters = good_units;
photo_tagg_data.mean_latency = good_units;
photo_tagg_data.reliability = good_units;
photo_tagg_data.first_over_all = good_units;
photo_tagg_data.background_fr_bef_led = good_units;
photo_tagg_data.mean_led_fr = good_units;
photo_tagg_data.led_over_background = good_units;
photo_tagg_data.led_vs_back_reject_null = good_units;
photo_tagg_data.led_vs_back_p = good_units;
photo_tagg_data.vis_driven_flag = good_units;
photo_tagg_data.led_driven_flag = good_units;
photo_tagg_data.est_width = good_units;
photo_tagg_data.led_driven_ttest = good_units;
photo_tagg_data.mean_latency_3bins = good_units;


%could add a filter to look only at labels that are not noise or only good

%% This loop will plot the psth, led histogram and sf tuning for each unit in
%labels
for j = 1:num_labels
    
    photo_tagging_param = plot_led_unit(spikes, good_units(j), 'sf', flag_plot);
    saveas(gcf, [num2str(good_units(j)) '.jpeg']);
    
    photo_tagg_data.jitters(j) = photo_tagging_param.jitter;
    photo_tagg_data.mean_latency(j) = photo_tagging_param.mean_latency;
    photo_tagg_data.reliability(j) = photo_tagging_param.reliability;
    photo_tagg_data.first_over_all(j) = photo_tagging_param.first_over_all;
    photo_tagg_data.led_over_background(j) = photo_tagging_param.led_over_background;
    photo_tagg_data.background_fr_bef_led(j) = mean(photo_tagging_param.background_fr_led_pulse);
    photo_tagg_data.mean_led_fr(j) = mean(photo_tagging_param.led_pulse_fr);
    photo_tagg_data.led_driven_ttest(j) = photo_tagging_param.led_driven_ttest_h;
    
    photo_tagg_data.vis_driven_flag(j) = is_visual_driven(spikes, good_units(j));
    [photo_tagg_data.led_driven_flag(j) photo_tagg_data.mean_latency_3bins(j)]= shuffle_led_spiketimes_10000(spikes, good_units(j));
    photo_tagg_data.est_width(j) = calc_jitters(spikes, good_units(j));
    
    %save figure
    saveas(gcf, ['cluster ' num2str(good_units(j)) '.jpeg']);
    
end

%%
%Calculate mean of mean and STD of means
 photo_tagg_data.mean_jitters =  mean(photo_tagg_data.jitters(~ isnan(photo_tagg_data.jitters)));
 photo_tagg_data.std_jitters =  std(photo_tagg_data.jitters(~ isnan(photo_tagg_data.jitters)));
 photo_tagg_data.grand_mean_latency = mean(photo_tagg_data.mean_latency(~ isnan(photo_tagg_data.mean_latency)));
 photo_tagg_data.std_mean_latency = std(photo_tagg_data.mean_latency(~ isnan(photo_tagg_data.mean_latency)));
 photo_tagg_data.mean_reliability =  mean(photo_tagg_data.reliability(~ isnan(photo_tagg_data.reliability)));
 photo_tagg_data.std_reliability =  std(photo_tagg_data.reliability(~ isnan(photo_tagg_data.reliability)));
 photo_tagg_data.mean_led_over_background =  mean(photo_tagg_data.led_over_background(~ isnan(photo_tagg_data.led_over_background)));
 photo_tagg_data.std_led_over_background =  std(photo_tagg_data.led_over_background(~ isnan(photo_tagg_data.led_over_background)));
 
 

 %% SAVE 
  
save([spikes.kilosort_path 'photo_tagg_data.mat'], 'photo_tagg_data');

%% Plot

figure;

subplot(1, 6, 1);
histogram(photo_tagg_data.mean_latency*1000, 15);
xlabel('Mean Latency (ms)');
ylabel('# units');
mean_line_val = photo_tagg_data.grand_mean_latency*1000;
hold on;
line([mean_line_val, mean_line_val], ylim, 'LineWidth', 2, 'Color', 'r');

subplot(1, 6, 2);
histogram(photo_tagg_data.jitters*1000, 15);
xlabel('Jitters (ms)');
ylabel('# units');
mean_line_val = photo_tagg_data.mean_jitters*1000;
hold on;
line([mean_line_val, mean_line_val], ylim, 'LineWidth', 2, 'Color', 'r');

subplot(1, 6, 3);
histogram(photo_tagg_data.reliability*100, 15);
xlabel('Reliability (%)');
ylabel('# units');
mean_line_val = photo_tagg_data.mean_reliability*100;
hold on;
line([mean_line_val, mean_line_val], ylim, 'LineWidth', 2, 'Color', 'r');

subplot(1, 6, 4);
histogram(photo_tagg_data.led_over_background, 15);
xlabel('LED fr / background fr ');
ylabel('# units');
mean_line_val = photo_tagg_data.mean_led_over_background;
hold on;
line([mean_line_val, mean_line_val], ylim, 'LineWidth', 2, 'Color', 'r');


subplot(1, 6, 5);
scatter(photo_tagg_data.jitters*1000, photo_tagg_data.mean_latency*1000)
xlabel('Jitter (ms)')
ylabel('Mean Latency (ms)')
title('Jitter vs Latency')

figure;
histogram(photo_tagg_data.est_width, 20);
xlabel('Estimated histogram width')
ylabel('# units')

%% ARBATRARY CUTOFF FOR PHOTO-TAGGED CELLS 
%used observation of what I imagined a good led-tied cluster would look
%like and then worked their way back
%used jitter < 2.2
%mean latency < 10ms

%photo_tagg_data.led_tagged_units = photo_tagg_data.jitters*1000 < 2.2 & photo_tagg_data.mean_latency*1000 < 10;
%photo_tagg_data.labels(led_tagged_units)

