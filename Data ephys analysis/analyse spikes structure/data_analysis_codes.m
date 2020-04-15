close all;

filtered_spikes = spikes.labels;
filtered_spikes((filtered_spikes(:,2)==0),:) = [];

%assign = filtered_spikes(:,1);
assign = 103;

tuning_param = 1; % 0 = ori; 1 = spat; 2 = temp;

stim_cond_list = unique(spikes.vs_params(:, tuning_param));

%plot when visual stimulation is displayed
x = [0, spikes.vs_values(1,1), spikes.vs_values(1,1), spikes.vs_values(3,1) + spikes.vs_values(1,1), spikes.vs_values(3,1) + spikes.vs_values(1,1), spikes.vs_values(7,1)];
y = 2*[1, 1, 1.5, 1.5, 1, 1];

ymax = 20;

psth_unit(spikes, assign, stim_cond_list(1), [0 ymax]);
hold on;
plot(x, y);
psth_unit(spikes, assign, stim_cond_list(2), [0 ymax]);
hold on;
plot(x, y);
psth_unit(spikes, assign, stim_cond_list(3), [0 ymax]);
hold on;
plot(x, y);
psth_unit(spikes, assign, stim_cond_list(4), [0 ymax]);
hold on;
plot(x, y);
psth_unit(spikes, assign, stim_cond_list(5), [0 ymax]);
hold on;
plot(x, y);
psth_unit(spikes, assign, [],[]); %average over all spatial freq
hold on;
plot(x, y);

%data_analysis_ori_spontaneousfix(spikes, assign);
data_analysis_spatial_tuning(spikes, assign);
%data_analysis_temporal_tuning(spikes, assign);