%% SF

%psth per stimcond
load('E:\extarcellular\Andreanne\KilosortData\vgt-0024-m0-t4-sf\spikes.mat')
stim_cond = unique(spikes.vs_params(:,1));
num_stim = length(stim_cond);

led_unit = 259;
non_led_unit = 276;
units = [led_unit non_led_unit];
num_units = length(units);

colour = {'c', 'k'};
ymax = 8;

h = figure;
hold on

counter = 1;
for i = 1:length(units) %each unit
    for j = 1:num_stim
        
        filter_spikes = filtspikes(spikes, 0, 'assigns', units(i), 'stimcond', stim_cond(j));
        subplot(num_units, num_stim, counter)
        psth_bl(filter_spikes, [], [], [],colour{i})
        set(gca, 'FontSize', 10)
        counter = counter + 1;
        ylim([0 ymax])
        
        if i == 1
            title(num2str(stim_cond(j))) 
        end
    end
%     counter = 1; %if want overlaping
end





%% sf tuning 
figure;
hold on
%subplot(2,2,2)
[stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, non_led_unit, 'sf', 0, 1);
tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;
tuning_fr_evoked = tuning_fr_evoked/max(tuning_fr_evoked);
plot(stim_cond_list, tuning_fr_evoked);

[stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, led_unit, 'sf', 0, 1);
tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;
tuning_fr_evoked = tuning_fr_evoked/max(tuning_fr_evoked);
plot(stim_cond_list, tuning_fr_evoked);

xlabel('spatial frequency');
ylabel('firing rate (spikes/s)');
set(gca, 'FontSize', 20)
legend({'non led-driven unit', 'led-driven unit'}, 'FontSize', 20)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);


%% broad, narrow, spontaneous
load('E:\extarcellular\Andreanne\KilosortData\vgt-0009-m0-t2-ori\spikes.mat')
led_unit = 129; 
broad_unit = 56; 
non_led_unit = 265;

ymax = 850;

figure;
%subplot(1,3,1)
[length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, non_led_unit, 1, [0 25], 1, 1);
ylim([0 ymax])
set(gca, 'FontSize', 15)
set(gcf,'renderer','Painters')
figure;
%subplot(1,3,2)
[length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, led_unit, 1, [0 25], 1, 1);
ylim([0 ymax])
set(gca, 'FontSize', 15)
set(gcf,'renderer','Painters')
figure;
%subplot(1,3,3)
[length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, broad_unit, 1, [0 25], 1, 1);
ylim([0 ymax])
set(gcf,'renderer','Painters')
set(gca, 'FontSize', 15)

%set(gca, 'FontSize', 20)
%legend({'non led-driven unit', 'led-driven unit'}, 'FontSize', 20)


%% DS

%led hist
% load('E:\extarcellular\Andreanne\KilosortData\vgt-0024-m0-t5-ori\spikes.mat')
% led_unit = 9; 
% non_led_unit = 1;

load('E:\extarcellular\Andreanne\KilosortData\vgt-0024-m0-t4-ori\spikes.mat')
led_unit = 3; 
non_led_unit = 2;

%psth per stimcond
stim_cond = unique(spikes.vs_params(:,10));
num_stim = length(stim_cond);


units = [led_unit non_led_unit];
num_units = length(units);

colour = {'c', 'k'};
ymax = [10 10];

h = figure;
hold on

counter = 1;
for i = 1:length(units) %each unit
    for j = 1:num_stim
        
        filter_spikes = filtspikes(spikes, 0, 'assigns', units(i), 'stimcond', stim_cond(j));
        subplot(num_units, num_stim, counter)
        psth_bl(filter_spikes, [], [], [],colour{i})
        set(gca, 'FontSize', 10)
        counter = counter + 1;
        ylim([0 ymax(i)])
        
        if i == 1
            title(num2str(stim_cond(j))) 
        end
    end
%     counter = 1; %if want overlaping
end

% %% 
% h = figure;
% hold on
% %subplot(2,2,3)
% hold on
% [length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, non_led_unit, 1, [0 25], 1, 1);
% [length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, led_unit, 1, [0 25], 1, 1);
% 
% % set(gca, 'FontSize', 20)
% % legend({'non led-driven unit', 'led-driven unit'}, 'FontSize', 20)

%% ori tuning
figure;

[stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, non_led_unit, 'ori', 0, 1);
tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;


if any(tuning_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
    tuning_fr_evoked = tuning_fr_evoked-min(tuning_fr_evoked);
    fprintf('tuning_fr_evoked contains negative value, force minimum to be zero\n');
end
tuning_fr_evoked = tuning_fr_evoked/max(tuning_fr_evoked);

theta = stim_cond_list/180*pi;
theta = [theta; theta(1)];
tuning_fr_evoked = [tuning_fr_evoked; tuning_fr_evoked(1)];
polar(theta, tuning_fr_evoked);

hold on


[stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, led_unit, 'ori', 0, 1);
tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;
if any(tuning_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
    tuning_fr_evoked = tuning_fr_evoked-min(tuning_fr_evoked);
    fprintf('tuning_fr_evoked contains negative value, force minimum to be zero\n');
end
tuning_fr_evoked = tuning_fr_evoked/max(tuning_fr_evoked);
theta = stim_cond_list/180*pi;
theta = [theta; theta(1)];
tuning_fr_evoked = [tuning_fr_evoked; tuning_fr_evoked(1)];
        %tuning_fr_total = [tuning_fr_total; tuning_fr_total(1)];
polar(theta, tuning_fr_evoked);
        %polar(theta, tuning_fr_total);

% set(gca, 'FontSize', 20)
% legend({'non led-driven unit', 'led-driven unit'}, 'FontSize', 15)
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);

%% SF POPULATION CUMMULATIVE PLOT
load('E:\extarcellular\Andreanne\data_analysis\photo-tag\population tuning\sf_summary_led24-Feb-2020.mat')
temp_low = preferred_stim_compile_low(:,2);
temp_high = preferred_stim_compile_high(:,2);

temp_low(abs(temp_low - 0.0400)<0.0000001) = 1;
temp_low(abs(temp_low - 0.0800)<0.0000001) = 2;
temp_low(abs(temp_low - 0.1600)<0.0000001) = 3;
temp_low(abs(temp_low - 0.3200)<0.0000001) = 4;
temp_low(abs(temp_low - 0.4500)<0.0000001) = 5;
temp_high(abs(temp_high - 0.0400)<0.0000001) = 1;
temp_high(abs(temp_high - 0.0800)<0.0000001) = 2;
temp_high(abs(temp_high - 0.1600)<0.0000001) = 3;
temp_high(abs(temp_high - 0.3200)<0.0000001) = 4;
temp_high(abs(temp_high - 0.4500)<0.0000001) = 5;

figure;
hold on
ecdf(temp_low)
ecdf(temp_high)
        

