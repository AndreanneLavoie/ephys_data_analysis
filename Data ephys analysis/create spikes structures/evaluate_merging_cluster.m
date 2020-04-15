function evaluate_merging_cluster(spikes, assigns, tuning, tagging_flag)
%INPUT
    %spikes construct
    %assigns (clusterID of interest)
    %tuning parameter used during experiment ie: 'sf', 'tf', 'ori'
    %tagging_flag whether optogenetic tagging is used(1 if used, 0 if not)
    
%OUTPUT
    %refractory period violation for each cluster
    %tuning plot for each cluster
    %refractory period violation for merged super-cluster
    %tuning plot for merged super-cluster

%% Calculate subplot parameters based on number of assigns that will be compared
num_assigns = length(assigns);
subplot_x = 4; %4 rows of graphs when plotting; row #1 tuning; row #2 trial psth; #3 led hist gfirst #4 led hist all
subplot_y = num_assigns + 1; 
num_plots = subplot_x * subplot_y; %two plot per assign plus two plot for merged assigns 

bin_psth = 50; %ms
bin_LED = 1; %ms
bin_limit = [0 25]; %ms

% if num_assigns > 6       
%     subplot_x = 2;
%     subplot_y = ceil(num_plots/2);
% else
%     subplot_x = 1;
%     subplot_y = num_plots;
% end

% led = 1;

%does led filed exist
if  tagging_flag
    led=unique(spikes.led);%in tagging experiments, since led appear outside of VS, all trials can be used in tuning
    led(isnan(led))=[];%remove nan
else
    led=0;   % in optogenetic manipulation, only led off trials will be used in tuning.
end

%% Caluculate refractory period violation for each cluster (assign)

%returns a [assigns x 5] matrix with colums representing [assign, numVio, numspk, percent_vio, consider_bad];
ref_period_vio_ind_assign = checkRefractoryPeriodViolations_KiloSort_AL(spikes, assigns, 'select'); 

%returns a [1 x 4] vector with colums representing [numVio, numspk, percent_vio, consider_bad];
ref_period_vio_merge = checkRefractoryPeriodViolations_KiloSort_AL(spikes, assigns, 'merge'); 

%% Caluculate tuning curve for each cluster (assign)
h = figure;
set(h,'position', [50, 50, 1800, 900]);

for j=1:num_assigns
    %plot tuning curve for each assign
    subplot(subplot_x, subplot_y, j);
    [stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, assigns(j), tuning, 0,led);
    tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;
    switch tuning
    case 'sf'
        plot(stim_cond_list, tuning_fr_evoked);

    case 'tf'
        plot(stim_cond_list, tuning_fr_evoked);

    case 'ori'
        if any(tuning_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
            tuning_fr_evoked = tuning_fr_evoked-min(tuning_fr_evoked);
            fprintf('tuning_fr_evoked contains negative value, force minimum to be zero\n');
        end
        theta = stim_cond_list/180*pi;
        theta = [theta; theta(1)];
        tuning_fr_evoked = [tuning_fr_evoked; tuning_fr_evoked(1)];
        %tuning_fr_total = [tuning_fr_total; tuning_fr_total(1)];
        polar(theta, tuning_fr_evoked);
        %polar(theta, tuning_fr_total);
    end
    
    title(['Cluster ID:' num2str(assigns(j)) ' v%:' num2str(ref_period_vio_ind_assign(j, 4), '%6.2f') '('  num2str(ref_period_vio_ind_assign(j, 3), '%6.0f') ')']);
    
    
    %plot trial psth for each assign
    
    filt_spikes_unit = filtspikes(spikes, 0, 'assigns', assigns(j),'led',led); %filter by unit
    subplot(subplot_x, subplot_y, j + subplot_y); %the jth column of second row;
    title(['Cluster ID:' num2str(assigns(j)) ' LED ' num2str(led)  ' psth']);
    filt_unit_psth = psth_bl(filt_spikes_unit, bin_psth);
    %photo_tagging_param = plot_led_unit(spikes, assigns(j), tuning, 1)

    
    if tagging_flag
        subplot(subplot_x, subplot_y, j + subplot_y*2);
        [length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, assigns(j),bin_LED, bin_limit, 1, 1);
        title('LED 1st  hist')
    
        %led all hist
        subplot(subplot_x, subplot_y, j + subplot_y*3);
        [length_led_all, led_spiketimes_all, led_hist_all] = make_first_led_spiketimes(spikes, assigns(j), bin_LED, bin_limit, 0, 1);
        title('LED all hist');
        
    else
        subplot(subplot_x, subplot_y, j + subplot_y*2);
        filt_spikes_unit = filtspikes(spikes, 0, 'assigns', assigns(j),'led',1); %filter by unit
        if ~isempty(filt_spikes_unit.spiketimes)
            title(['Cluster ID:' num2str(assigns(j)) ' LED on psth']);
            filt_unit_psth = psth_bl(filt_spikes_unit, bin_psth);
        end
    end
   
end
    
%% Caluculate refractory period violation and tuning curve MERGED cluster

%MERGED tuning curve
subplot(subplot_x, subplot_y, subplot_y);
[stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, assigns, tuning, 0,led);
tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;

switch tuning
case 'sf'
    plot(stim_cond_list, tuning_fr_evoked);

case 'tf'
    plot(stim_cond_list, tuning_fr_evoked);

case 'ori'
    if any(tuning_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
        tuning_fr_evoked = tuning_fr_evoked-min(tuning_fr_evoked);
        fprintf('tuning_fr_evoked contains negative value, force minimum to be zero\n');
    end
    theta = stim_cond_list/180*pi;
    theta = [theta; theta(1)];
    tuning_fr_evoked = [tuning_fr_evoked; tuning_fr_evoked(1)];
    %tuning_fr_total = [tuning_fr_total; tuning_fr_total(1)];
    polar(theta, tuning_fr_evoked);
    %polar(theta, tuning_fr_total);
end

title(['Merged cluster:'    ' vio %:   ' num2str(ref_period_vio_merge(4), '%6.2f') '(' num2str(sum(ref_period_vio_ind_assign(:, 3))) ')' ]);


%MERGED psth
filt_spikes_unit = filtspikes(spikes, 0, 'assigns', assigns, 'led', led);    
subplot(subplot_x, subplot_y, subplot_y*2);
merge_psth = psth_bl(filt_spikes_unit, bin_psth);
title(['PSTH Merged LED: ' num2str(led)]);

%MERGED led  histogram
if tagging_flag
    %led 1st hist
    subplot(subplot_x, subplot_y, subplot_y*3);
    [length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, assigns, bin_LED, bin_limit, 1, 1);
    title('Merged LED_tagging 1st hist');
    
    %led all hist
    subplot(subplot_x, subplot_y, subplot_y*4);
    [length_led_all, led_spiketimes_all, led_hist_all] = make_first_led_spiketimes(spikes, assigns, bin_LED, bin_limit, 0, 1);
    title('Merged LED_tagging all hist');
else
    subplot(subplot_x, subplot_y, subplot_y*3);
    filt_spikes_unit = filtspikes(spikes, 0, 'assigns', assigns,'led',1); %filter by unit
    if ~isempty(filt_spikes_unit.spiketimes)
        title(['Merged LED on psth']);
        filt_unit_psth = psth_bl(filt_spikes_unit, bin_psth);
    end
end


end