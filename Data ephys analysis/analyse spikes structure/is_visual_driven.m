function [visual_activity_flag, direction_selectivity_flag] = is_visual_driven(spikes, assign, tuning)
%INPUT
%spikes construct
%assign
%tuning

%OUTPUT 
% visual_activity_flag (logic value)
% direction_selectivity_flag (temporal-nasal directions)
% DSI_over_cutoff


%% timing params
winstart = spikes.vs_params(1, 9)/1000; %trig time in sec
winwidth = spikes.vs_params(1, 7)/1000; %vis stim duration
stimTF = spikes.vs_params(1,2); %Hz
stimPeriod = 1/stimTF; %in sec

%set window ranges for baseline


num_trials = max(spikes.trials);

switch tuning 
    case 'sf'
        %set window temp->nasal and nasal->temp directions
        winBase=[0; winstart];
        winTN = winstart+[[0:stimPeriod:(winwidth-0.0001) (3/4*stimPeriod):stimPeriod:winwidth]; [(1/4*stimPeriod):stimPeriod:winwidth stimPeriod:stimPeriod:winwidth]];%T windwon of temponasal direction; each col defines one subwindow
        winNT = winstart+[(1/4*stimPeriod):stimPeriod:winwidth; (3/4*stimPeriod):stimPeriod:winwidth];%T windwon of nasotemperal direction 

        %% initialize output var for loop over all trials
        visual_activity = zeros(num_trials,3); %first row baseline, second row temp->nas, third row nas->temp;

        for k=1:num_trials
            filtered_spikes = filtspikes(spikes, 0, 'assigns', assign, 'trials', k);

            %baseline
            [baseline_fr temp]=computeFR_multiWin(filtered_spikes,winBase);
            [TN_fr temp]=computeFR_multiWin(filtered_spikes,winTN);
            [NT_fr temp]=computeFR_multiWin(filtered_spikes,winNT);

            visual_activity(k, :) = [baseline_fr TN_fr NT_fr];

        end

        %% paired t-test

        %comparing baseline to temp_nasal_dir
        [TN_h, TN_p] = ttest(visual_activity(:, 1), visual_activity(:, 2));

        %comparing baseline to nasal_temp_dir
        [NT_h, NT_p] = ttest(visual_activity(:, 1), visual_activity(:, 3));

        %determine if unit is direction selective
        [DS_h, DS_p] = ttest(visual_activity(:, 2), visual_activity(:, 3));
        visual_activity_mean = mean(visual_activity);

        %% check for nans
        TN_nan_flag = isnan(TN_h);
        NT_nan_flag = isnan(NT_h);
        DS_nan_flag = isnan(DS_h);

        if NT_nan_flag
            NT_h = 0;
        end

        if TN_nan_flag
            TN_h = 0;
        end

        if DS_nan_flag
            DS_h = 0;
        end


        %% determine output

        %output 1 if any visual driven activity (either/both TN and NT directions)
        if TN_h | NT_h
            visual_activity_flag = 1;
        else
            visual_activity_flag = 0;
        end

        %output 1 if any direction_selectivity_flag for TN direction only
        if DS_h & visual_activity_mean(2) > visual_activity_mean(3)
            direction_selectivity_flag = 1;
            rev_direction_selectivity_flag = 0;
        elseif DS_h & visual_activity_mean(2) < visual_activity_mean(3)
            direction_selectivity_flag = 0;
            rev_direction_selectivity_flag = 1;
        else
            direction_selectivity_flag = 0;
            rev_direction_selectivity_flag = 0;
        end
        
    case 'ori'
        %set window ranges for vis_stim
        winBase=[0 winstart];
        winStim=[winstart winstart+winwidth];
        
        % compute baseline fr for each trial (third var output)
        
        visual_activity = zeros(num_trials,2); %first row baseline, second stim
        for k=1:num_trials
            filtered_spikes = filtspikes(spikes, 0, 'assigns', assign, 'trials', k);

            [baseline_fr, baseline_fr_sem, baseline_fr] = computeFR(filtered_spikes, winBase);
            [stim_fr, stim_fr_sem, stim_fr] = computeFR(filtered_spikes, winStim);

            visual_activity(k, :) = [baseline_fr stim_fr];

        end
        
        %comparing baseline to stim fr
        [VDA_h, VDA_p] = ttest(visual_activity(:, 1), visual_activity(:, 2));
        VDA_nan_flag = isnan(VDA_h);
        
        if VDA_nan_flag
            VDA_h = 0;
        end
        
        %output 1 if any visual driven activity (either/both TN and NT directions)
        if VDA_h
            visual_activity_flag = 1;
        else
            visual_activity_flag = 0;
        end
        
        direction_selectivity_flag = nan;
end
end

