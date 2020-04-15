function spikes = is_visual_driven_spikes(spikes, tuning)
%% this function calls is_visual_driven for all good spikes and saves as new field in spikes construct
%INPUT
%spikes construct


%OUTPUT 
%spikes with new .is_visual_driven field [cluster is_vis_driven]; logic value 


% identify good units
good_units = spikes.labels(spikes.labels(:,2) == 2);
num_units = length(good_units);

visual_activity = nan(num_units, 2);

    for i = 1:num_units
        [visual_activity_flag, direction_selectivity_flag] = is_visual_driven(spikes, good_units(i), tuning);
        visual_activity(i,:) = [good_units(i) visual_activity_flag];
    end
    
    spikes.is_visual_driven = visual_activity;

    
end