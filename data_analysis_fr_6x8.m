%calculate the firing rate for each receptive field 

%receptive field parameters;
num_col = 8;
num_row = 6;

%time window parameters
twin_start = 0;
twin_stop = 0.3; %in sec
twin_baseline = -0.3;

%generate recpetive field index matrix (of size row * col); %on = white;
%black = off
stim_ind_map = 0:(num_col*num_row -1);
stim_ind_map = reshape(stim_ind_map, 6, 8);
stim_ind_map = [stim_ind_map; stim_ind_map + 48];

%place holder
firing_rates = stim_ind_map; 

%first do white only; use filter function to extract all the stim cond from
%
for i = 1:size(stim_ind_map, 1)
    for j =  1:size(stim_ind_map, 2)
    
        filtered_spikes = filtspikes(spikes, 0, 'stimcond', stim_ind_map(i,j));
        
        %compute firing rate for each stim cond and substract baseline
        %firing rate
        temp_fr = computeFR(filtered_spikes, [twin_start twin_stop]) - computeFR(filtered_spikes, [twin_baseline twin_start]);
        firing_rates(i, j) = temp_fr;
    end
end

plottingcolorRF3AL(firing_rates, 1, 1, 1, 0)

