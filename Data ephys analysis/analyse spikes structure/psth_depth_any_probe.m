function spikes = psth_depth_any_probe(spikes, tuning)
%% INPUT : spikes construct

%OUTPUT:   function will generate muliple new spikes fields:
    %spikes.channel - associates a channel of origine for each spiketimes
    %spikes.unique_depths - unique depth values for the probe
    %spikes.notdtn_range [num_shanks x 2] 1st column: min val, 2nd column: max val; 
%       NOTE: spikes.notdtn_range uses DSI with cut of 0.3 to determine
%       wether visual driven activity is direction selective
    % spikes.chan_depth_map_index [channelID depth_range_index in_notdtn]

% PLOT: also directly save psth_depth plot to correct kilosort path

%% create spikes.channel field
spikes.channel = spikes.assigns;
clusters = unique(spikes.assigns);
for i = 1:length(unique(spikes.assigns))
    spikes.channel(spikes.channel == clusters(i)) = spikes.clusterinfo.channel(i);
end

load('E:\extarcellular\Andreanne\chanMaps\chanMapA2x32-Poly5.mat');
spikes.unique_depths = flip(unique(ycoords)); %flip so that 0 = bottom 
num_depths = length(spikes.unique_depths);

%parameters
channels = chanMap - 1;  %convert to 0 based Intan channels
numchannel = length(channels);
shanks = kcoords; %1 = right; 2=left
num_shanks = length(unique(shanks));

%convert depth to channel number so that 1 == top channel to #of channel
% spikes.channelID = (numchannel - 1) - spikes.channel_depth;

%scaling the vertical axis to compare amplitude of activity/direction selectivity
nsum=zeros(num_depths, num_shanks);

%initializing vda and ds info
visual_activity_flag = zeros(num_depths, num_shanks);
direction_selectivity_flag = zeros(num_depths, num_shanks);
DSI = zeros(num_depths, num_shanks);

% 2 - plot psth per depth and shank
figure;
for i=1:num_depths
    
    %plot tunning curve according to depth
    %data_analysis_spatial_tuning(spikes, assign);
    
    for j=1:num_shanks
        
        %grab the identity of all channels of given depth and shank
        channels_per_depth = channels(ycoords == spikes.unique_depths(i) & kcoords == j);
        
        %identify all clusters from those channels (Assigns)
        index = nan(length(spikes.clusterinfo.channel),length(channels_per_depth));
        for n= 1:length(channels_per_depth)
            index(:, n) = spikes.clusterinfo.channel == channels_per_depth(n);
        end
        index = sum(index, 2);
        Assigns = spikes.labels(logical(index), 1);
        
        %filter spikes for those assigns
        spikes_per_depth = filtspikes(spikes, 0, 'assigns', Assigns );
        spikes_per_depth.sweeps.trials = 1: size(spikes.vs_params, 1);
        
        %plot spikes 
        subplot(num_depths, num_shanks,((i-1)*num_shanks + j));
        [hPsth, hAxes, n, centers, edges] = psth_bl(spikes_per_depth);
        nsum(i,j)=max(n);
        axis tight off;
        
        %EVALUATE DIRECTION SELECTIVITY
        if isnan(Assigns)
            continue
        else
            [visual_activity_flag(i,j), direction_selectivity_flag(i,j), DSI(i,j)] = is_visual_driven(spikes_per_depth, Assigns, tuning, 0.3);
        end
    end
end

%loop to set scaling for y
for i=1:num_depths
    ymax=max(nsum(i,:));
    for j=1:num_shanks
        subplot(num_depths, num_shanks,((i-1)*num_shanks + j));
        if ymax>0
            ylim([0 ymax]);
        end
    end
end

suptitle('RIGHT SHANK                                         LEFT SHANK')
savefig([spikes.kilosort_path 'psth_depth'])

%determine NOT/DTN range 
in_notdtn = DSI & visual_activity_flag;
spikes.notdtn_range = nan(num_shanks, 2); %eaach row is one shank; %column 1: min index, column 2: max index
for p = 1:num_shanks
    try
        spikes.notdtn_range(p, 1) = find(in_notdtn(:,p), 1, 'first');
    catch ME
        if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
            msg = ['Unable to perform assignment because the size of the left side is 1-by-1 and the size of the right side is 0-by-1.'];
            spikes.notdtn_range(p, 1) =0;
        end
    end
    
    try    
        spikes.notdtn_range(p, 2) = find(in_notdtn(:,p), 1, 'last');
        catch ME
        if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
            msg = ['Unable to perform assignment because the size of the left side is 1-by-1 and the size of the right side is 0-by-1.'];
            spikes.notdtn_range(p, 2) =0;
        end
    end
end

%store mapping between channel and depth
spikes.chan_depth_map_index = zeros(numchannel, 3);
spikes.is_notdtn = spikes.channel;
for t = 1:numchannel
    shank = kcoords(t);
    depth_ind = find(spikes.unique_depths == ycoords(t));
    in_notdtn = depth_ind >= spikes.notdtn_range(shank, 1) & depth_ind <= spikes.notdtn_range(shank, 2);
    spikes.chan_depth_map_index(t, :) = [channels(t) depth_ind in_notdtn];
    spikes.is_notdtn(spikes.channel == spikes.chan_depth_map_index(t,1)) = spikes.chan_depth_map_index(t,3);
end

%save spike construct in kilosort path
save([spikes.kilosort_path 'spikes.mat'], 'spikes');
end

