function [spikes] = psth_depth_any_probe_liu(spikes, tuning, probe)
%% INPUT : 
%spikes construct after opening phy and saving clusterinfo.tsv 
%tuning ('sf', 'ori' or 'tf'),
%probe ('chanMapA1x32Poly3.mat' or 'chanMapA2x32-Poly5.mat')

%OUTPUT:   function will generate muliple new spikes fields:
    %spikes.channel (if not already generated) - associates a channel of origine for each spiketimes
    %unique_depths - unique depth values for the probe
    %notdtn_range [num_shanks x 2] 1st column: min val, 2nd column: max val; 
%       NOTE: spikes.notdtn_range uses DSI with cut of 0.3 to determine
%       wether visual driven activity is direction selective
    % chan_depth_map_index [channelID depth_range_index in_notdtn]

% PLOT: also directly save psth_depth plot to correct kilosort path

%% create spikes.channel field
load(['E:\extarcellular\Andreanne\chanMaps\' probe])
if ~ isfield(spikes, 'channel')
    spikes.channel = spikes.assigns;
    clusters = spikes.clusterinfo.id;
    for i = 1:length(clusters)
        spikes.channel(abs(spikes.assigns - clusters(i)) < 1e-6) = Intan_chan_num(chanMap(spikes.clusterinfo.channel(i)+1)); %both phy channel info and Intan channels are 0 based, but need to be 1 based for indexing
    end
end
spikes.probe = probe;

unique_depths = flip(unique(ycoords)); %flip so that 0 = bottom 
num_depths = length(unique_depths);

%parameters
channels = Intan_chan_num;  %convert to 0 based Intan channels
numchannel = length(channels);
shanks = kcoords; %1 = right; 2=left
num_shanks = length(unique(shanks));

%convert depth to channel number so that 1 == top channel to #of channel
% spikes.channelID = (numchannel - 1) - spikes.channel_depth;

%scaling the vertical axis to compare amplitude of activity/direction selectivity
nsum=zeros(num_depths, num_shanks);

%initializing vda and ds info
% visual_activity_flag = zeros(num_depths, num_shanks);
% direction_selectivity_flag = zeros(num_depths, num_shanks);
% DSI = zeros(num_depths, num_shanks);

% 2 - plot psth per depth and shank
figure('Renderer', 'painters', 'Position', [300 500 1400 300]);
%to maximize use of space: row = 1 shank; column = 1 depth
for j=1:num_shanks
    
    %plot tunning curve according to depth
    %data_analysis_spatial_tuning(spikes, assign);
    
    for i=1:num_depths
        
        %grab the identity of all channels of given depth and shank
        channels_per_depth = channels(chanMap(abs(ycoords-unique_depths(i))<1e-6 & kcoords == j));
        
        %filter spikes for those clusters from those channels
        spikes_per_depth = filtspikes(spikes, 0, 'channel', channels_per_depth);
        spikes_per_depth.sweeps.trials = 1: size(spikes.vs_params, 1);
        
        %plot spikes 
        subplot(num_shanks,num_depths,((j-1)*num_depths + i));
        switch tuning
            case 'sf'
                [hPsth, hAxes, n, centers, edges] = psth_bl(spikes_per_depth);
                nsum(i,j)=max(n);
                axis tight off;
                title(num2str(i)); %title(num2str(channels_per_depth));
                if i == 1
                    ylabel(['Shank ' num2str(j)]);
                end
            case 'ori'
                try
                    [stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, unique(spikes_per_depth.assigns), 'ori', 0,1);
                    tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;
                    if any(tuning_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
                        tuning_fr_evoked = tuning_fr_evoked-min(tuning_fr_evoked);
                        fprintf('tuning_fr_evoked contains negative value, force minimum to be zero\n');
                    end
                    theta = stim_cond_list/180*pi;
                    theta = [theta; theta(1)];
                    tuning_fr_evoked = [tuning_fr_evoked; tuning_fr_evoked(1)];
                    polarplot(theta, tuning_fr_evoked);
                    Ax=gca;
                    %Ax.FontSize = 5;
                    Ax.ThetaGrid = 'off';
                    Ax.RGrid = 'off';
                    Ax.RTickLabel = []; 
                    Ax.ThetaTickLabel = [];
                    
                    if i == num_depths & j ~= 1
                         title(['Shank ' num2str(j) ]); 
                    end
                    
                    if j == 1
                        title(num2str(i)); %title(num2str(channels_per_depth));
                    end

                catch ME %if not much activity in section of probe which messes up compute fr, it will skip it
                    if strcmp(ME.identifier,'Error using subplot (line 327).\n Index exceeds number of subplots.')
                        continue
                    elseif strcmp(ME.identifier,'The logical indices contain a true value outside of the array bounds.')
                        continue
                    end
                
                end

%         %EVALUATE DIRECTION SELECTIVITY
%         if isnan(unique(spikes_per_depth.assigns))
%             continue
%         else
%             [visual_activity_flag(i,j), direction_selectivity_flag(i,j), DSI(i,j)] = is_visual_driven(spikes_per_depth, unique(spikes_per_depth.assigns), tuning, 0.3);
%         end
    end
    end


%loop to set scaling for y
% switch tuning
%     case 'sf'
%         for i=1:num_depths
%             ymax=max(nsum(i,:));
%             for j=1:num_shanks
%                 subplot(num_depths, num_shanks,((i-1)*num_shanks + j));
%                 if ymax>0
%                     ylim([0 ymax]);
%                 end
%             end
%         end
% end
end

suptitle('Depth')

%save with different name depending if its psth or tuning curve
switch tuning
    case {'sf', 'tf' }
        savefig([spikes.kilosort_path 'psth_depth'])
    case 'ori'
        savefig([spikes.kilosort_path 'tuning_depth'])
end

%determine NOT/DTN range 
% in_notdtn = DSI & visual_activity_flag;
% notdtn_range = nan(num_shanks, 2); %eaach row is one shank; %column 1: min index, column 2: max index
% for p = 1:num_shanks
%     try
%         notdtn_range(p, 1) = find(in_notdtn(:,p), 1, 'first');
%     catch ME
%         if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
%             msg = ['Unable to perform assignment because the size of the left side is 1-by-1 and the size of the right side is 0-by-1.'];
%             notdtn_range(p, 1) =0;
%         end
%     end
%     
%     try    
%         notdtn_range(p, 2) = find(in_notdtn(:,p), 1, 'last');
%         catch ME
%         if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
%             msg = ['Unable to perform assignment because the size of the left side is 1-by-1 and the size of the right side is 0-by-1.'];
%            notdtn_range(p, 2) =0;
%         end
%     end
% end
% 
% %store mapping between channel and depth
% chan_depth_map_index = zeros(numchannel, 3);
% spiketime_is_notdtn = spikes.channel;
% for t = 1:numchannel
%     shank = kcoords(t);
%     depth_ind = find(unique_depths == ycoords(t));
%     in_notdtn = depth_ind >= notdtn_range(shank, 1) & depth_ind <= notdtn_range(shank, 2);
%     chan_depth_map_index(t, :) = [channels(t) depth_ind in_notdtn];
%     spiketime_is_notdtn(spikes.channel == chan_depth_map_index(t,1)) = chan_depth_map_index(t,3);
% end
% 
% %save spike construct in kilosort path
save([spikes.kilosort_path 'spikes.mat'], 'spikes');
end

