%%% psth by channel depth 
%path = '';

%load spikes structure 
%spikes = load(path);

%parameters
channelsID = 1:16;  %GRAB ACTUAL NUM OF CHANNELS FROM INFO FILE; all channel numbers %%0 is BOTTOM
numchannel = length(channelsID);
stimcond = spikes.labveiw_param.unipara;
stimnum = length(stimcond);

%convert depth to channel number so that 1 == top channel to #of channel
spikes.channelID = (numchannel - 1) - spikes.channel_depth;

%filter out noise 
filtered_spikes = spikes.labels;
filtered_spikes((filtered_spikes(:,2)==0),:) = [];

%scaling the vertical axis to compare amplitude of activity/direction selectivity
nsum=zeros(numchannel, stimnum);

figure;
for i=1:numchannel
    
    %plot tunning curve according to depth
    %data_analysis_spatial_tuning(spikes, assign);
    
    for j=1:stimnum
        %filter spikes
        spikes_per_depth_cond = filtspikes(spikes, 0, 'assigns', filtered_spikes(:,1), 'channelID', i, 'stimcond', stimcond(j) );
        spikes_per_depth_cond.sweeps.trials = 1: size(spikes.vs_params, 1);
        
        %plot spikes 
        subplot(numchannel, stimnum,((i-1)*stimnum + j));
        [hPsth, hAxes, n, centers, edges] = psth_bl(spikes_per_depth_cond);
        nsum(i,j)=max(n);
        axis tight off;
    end
end

%loop to set scaling for y
for i=1:numchannel
    ymax=max(nsum(i,:));
    for j=1:stimnum
        subplot(numchannel, stimnum,((i-1)*stimnum + j));
        if ymax>0
            ylim([0 ymax]);
        end
    end
end
