function rem_channels = depth2chan(spikes, depth_by_shank_matrix, include_mediocre_flag)
    
%% grab all Intan channel numbers that are at a particular depth and shank 

%INPUT: 
%1) spikes construst after phy and psth_depth_any_probe_liu codes have been
%            run (must have a spikes.probe field to pull out chanMap info)
%2) matrix [(depth shank mediocre)x n_excluded-chan] that need to excluded in data analysis

%OUTPUT:
%vector [channel x 1]; all channels at all depths/shanks 

%%
    if nargin<3
       include_mediocre_flag = 0;
    end

    load(['E:\extarcellular\Andreanne\chanMaps\' spikes.probe])
    unique_depths = flip(unique(ycoords));
    
    size_depth_by_shank = size(depth_by_shank_matrix);
    
    rem_channels = [];
   
    for n = 1: size_depth_by_shank(1) %for each row, which is each depths
        
        depth = depth_by_shank_matrix(n, 1);
        shank = depth_by_shank_matrix(n, 2);
        mediocre = depth_by_shank_matrix(n, 3);
        
        if include_mediocre_flag & mediocre == 1 %do not add these channels to list that will be removed
            continue
        else
            channels_per_depth = Intan_chan_num(chanMap(abs(ycoords-unique_depths(depth))<1e-6 & kcoords == shank));
            rem_channels = [rem_channels channels_per_depth];
        end
        
        
    end
    
end