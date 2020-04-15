% omits ops.fbinary, ops.chanMap, ops.NchanTOT, ops.Nchan, ops.Nfilt
% no need to change
ops.GPU                 = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)		
ops.parfor              = 0; % whether to use parfor to accelerate some parts of the algorithm		
ops.verbose             = 1; % whether to print command line progress		
ops.showfigures         = 1; % whether to plot figures during optimization		
		
ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys' %%%data type of the input no need to change for one file per channel Intan recording		
% ops.fbinary             = 'E:\Documents\DATA\Satoru\exp_id\example_short.dat'; % will be created for 'openEphys' %%% made in mergeIntanOFPCSystem2.m		
ops.fproc               = fullfile(savePath, 'temp_wh.dat'); % 'E:\Documents\DATA\Satoru\temp_wh.dat'; % residual from RAM of preprocessed data		
ops.root                = savePath; % 'E:\Documents\DATA\Satoru'; % 'openEphys' only %%%where kilosort data (the concatenated file ops.fbinary) is saved (experiment specific folders)			
ops.fs                  = 30000;        % sampling rate		(omit if already in chanMap file) %%% change the sampling rate of intan recording to 30k
% ops.NchanTOT            = 32;           % total number of channels (omit if already in chanMap file) %%% made in mergeIntanOFPCSystem2.m	
% ops.Nchan               = 32;           % number of active channels (omit if already in chanMap file) %%% made in mergeIntanOFPCSystem2.m
% ops.Nfilt               = 128;           % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32) %%% made in mergeIntanOFPCSystem2.m    		
ops.nNeighPC            = 16; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)		
ops.nNeigh              = 24; % visualization only (Phy): number of neighboring templates to retain projections of (16)		
		% the above two parameters determine how spikes are highlighted in phy tracer veiwer
        
% options for channel whitening		%%%determine noise
ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below) %%% 'full' uses template matching, 'noSpikes' uses standard deviation; 'full' works to remove electrical noise crossing all channels, but not 'noSpikes' 	
ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1) %%% 1 means every batch of data is for whitening		
ops.whiteningRange      = Inf; %32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)		
		
% define the channel map as a filename (string) or simply an array		
% ops.chanMap             = 'C:\DATA\Spikes\Piroska\chanMap.mat'; % make this file using createChannelMapFile.m	 %%%loaded in mergeIntanOFPCSystem2.m	
ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). 		
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if a chanMap file		
		
% other options for controlling the model and optimization		
ops.Nrank               = 3;    % matrix rank of spike template model (3)  %%% no need to change, for spike detection
ops.nfullpasses         = 12; %8%6;    % number of complete passes through data during optimization (6) %%% 8 is the number of iterations for template matching; the more the better, but more time consuming		
ops.maxFR               = 50000; %20000  % maximum number of spikes to extract per batch (20000) %%% translate firing rate to this number(FR*batch time). If the result has # of spikes close to this number, increase it. make a relatively bigger number to avoid loss of spike detection	
ops.fshigh              = 300;   % frequency for high pass filtering		
% ops.fslow             = 2000;   % frequency for low pass filtering (optional)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection	%%% no need to change	
ops.scaleproc           = 200;   % int16 scaling of whitened data		                            %%% no need to change
ops.NT                  = 32*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory) 	%%% no need to change
% for GPU should be multiple of 32 + ntbuff		
		
% the following options can improve/deteriorate results. 		
% when multiple values are provided for an option, the first two are beginning and ending anneal values, 		
% the third is the value used in the final pass. 		
ops.Th               = [4 14 14];    % threshold for detecting spikes on template-filtered data ([6 12 12]) %%% first number is the thredhold (folds in standard deviation) in the first interation; the second is for the ending; the third for the final pass; the smaller first value, the more spikes to be detected.	
ops.lam              = [10 35 35];   %%%these number determines the penalty of the similarity to the spike template; Penalty: each spike is given a score based on similarity to a template (the bigger score, the worse similarity); the values here will be multiplied with the similarity score to determine penalty. Bigger number means more stringency. For thalamus where spikes from different units look similar, [10 35 35]; for cortex [7 30 30];  % large means amplitudes are forced around the mean ([10 30 30])	
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)		
ops.momentum         = 1./[20 200];  % start with high momentum and anneal (1./[20 1000])	%%% read the paper for this parameter	
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)	%%% means yes	
ops.mergeT           = 0.1;           % upper threshold for merging (.1)	%%% smaller number means more splitting and less merging	
ops.splitT           = 0.1;           % lower threshold for splitting (.1)	%%% preferrs more splitting than merging	
		
% options for initializing spikes from data		
ops.initialize      = 'fromData'; % 'fromData' or 'no'r		%%%'fromData' means real spikes; 'no' means a model spike 
ops.spkTh           = -4;      % spike threshold in standard deviations (4) -6??	%%% the sign means nothing	
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])	%%% number of time(3), channel(1) to find the peak and valley of spike waveform
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])		
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5]) %%% to avoid double accounting of spikes in template gen, mask the neighoring channels.		
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 50000;   % maximum "unique" spikes to consider (10000)	%%% number of spikes randoms picked up to make template
		
% load predefined principal components (visualization only (Phy): used for features)		
dd                  = load('PCspikes2.mat'); % you might want to recompute this from your own data	% use the original file.	
ops.wPCA            = dd.Wi(:,1:7);   % PCs 		
		
% options for posthoc merges (under construction) % This is an optional step after automatic sorting, but before manual merging/splitting. But Satoru was not satified with its performance.		
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)		
ops.epu     = Inf;		
		
ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.
