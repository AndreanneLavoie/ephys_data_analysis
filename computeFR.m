function varargout = computeFR(spikes,twindow)
% function fr = computeFR(spikes,twindow)
% Computes firing rate (spikes/sec) in window defined by twindow (units,
% sec)
%
% INPUT
%   spikes: spikes struct
%   twindow: 
%
% OUTPUT
%   varargout{1} = fr_avg: Average firing in window
%   varargout{2} = fr_sem: Std error of the mean
%   varargout{3} = fr: all firing rates on trial-by-trial basis

% Created: 5/16/10 - SRO
% Modified: 10/25/10 - SRO (compute mean and sem)

% Compute window size (sec)
winsize = diff(twindow);

% Compute number of trials
ntrials = length(spikes.sweeps.trials);

% Loop through trials
for i = 1:ntrials
    trial = spikes.sweeps.trials(i);
    tempspikes = filtspikes(spikes,0,'trials',trial);
    fr(i) = (sum(tempspikes.spiketimes >= twindow(1) & ...
    tempspikes.spiketimes <= twindow(2)))/winsize;
end

if exist('fr','var')
    % Compute average firing and SEM
    fr_avg = mean(fr);
    fr_sem = std(fr)/sqrt(ntrials);
else
    fr_avg = nan;
    fr_sem = nan;
    fr = nan;
end

% Outputs
varargout{1} = fr_avg;
varargout{2} = fr_sem;
varargout{3} = fr;

