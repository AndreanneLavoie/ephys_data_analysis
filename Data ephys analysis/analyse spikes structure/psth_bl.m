function [varargout] = psth_bl(spikes,binsize,hAxes,bsmooth,colorcode,linwid,plotting)
% function [varargout] = psth(spikes,binsize,hAxes,bsmooth,colorcode,linwid)
%
% INPUTS
%   spiketimes:
%   binsize:in ms
%   hAxes:
%   bsmooth:
%
% OUTPUTS
%   varargout{1} = hPsth;
%   varargout{2} = hAxes;
%   varargout{3} = n;
%   varargout{4} = centers;
%   varargout{5} = edges;

% Created:  3/14/10 - SRO
% Modified: 5/14/10 - SRO
%           6/8/10 - SRO

%Modified by BL 6/25/13
if nargin < 2
    binsize = 50;
    hAxes = gca;
    bsmooth = 0;
end

% Use current axes if hAxes not supplied
if nargin < 3
    hAxes = gca; %
    bsmooth = 0;
end

if nargin < 4
    bsmooth = 0;
end

if nargin <5
    colorcode=[0 0 0];
end

if nargin < 6
    linwid=1.5;
end

if nargin<7
    plotting=1;
end

if isempty(binsize)
    binsize = 50;
end
if isempty(hAxes)
    hAxes=gca;
end
if isempty(bsmooth)
    bsmooth = 0;
end
if isempty(colorcode)
    colorcode=[0 0 0];
end
if isempty(linwid)
    linwid=1.5;
end

% Set duration and number of trials](in s)
start_time = 0;
vs_start_time = ((spikes.vs_params(1,9)/1000 + spikes.vs_params(1,6))/1000); %trigger time + t_bef
duration = spikes.vs_params(1,7)/1000; 
t_aft = (spikes.vs_params(1,8))/1000;  
led_time = (spikes.vs_params(1,13)) * (spikes.vs_params(1,16))/1000;
end_time = vs_start_time + duration + t_aft + led_time; %start_time is neg so will be substracted from duration

numtrials = length(1: size(spikes.vs_params, 1));

% Set spiketimes
spiketimes = spikes.spiketimes;

% Convert binsize from ms to s
binsize = binsize/1000;

% Get counts
edges = start_time:binsize:end_time;
n = histc(spiketimes,edges);
n = n/numtrials/binsize;

if all(isnan(n))
    n = 0;
end

% Compute center of bins
centers = edges + diff(edges(1:2))/2;

% Last point of n contains values falling on edge(end) -- usually zero
if plotting
    if bsmooth
        hPsth = line('XData',centers(1:end-1),'YData',smooth(n(1:end-1),3),...
            'Parent',hAxes,'LineWidth',linwid,'Color',colorcode);
    else
        hPsth = line('XData',centers(1:end-1),'YData',n(1:end-1),...
            'Parent',hAxes,'LineWidth',linwid,'Color',colorcode);
    end

    % Set default axes properties
    if sum(n) == 0
        maxN = 0.1;
    elseif isempty(n)
        maxN = 0.1;
    else
        maxN = max(n);
    end
    axis([start_time end_time 0 (maxN+5)]);
    set(hAxes,'TickDir','out','FontSize',8);
    xlabel(hAxes,'seconds');
    ylabel(hAxes,'spikes/s');
else
    hPsth=[];
end
% Outputs
varargout{1} = hPsth;
varargout{2} = hAxes;
varargout{3} = n;
varargout{4} = centers;
varargout{5} = edges;