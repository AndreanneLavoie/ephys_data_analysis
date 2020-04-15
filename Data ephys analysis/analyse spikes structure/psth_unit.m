function psth_unit(spikes, unitID, bin, stimcond, yrange)

%INPUT: spikes structure & unit ID number on which psth will be run

%OUTPUT: psth plot
pretrig = spikes.vs_params(1, 9)/1000;
time_bef_stim = spikes.vs_params(1, 6)/1000;
duration =  spikes.vs_params(1, 7)/1000;
time_aft_stim =  spikes.vs_params(1,8)/1000;
vs_start = pretrig + time_bef_stim;  % absolute time at which stim starts (trig_time + time_bef)
stop = vs_start + duration + time_aft_stim; %in sec

if isempty(stimcond)
    spikes_unit = filtspikes(spikes, 0, 'assigns', unitID(:,1));
else 
    spikes_unit = filtspikes(spikes, 0, 'assigns', unitID(:,1), 'stimcond', single(stimcond));
end

h = figure;
psth_bl(spikes_unit, bin);

%set location of plot
if ~ isempty(yrange)
    axis([0, stop, yrange(1), yrange(2)])
end

set(h, 'position', [1,1,500, 500]);
end

