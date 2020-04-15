figure;

hold on;

assign = 206;

psth_bl(filtspikes(spikes, 0, 'led', 0, 'assigns', assign));
psth_bl(filtspikes(spikes, 0, 'led', 1, 'assigns', assign), [], [], [], 'r');