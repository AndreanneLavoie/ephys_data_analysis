%abs_spiketimes = [2.2 3.3 4.4 5.5 6.6 7.7 8.8 9.9];
%vstrigger_times = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
%vs_stim_order = [30 60 90 120 150 180 210 240 270];

function [rel_time, stim_cond] = make_spike_rel_and_cond(abs_spiketimes, vstrigger_times, vs_stim_order)
%this function takes as INPUT: absolute spike times, vs stimuli trigger
%times, and the absolute order of visual stim types
%RETURNS: relative spike times and spike stimuli condition

rel_time = NaN(1, length(abs_spiketimes));%each spike should have a time dif between closest stim_time
stim_cond = NaN(1, length(abs_spiketimes));%the condition of vs stimuli nearest, given spike 

for i = 1:length(abs_spiketimes)
    
    dif = vstrigger_times - abs_spiketimes(i);
        
    for n = 1:(length(dif))
            
        if dif(n) < 0 & dif(n+1) > 0
            
            rel_time(i) = dif(n); %create a list of relative spike times
            stim_cond(i) = vs_stim_order(n);
        end    
    end       
end

%rel_time
%stim_cond