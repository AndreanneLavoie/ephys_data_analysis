function rel_times = rel_spiketimes(abs_times, vs_trigger)
%This function will output the relative times between trigger and spike 
%for each spike

rel_times = [];
dif_list = [];

for T = 1:length(abs_times)
    for t = 1:length(vs_trigger)
        dif = vs_trigger(t)-abs_times(T);
        dif
        if dif < 0
            dif_list = [dif_list, dif];
            dif_list
        end
    
    end
    %last value (smallest) that is not 0 will be assigned to each spike
    rel_times
    rel_times = [rel_times, dif_list(end)];
end
end
%TEST Ts are not producing negative and therefore dif_list is empty producing index error at the end
