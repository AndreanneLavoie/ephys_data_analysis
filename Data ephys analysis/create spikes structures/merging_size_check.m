%check dimensions of time.dat, analague and spikestimes format after merging by comparing to non merged
%data

merge_path = 'E:\extarcellular\Andreanne\KilosortData\vgt-0020-m1-t1-sfANDori';
single_path = 'E:\extarcellular\Andreanne\IntanData\vgt-0020-m1-t1-sf_ephys_191205_111358';
single_path2 = 'E:\extarcellular\Andreanne\IntanData\vgt-0020-m1-t1-ori_ephys_191205_121335';
single_kilo_path = 'E:\extarcellular\Andreanne\KilosortData\vgt-0020-m1-t1-sf';
single_kilo_path2 = 'E:\extarcellular\Andreanne\KilosortData\vgt-0020-m1-t1-ori';

%check time.dat format
t_merge = importIntanTime('time.dat', merge_path, 30000);
t_single = importIntanTime('time.dat', single_path, 30000);
t_single2 = importIntanTime('time.dat', single_path2, 30000);

t_merge_size = size(t_merge);
t_single_size = size(t_single);
t_single_size2 = size(t_single2);
single_addition = t_single_size(1) + t_single_size2(1);

if t_merge_size(1) - single_addition ~= 0
   fprintf('ERROR: time.dat merging incorrect, check size \n')
end


%check analogue files formats
analogue_merge = importIntanAnalog('board-ADC-01.dat', merge_path);
analogue_single = importIntanAnalog('board-ADC-01.dat', single_path);
analogue_single2 = importIntanAnalog('board-ADC-01.dat', single_path2);

t_merge_size = size(analogue_merge)
t_single_size = size(analogue_single)
t_single_size2 = size(analogue_single2)

%check spiketimes files formats
abs_spiketimes_merge = double(readNPY([merge_path '\spike_times.npy']))/double(30000); %absolute time point for each spike   
abs_spiketimes_single = double(readNPY([single_kilo_path '\spike_times.npy']))/double(30000); %absolute time point for each spike   
abs_spiketimes_single2 = double(readNPY([single_kilo_path2 '\spike_times.npy']))/double(30000); %absolute time point for each spike   

t_merge_size = size(abs_spiketimes_merge)
t_single_size = size(abs_spiketimes_single)
t_single_size2 = size(abs_spiketimes_single2)
single_addition = t_single_size(1) + t_single_size2(1)

if t_merge_size(1) - single_addition ~= 0
   fprintf('ERROR: spike_times.npy merging incorrect, check size \n')
end
%%RESOLVED: it actually makes sense that merge has a much bigger number of
%%spikes because 