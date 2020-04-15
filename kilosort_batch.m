%%%runKiloSort for many folders with same experiment name

%Assumes raw data is saved in 'E:\extarcellular\Andreanne\IntanData'
%Will save sorted data in 'E:\extarcellular\Andreanne\KilosortData' with
%same experiment name 

%also assumes electrode probe used is 'A2x32-Poly5' to change: specify in
%runKiloSort4Intan_many_files (tha variable is called amp)

exp_name = 'c57-0010-m3';

exp_folders = extract_many_folder_names(exp_name);
num_exp = length(exp_folders);

for ind = 1:num_exp
    
    current_folder = exp_folders{ind};
    runKiloSort4Intan2_many_files(current_folder);
    
end