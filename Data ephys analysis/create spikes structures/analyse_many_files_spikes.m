%run spikes creation for all files in particular path with specified
%foldernames
exp_name = 'c57-0012-m0';
data_path = {'E:\extarcellular\Andreanne\IntanData\'};
folder_names = {'c57-0012-m0-t11-sf1_ephys_190627_183510', 'c57-0012-m0-t10-sf1_ephys_190627_181730', 'c57-0012-m0-t9-sf1_ephys_190627_175747', 'c57-0012-m0-t8-sf1_ephys_190627_173459', 'c57-0012-m0-t7-sf1_ephys_190627_170830', 'c57-0012-m0-t6-sf1_ephys_190627_164223', 'c57-0012-m0-t5-sf1_ephys_190627_153824', 'c57-0012-m0-t4-tf1_ephys_190627_143428', 'c57-0012-m0-t4-sf1_ephys_190627_134451', 'c57-0012-m0-t3-sf1_ephys_190627_123742', 'c57-0012-m0-t2-tf1_ephys_190627_120954', 'c57-0012-m0-t2-sf2_ephys_190627_111907','c57-0012-m0-t2-sf1_ephys_190627_111837', 'c57-0012-m0-t1-tf1_ephys_190627_102455', 'c57-0012-m0-t1-sf1_ephys_190627_093630'}; 

% all_folders_in_path = dir(fullfile(data_path, '*'));
% folder_names = nan(length(all_folders_in_path), 1);

% for i = 1:length(all_folders_in_path)
%    
%     if length(all_folders_in_path(i).name) > 10
%         if all_folders_in_path(i).name(1:length(exp_name)) == exp_name
%             folder_names(i) = all_folders_in_path(i).name; 
%         end
%     end
% end

destination_path = 'E:\extarcellular\Andreanne\data_analysis\NOTDTN recording practice';
%[file, path, filterindex] = uigetfile('', 'Select Data Files', 'MultiSelect', 'on');
%selpath = uigetdir('E:\extarcellular\Andreanne\IntanData') %no multiselect mode :(

for i = 1:length(folder_names)
    
%     if isnana(folder_names)
%         continue
%     end
    
    %figure;
   
    full_path = [data_path{1,1} folder_names{1,i}];
    spikes = create_spikes(full_path);
    psth_depth
    
    %mkdir(['E:\extarcellular\Andreanne\data_analysis\NOTDTN recording practice\' folder_names{1,i} '\']);
    %savefig([destination_path '\' folder_names{1,i} '\psth_depth.fig']);
    
    %close all;
    
end
