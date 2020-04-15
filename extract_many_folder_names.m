function exp_folders = extrcat_many_folder_names(exp_name)
%%% raw data must be saved in IntanData folder
%INPUT: exp_name string which is the begining of raw data folder name. e.g.
%= 'c57-0010-m3'; for folder 'c57-0010-m3-t3-sf-red-out_ephys_190911_125159'
%OUTPUT: varibale containing cells with all folders containing that
%exp_name inside IntanData directory

exp_name_len = length(exp_name);

d = dir('E:\extarcellular\Andreanne\IntanData');
d = d(~ismember({d(:).name},{'.','..'})); %removes '.' and '..' folders

%initialize var 
exp_folders = cell(length(d),1);

index = 1;
for i=1:length(d)
    if strcmp(d(i).name(1:exp_name_len), exp_name)
        exp_folders{index} = d(i).name;
        index = index + 1;
    end
end

%remove extra empty cells
exp_folders = exp_folders(~cellfun('isempty',exp_folders));


