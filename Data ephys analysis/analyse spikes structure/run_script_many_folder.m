%%find folder name and run script

ParentFolder='E:\extarcellular\Andreanne\IntanData\';
exp_name = 'c57-0014-m3';
Fig_path = 'E:\extarcellular\Andreanne\data_analysis\NOTDTN recording practice\';

exp_name_len = length(exp_name);
AllFile=dir(fullfile(ParentFolder,'.'));
FolderOnly=AllFile([AllFile.isdir]);

FolderNames = {};

for k=1:length(FolderOnly)
    
    if length(FolderOnly(k).name) > exp_name_len - 1 
        
        if all(FolderOnly(k).name(1:exp_name_len) == exp_name)
            
           FolderNames{end+1} = FolderOnly(k).name
%            spikes = create_spikes(fullfile(ParentFolder,FolderNames{k}));
%            mkdir(Fig_path, FolderNames{k});
%            psth_depth;
%            savefig(fullfile(Fig_path, FolderNames{k}), 'psth_depth_sf.fig');
        end
        
    end
  
end

%remove bad trial
% FolderNames{1,1} = FolderNames{1,15:end};
% j = 13;

%for j=1:length(FolderNames)
    
%     spikes = create_spikes(fullfile(ParentFolder,FolderNames{1,j}));
%     mkdir(Fig_path, FolderNames{1,j});
%     psth_depth;
%     savefig([fullfile(Fig_path, FolderNames{1,j}) '\psth_depth_sf.fig']);
%     close all;
  
%end