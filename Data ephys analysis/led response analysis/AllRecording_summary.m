%% This function recordings in kilosort folder and create a summary
function AllRecording_summary(user,mode)
%%mode
%'all': updates all recordings
%'new': only run the analysis for new ones and updates the info to kilosort
%sumary.xls and summary.mat
%by default, the function updates new recordingds to summary.mat and kilosort summary.xls; 
%but if these two files do not exist in the folder, the function will run on all of the recordings

%%checklist outputs
%kilosort summary.xls: 
%(1)recording directory (2)	whether phy2
%has run and why not (3)whether create_spikes has run and why not
%(4)whether cluster analysis has run and why not (5)Intan channel number (6)notes
%kilosort_summary.mat: 
%kilosort_sumarylist: complete version of kilosort summary.xls with all componenets of error messages
%save in column 6
%kilosort_checklist: a logical version of summarylist for futher analysis
%(1)recording directory (2)	whether phy2 has run (3)whether create_spikes has run (4)whether cluster analysis has run

%%the checklist go throughs each folder
%1)try to update spikes construct if it exist
%2)check if tagging is used in the recordings and run cluster_analysis
%3)try to run create_spikes if spikes construct is absent: 
%(1) kilosort folder name is used to compute intan folder name. If multiple
%intan folder matches the name, it notes in kilosort checklist summary.xls,
%and create_spikes needs to run manually
%(2) only an overall psth will be plotted if create_spikes can run as
%cluster_info.tsv needed for cluster_analysis require manual saving in
%phy2. You can use "phy2_batch_save.m" to go through kilosort folder and
%open phy2 in the recordings without cluster_info.tsv.
%4)a "try-catch-continue" syntax is used to catch and note any possible
%errors occur when running create_spikes that I'm not aware of

%%to remove speicifc mice / recordings from running, add in exclude_mice /
%%exclude_recordings
if nargin<2
    mode = 'new';%only updates the new recordings by default
end

switch user
    case 'AL'
        kilosort_path='E:\extarcellular\Andreanne\KilosortData';%parent folder of all the kilosort data
        Intan_path='E:\extarcellular\Andreanne\IntanData';%parent folder of all the intan data
        electrode_config ='A2x16-Poly2';%the probe configuration of the recordings to be analyzed
    case 'HH'
        kilosort_path='E:\extarcellular\Hayley\KilosortData';%parent folder of all the kilosort data
        Intan_path='E:\extarcellular\Hayley\IntanData';%parent folder of all the intan data
        electrode_config ='A1x32-Poly3'; %the probe configuration of the recordings to be analyzed
end

%% grab all recordings in kilosort folder
cd(kilosort_path)%open parent folder

%if no summary.mat exist, the function will run on all of the recordings
if ~isfile('kilosort_summary.mat')
    mode = 'all';
end

%grab all NOT-DTN recordings: vgt, c57
all_folders = dir; %information for all folders in current directory
recording_list=string({all_folders([all_folders.isdir]).name})';%the names for all folders in the directory
recording_list = recording_list(startsWith(recording_list,["vgt-","c57-"]));%list folder names for all the recodings starts with 'vgt' or 'c57'

%these are the mice runs that need to be exclude from the list
exclude_mice= {'c57-0007-m0-'...
                };
%remove the recordings from exclude_mice
for k = 1:length(exclude_mice)
    recording_list(startsWith(recording_list,exclude_mice(k)))=[];
end

%these are the recordings that need to be exclude from the list
exclude_recordings={};         
%remove the recordings in exclude_recordsings from list
for k = 1:length(exclude_recordings)
    recording_list(startsWith(recording_list,exclude_recordings(k)))=[];
end

%% different update mode
switch mode
    % run on all the reocrdings in the folder
    case 'all'
        recording_list=recording_list;
    % update new recordings to list    
    case 'new'
         updatelist={};
         load('kilosort_summary.mat')
         for i = 1:length(recording_list)
             if ~strcmp(kilosort_checklist(:,1),recording_list(i))
                 updatelist(end+1)={recording_list(i)};
             end
         end
         recording_list=updatelist;
end
%% Check if phy2, create_spikes and led analyses have run
checklist=cell(length(recording_list),4);
%(1)recording directory (2) whether phy/phy2 has run (3) whether
%create_spikes has run (4) whether cluster analysis has run 
summarylist=cell(length(recording_list),7);
%(1)recording directory (2) whether phy/phy2 has run (3) whether
%create_spikes has run (4) whether cluster analysis has run (5)Intan
%channel number (6)notes (7) error message
for j = 1: length(recording_list)
    cd ([kilosort_path '\' recording_list{j}])
    checklist(j,1:4) = {recording_list{j},isfile('cluster_info.tsv'),isfile('spikes.mat'),isdir('cluster_analysis')};  
end
summarylist(:,1:4)=checklist(:,1:4);

%go through the checklist and open phy2 for any recordings that don't have cluster_info.tsv file
%manually save and close phy2 the loop will proceed to the next iteration
for m=1: size(checklist,1)
    if double(checklist{m,2})~=1
     [~,~]=system(['e:&&cd ' kilosort_path '\' recording_list{m} '&&activate phy2&&phy template-gui params.py']);
     %if the command is not recognized by cmd, check if anaconda has been added to the system variable
     %to add anaconda see E:\extarcellular\Andreanne\code --- Manual --- Installation --- 6.	To run anaconda from cmd
    end
end 

%update checklist
for n = 1: length(recording_list)
    cd ([kilosort_path '\' recording_list{n}])
    checklist(n,1:4) = {recording_list{n},isfile('cluster_info.tsv'),isfile('spikes.mat'),isdir('cluster_analysis')};  
end
summarylist(:,2:4)=checklist(:,2:4);

%% try run/udpate create_spikes & cluster_analysis if phy/ phy2 has run
for i = 1:length(recording_list) 
    %%try update create_spikes
    if double(checklist{i,2})==1 %phy2 has run
        if double(checklist{i,3})==1%if create_spikes has run
            %convert char to number with double format
            load([kilosort_path '\' recording_list{i} '\' 'spikes.mat'])
             if size(spikes.vs_params,2)>16%if the vs_params is the same as our current version
               if isfield(spikes,'raw_data_path')%if spikes structure has Intan path field
                  try spikes=create_spikes(spikes.raw_data_path);%update create_spikes
                  catch ME
                      summarylist{i,3}='error in updating creat_spikes';
                      summarylist{i,7}=ME;
                      % continue
                  end
               else%need to run create_spikes from intan folder
                   checklist{i,3}=0;
                   summarylist{i,3}=0;
               end
             else                    
               summarylist{i,3}='vs params is different from current version';
             end
       %%try re-run create_spikes from intan path
        else 
              cd(Intan_path)
              Intan_data=dir(['*' recording_list{i} '*']);%find corresponding Intan folder
             if ~(size(Intan_data,1)>1)%if only 1 Intan folder for that recording location
                cd([Intan_path '\' Intan_data.name])
                summarylist{i,5}=size(dir('amp-A*.dat'),1);%put in Intan channel numbers  
                try %trying running create_spikes and catch any possible error without stopping the loop 
                    spikes=create_spikes([Intan_path '\' Intan_data.name]);
                catch ME
                    summarylist{i,3}= ME.message;%display the error message in the second column to explain why create_spikes can't run
                    summarylist{i,7}= ME;
                    continue
                end
                checklist{i,3}=1;%update check list
                summarylist{i,3}=1;%update summary list 
            else
                summarylist{i,3}='multiple intan exist for the recording';
            end
         end
    else
       summarylist{i,2}='need to run phy2';
    end
end

%% save an overall psth image in each recording folder
for i=1:length(recording_list)
    cd([kilosort_path '\' recording_list{i}])%save psth in recording folder
    if isfile('spikes.mat')
        load([kilosort_path '\' recording_list{i} '\' 'spikes.mat'])
        bin_psth = 50;
        figure
        psth = psth_bl(spikes, bin_psth); %plot overall psth for recordings without cluster_info
        title('PSTH')
        saveas (gcf, 'overall PSTH.jpeg')
        saveas (gcf, 'overall PSTH.fig')
    end
end  

%% try run cluster_analysis
for i=1:length(recording_list) 
    tagging_flag=0; %1 tagging is used, 0 tagging is not used
    if double(checklist{i,3})== 1
       load([kilosort_path '\' recording_list{i} '\' 'spikes.mat'])
       if all(spikes.vs_params(:,17)) %if all trials have led
           if abs(spikes.vs_params(1,13)-spikes.vs_params(1,12))<spikes.vs_params(1,12)*0.2 % 20% difference for trailing off
              summarylist{i,4}='a single pulse of led is apply after vs stimulus';
           else %tagging is used in the recording, thus do the tagging analysis
               tagging_flag = 1;
              if sum(spikes.labels(:,2)==3)>10%if more than 10 unsorted clusters exist, defined as merging not done
                 summarylist{i,6}='have not done merging';
              elseif sum(spikes.labels(:,2)==3)>0 && sum(spikes.labels(:,2)==3)<10
                 summarylist{i,6}='unsorted cluster exist';
              end
              spikes.labels(spikes.labels(:,2)==3,2)=2;%plot all unassigned cluster as well
                  
               kilosort_name=checklist{i,1};
               hyphen=strfind(kilosort_name,'-');%find hyphen position
               tuning =  kilosort_name(hyphen(4)+1:hyphen(4)+2);%extract tuning information from text after the last underscore                              
               if tuning == 'or'
                  tuning = 'ori';
                  summarylist{i,4}=1;
              elseif tuning == 'sf'
                  summarylist{i,4}=1;
              elseif tuning == 'tf'
                   summarylist{i,4}='tuning is sf';
              else
                   summarylist{i,4}='tuning is not sf,tf or ori';
               end
           end
        elseif all(~spikes.vs_params(:,17))%if non of the trials have led
           summarylist{i,4}='led is not used in the recording';
        else %led is on from some of the trials
           if spikes.vs_params(1,9)== spikes.vs_params(1,1)
              summarylist{i,4}='led is applied with vs stimulus';
           end
       %there are currently no other possibilities
       end
       if isempty(summarylist{i,5})
          summarylist{i,5}=size(spikes.intan_info.amplifier_channels,2);%put in Intan channel numbers
       end
       try spikes=cluster_analysis(spikes,[],tuning,electrode_config, tagging_flag);
        catch ME
           summarylist{i,4}= ME.message;
           summarylist{i,7}=ME;
           % continue
        end
        checklist{i,4}=1;%update checklist, summarylist will be used to take down info
    end
    fprintf([num2str(length(recording_list)-i) ' recordings left\n'])
end       
%% save checklist to kilosort folder
switch mode
    case 'all'
        kilosort_summarylist=summarylist;
        kilosort_checklist=checklist;
    case 'new' %updates result to kilosort_summarylist and kilosort_checklist
        kilosort_summarylist={kilosort_summarylist;summarylist};
        kilosort_checklist={kilosort_checklist;checklist};
end
cd(kilosort_path)
%remove old files
if isfile('summary.mat')
    delete ('summary.mat')
end
if isfile('kilosort summary.xls')
    delete('kilosort summary.xls')
end
%save updated files
T = cell2table(kilosort_summarylist);
header={'recording','phy2','create_spikes','cluster_analysis','Intan_channel_number','notes','ME_messsage'};
T.Properties.VariableNames = header;%rename variables in table
writetable(T,'kilosort summary.xls')
save('kilosort_summary.mat','kilosort_checklist','kilosort_summarylist')
end