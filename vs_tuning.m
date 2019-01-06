function [tuning_para, tuning_val] = vs_tuning(fileName, filedir)
%This function will take .npy file and extract the value of changing tuning
%variable;
% NOTE: requires NPY to matlab package which can be obtained from:
% https://github.com/kwikteam/npy-matlab

FullfileName = fullfile(filedir, fileName);

%import all 11 columns from python vs params variable 
vs = readNPY(FullfileName);

%find the column that was used as tuning paramerter; ie the one that is
%actally changing within the column;

%Will compare the  first and second value to each other for every column
%and return the values that were used to tune and the parameter it
%corresponds to;
for i = 1:11
    if vs(1,i) ~= vs(2,i)
        tuning_val = vs(:,i);
        
        if i == 1
            tuning_para = 'spat';
        elseif i == 2
            tuning_para = 'temp';
        elseif i == 3
            tuning_para = 't_bef';
        elseif i == 4
            tuning_para = 't_dur';
        elseif i == 5
            tuning_para = 't_aft';
        elseif i == 6
            tuning_para = 'cont';
        elseif i == 7
            tuning_para = 'nan';
        elseif i == 8
            tuning_para = 'nan';
        elseif i == 9
            tuning_para = 'nan';
        elseif i == 10
            tuning_para = 'ori';
        elseif i == 11
            tuning_para = 'led';
        end
    end
end 
end

