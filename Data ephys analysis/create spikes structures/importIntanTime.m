function timestamps = importIntanTime(time_file_name, filedir, sampling_rate)
%this functions import timestamps during recording 

%directory where data is stored *SHOULD THIS BE AN INPUT OF THE FUNCTION?*
if nargin < 3
    sampling_rate = 30000;
end
    
%filedir = 'C:\Users\liulabextracellular\Documents\intan test data\testAnalogInput_A16_ephys__181214_104251';


timeFullFileName = fullfile(filedir, time_file_name); 

%extract header information to pull out sampling frequency 
%function will automatically dump all variables in workspace
%read_Intan_RHD2000_file

%pull time.dat into matlab:
fileinfo = dir(timeFullFileName); 
num_samples = fileinfo.bytes/4; % int32 = 4 bytes 
fid = fopen(timeFullFileName, 'r'); 
timestamps = fread(fid, num_samples, 'int32'); 
fclose(fid); 
timestamps = timestamps / sampling_rate; %frequency_parameters.amplifier_sample_rate; % sample rate from header file ; % this number comes from sample rate from 
%header file info.rhd and can be read in matlab using 
%read_Intan_RHD2000_file function; 
end

