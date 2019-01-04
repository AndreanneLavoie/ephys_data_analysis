%import intan analog signals into mat lab

function v = importIntanAnalog(fileName, filedir)

    %filedir = 'C:\Users\liulabextracellular\Documents\intan test data\testAnalogInput_A16_ephys__181214_104251';
    FullfileName = fullfile(filedir, fileName);

    fileinfo = dir(FullfileName); 
    num_samples = fileinfo.bytes/2; % uint16 = 2 bytes 
    fid = fopen(FullfileName, 'r'); 
    v = fread(fid, num_samples, 'uint16'); 
    fclose(fid); 
    v = v * 0.000050354; % convert to volts (only valid if board mode == 0) 
end