function vs_timepoints = vs_time2(time_file_name, analog_file_name, filedir) %maybe instead take the filename and raw data dir path as input and call other functions inside
    %Convert analog trace to digital, based on thereshold of 0.3V
    %Return each time point in which the digital signal goes from 0 to 1;

    %import raw data into matlab, includes time and one analog input (vs timing, led timing, behviour, camera timing) 
    t = importIntanTime(time_file_name, filedir); %ADD DIR
    v = importIntanAnalog(analog_file_name, filedir); %ADD DIR

    v01 = NaN(1,length(v));

    %assigh either 0 or 1 to elements of v01 for every value in v 
    for i = 1:length(v)
        if v(i)> 0.3
            v01(i) = 1;
        else
            v01(i) = 0;
    end

    %vs_timepoints = for all elements in v01, if i in v10 != i+1 and i==0, save value of t[i+1];

    vs_timepoints = nan;
    
    %compare each neiboring pairs of vales in v01 to each other and find
    %time point corresponding to pattern '0 1'  
    for i = 1:(length(v01))
        if v01(i) ~= v01(i+1) & v01(i) == 0 
            if length(vs_timepoints) >= 1
                vs_timepoints = [vs_timepoints, t(i+1)]; %append if not empty 
            else %if timepoints is empty
               vs_timepoints = t(i+1); 
            end
        end
    end
end