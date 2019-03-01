function vs_timepoints = vs_time2(time_file_name, analog_file_name, filedir, threshold, digfilter_winT) %maybe instead take the filename and raw data dir path as input and call other functions inside
    %Convert analog trace to digital, based on thereshold of 0.3V
    %Return each time point in which the digital signal goes from 0 to 1;

    %import raw data into matlab, includes time and one analog input (vs timing, led timing, behviour, camera timing) 
    t = importIntanTime(time_file_name, filedir); %ADD DIR
    deltat=mean(diff(t)); %in sec
    digfilter_winsize=round(digfilter_winT/deltat);
    v = importIntanAnalog(analog_file_name, filedir); %ADD DIR
%     v01 = NaN(1,length(v));

    %assigh either 0 or 1 to elements of v01 for every value in v 
%     for i = 1:length(v)
%         if v(i)> threshold
%             v01(i) = 1;
%         else
%             v01(i) = 0;
%         end
%     end
    dig_threshold=v>threshold;
    %vs_timepoints = for all elements in v01, if i in v10 != i+1 and i==0, save value of t[i+1];
    digtemp=[];
    for i=1:digfilter_winsize
        digtemp=[digtemp circshift(dig_threshold, i*-1)];
    end
    dig_threshold=all(digtemp, 2);
    der1st=diff(dig_threshold);
    %compare each neiboring pairs of vales in v01 to each other and find
    %time point corresponding to pattern '0 1'  
%     for i = 1:(length(v01)-1)
%         if v01(i) ~= v01(i+1) & v01(i) == 0 
%             if length(vs_timepoints) >= 1
%                 vs_timepoints = [vs_timepoints, t(i+1)]; %append if not empty 
%             else %if timepoints is empty
%                vs_timepoints = t(i+1); 
%             end
%         end
%     end
der1st=[0; der1st];

vs_start=t((der1st>0.5));
vs_end =t((der1st<-0.5));

vs_timepoints=[vs_start vs_end vs_end-vs_start];
end