%%check if clusters have actual led response instread of high spontaneouse firing rate
function led_response_clusters=is_led_response(spikes, assigns, tuning)
% spikes: spikes construct for the recording
% assigns: cluster IDs, put "[]" to analyze all good clusters
% tuning: the type of tuning: 'sf', 'tf', 'ori'

%% modify inputs
if isempty(assigns) %if no unit assign is specified, check all good units
   assigns=spikes.labels(spikes.labels(:,2)==2,1); 
end

%get led value for filtspikes function
led=unique(spikes.led);%get led value for led-on trials
led(isnan(led))=[];
led(led<1e-8)=[];

%initiate output
led_response_clusters=[];

%% determine time stcture per trial 
led_start=(spikes.vs_params(1,11)-spikes.vs_params(1,9))/1000;%led start time relative to trigger (in s)
led_dur=spikes.vs_params(1,13)/1000;%led pulse duration
led_num=spikes.vs_params(1,16);%number of led pulses per trial

switch tuning
    case 'sf'
        vs_end= (spikes.vs_params(1,7)-spikes.vs_params(1,9))/1000;%vs ends time relative to trigger (in s)
    case 'tf'
    case 'ori'%this part my change as Andreanne changes the ehpys code
        vs_end=(spikes.vs_params(1,6)-spikes.vs_params(1,9))/1000;%vs ends time relative to trigger (in s)
end

%% calculate firing rates
for i=1:length(assigns)
    
    filt_spikes = filtspikes(spikes, 0, 'assigns', assigns(i),'led',led);
    if any(abs(filt_spikes.labels(:,1)-assigns(i))<1e-8) %chek if the cluster exist
        if filt_spikes.labels(filt_spikes.labels(:,1)==assigns(i),2)== 2 %check if the cluster is good

            trials=unique(filt_spikes.trials);
            firing_rates=zeros(length(trials),3); 
            %the second column for spontaneous firng rate, the thrid column for averaged led firing rates
            firing_rates(:,1)=trials;%first column is trail number
            %sort spikes into different trials
            for j=1:length(trials)
                filt_spikes_trials=filtspikes(filt_spikes,0,'trials', trials(j),'led',led);
                
                %calculate the spontaneous firing rate using the time after vs stimulus ends and before led on
                num_spikes_sp=sum(filt_spikes_trials.spiketimes<led_start & filt_spikes_trials.spiketimes>vs_end);
                firing_rate_sp=num_spikes_sp/(led_start-vs_end);
                
                %calculate the average firing rate during led stimuli
                num_spikes_led=sum(~isnan(filt_spikes_trials.led_spiketimes));
                firing_rate_led=num_spikes_led/led_dur/led_num;
                
                firing_rates(j,2:3)=[firing_rate_sp firing_rate_led];
            end
            
            %paired t test
            if ttest(firing_rates(:,2), firing_rates(:,3))
               fprintf (['cluster' num2str(assigns(i)) ' has led response\n'])
               led_response_clusters=[led_response_clusters;assigns(i)];%save the cluster ID for led responsive ones
            end
            
            %update in spikes
            spikes.led_response_cluster=led_response_clusters;
            
       else
            fprintf (['cluster' num2str(assigns(i)) 'is not a good unit\n'])
       end
        
    else
        fprintf (['cluster' num2str(assigns(i)) 'does not exit\n'])
        
   %save the spikes construct to its folder
   save([spikes.kilosort_path 'spikes.mat'], 'spikes');
end
end
