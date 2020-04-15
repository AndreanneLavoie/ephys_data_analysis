function spikes=cluster_analysis(spikes,assigns,tuning, electrode_config, tagging, T_led_disp,T_nat, binwidth_hist)
% spikes: spikes construct for the recording
% assigns: cluster IDs, put "[]" to analyze all good clusters
% tuning: 'sf', 'tf', 'ori'
% electrode_config: the configuration of recording electrode 
%                   'A1x16', 'A1x16-Poly2', 'A1x32-Poly2', 'A1x32-Edge', 
%                   'A2x32', 'A2x32-Poly5', 'A1x64-Poly2', 'A1x32', 
%                   'A1x32-Poly3', 'A2x16-Poly2'
% tagging: whether tagging is used, 1 if tagging is used, 0 if not
% T_led_disp: the displayed time window after led event for led response average
% T_nat: time window for natural spike average


if nargin<5
    tagging=1;
end

if nargin<6
    T_led_disp=25;
end

if nargin<7
   T_nat=16;
end

if nargin<8
    binwidth_hist=0.0002;
end

%% variables
sampling_rate=spikes.intan_info.frequency_parameters.amplifier_sample_rate;
T_led_trunctrace = T_led_disp/1000; %the time window (sec) of truncated trace
T_nat_half_trunctrace = T_nat/2/1000; %the time wrindow (sec) before / after naturally occurred spike peaks
numsamp_led_trunctrace= round(sampling_rate*T_led_trunctrace); %the number of samples display in led plot
led_period_temp=diff(spikes.led_abstime(:,1));
temp=abs(led_period_temp - spikes.vs_params(1,13)/1000)>spikes.vs_params(1,13)/1000*0.1;%remove the time differece between leds of different trials
led_period_temp(temp) = [];
sample_per_led_period = int32(mean(led_period_temp)*sampling_rate);%the number of samples in each led period
numsamp_nat_half_trunctrace= round(sampling_rate*T_nat_half_trunctrace); %the number of samples in half of the natural response truncated trace
bin_psth = 50;
T_hist = [0 25]; %20; % time window for calculating jitter, in ms
digfilter_winsize=2;%digfilter window for waveform corelation analysis
% 2 = tolerate one false value

%  save latency waveform enery and cross-correlation in spike construct
% spikes.latency   1st column: clusterID; 2nd column: jitter (ms); 3nd columns: mean latency (ms)
num_good_unit=length(spikes.labels(spikes.labels(:,2)==2));
latency = nan(num_good_unit,3);  
latency(:,1) = spikes.labels(spikes.labels(:,2)==2,1); 
% spikes.waveform_energy    1st column: clusterID; 2nd column: waveform energy of the averaged LED response; 3rd column: waveform energy of the average natural spikes.
waveform_energy =nan(num_good_unit,3);
waveform_energy(:,1) = spikes.labels(spikes.labels(:,2)==2,1); 
% spikes.cross_correlation    1st column: clusterID; 2nd column: cross-correlation between averaged LED response and average natural spikes;
cross_correlation = nan(num_good_unit,2);
cross_correlation(:,1) = spikes.labels(spikes.labels(:,2)==2,1); 
peak_dis = nan(num_good_unit,2);%(1)cluster id (2)peak distance
peak_dis(:,1) = spikes.labels(spikes.labels(:,2)==2,1);

%% load chanel map
switch electrode_config
    case 'A1x16'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA1x16.mat');
    case 'A1x16-Poly2'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA1x16Poly2.mat'); 
    case 'A1x32-Poly2'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA1x32Poly2.mat'); 
    case 'A1x32-Edge'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA1x32Edge.mat'); 
    case 'A2x32'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA2x32.mat'); 
    case 'A2x32-Poly5'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA2x32-Poly5.mat'); 
    case 'A1x64-Poly2'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA1x64Poly2.mat'); 
    case 'A1x32'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA1x32Classic.mat'); 
    case 'A1x32-Poly3'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA1x32Poly3.mat'); 
    case 'A2x16-Poly2'
        ChanelMap = load('E:\extarcellular\Andreanne\chanMaps\chanMapA2x16-Poly2.mat'); 
end
    
    
%% Analyzing led and natural response for by cluster
if isempty(assigns) %if no unit assign is specified, plot all good units and un-assigned units
   assigns=spikes.labels(spikes.labels(:,2)==2,1); 
end

%create saving directory, delete all previous images for update
cd (spikes.kilosort_path)
if isdir('cluster_analysis')
   cd('cluster_analysis')
   delete('*.jpeg')
else
   mkdir ('cluster_analysis')
end

for i=1:length(assigns)
    %% cluster check
     if any(abs(spikes.labels(:,1)-assigns(i))<1e-8) %chek if the cluster exist
        if spikes.labels(spikes.labels(:,1)==assigns(i),2)== 2 %check if the cluster is good
            %% read in Intan data
            %get intan channel info
            filt_spikes_unit_nat = filtspikes(spikes, 0, 'assigns', assigns(i));
            %the channel number in phy2
            unit_channel = filt_spikes_unit_nat.clusterinfo.channel(filt_spikes_unit_nat.clusterinfo.id==assigns(i));
            %convert to actual Intan channel number
            unit_channel = ChanelMap.Intan_chan_num(ChanelMap.chanMap(ChanelMap.chanMap==(unit_channel+1)));
            if unit_channel<10
               unit_channel = [num2str(0) num2str(unit_channel)];
            else
               unit_channel = num2str(unit_channel);
            end
            
            FullfileName = [spikes.raw_data_path '\' 'amp-A-0' unit_channel '.dat'];
            fileinfo = dir(FullfileName); 
            num_samples = fileinfo.bytes/2; % int16 = 2 bytes 
            fid = fopen(FullfileName, 'r'); 
            Chan_data = fread(fid, num_samples, 'int16'); 
            fclose(fid); 
            Chan_data = Chan_data * 0.195; %in uV
            
            %filter data using band-pass
            Chan_data = filtdata(Chan_data,sampling_rate,[],'band',[500 10000],[300 12000]); % Band-pass 0.5-10 kHz
           
            %% analyze spontaneous response
            nat_abs_spiketimes=filt_spikes_unit_nat.abs_spiketimes(filt_spikes_unit_nat.spiketimes<(filt_spikes_unit_nat.vs_params(1,11)/1000));%vs params in ms
            % remove spikes before first vs and after last vs end (last vs time + vs period)
            nat_abs_spiketimes = nat_abs_spiketimes(nat_abs_spiketimes>filt_spikes_unit_nat.vstiming(1,1) & nat_abs_spiketimes<(filt_spikes_unit_nat.vstiming(end,1)+mean(filt_spikes_unit_nat.vstiming(:,3))));
            %aligning naturally occurred spikes by peaks.
            reponse_nat=zeros(length(nat_abs_spiketimes), (2*numsamp_nat_half_trunctrace+1)); % numsamp_nat_half_trunctrace before the spike and numsamp_nat_half_trunctrace after the spike
            for m=1:length(nat_abs_spiketimes)
                sample_nat_spike_ind=round(nat_abs_spiketimes(m)*sampling_rate);%get peak sample index
                reponse_nat(m,:)=Chan_data((sample_nat_spike_ind-numsamp_nat_half_trunctrace):(sample_nat_spike_ind+numsamp_nat_half_trunctrace))';
                %the max abs value always at the middle of each line
            end
            average_reponse_nat=mean(reponse_nat,1); %average recording for each relative time point (column by column), in uV                
            
            % calculate the time distance between the largest positive and largest negative peak
            [peak_positive,nat_positive_ind]=max(average_reponse_nat);
            [peak_negative,nat_negative_ind]=min(average_reponse_nat);
            %(1)cluster id (2)peak distance
            peak_dis(peak_dis(:,1)==assigns(i),2)=abs(nat_positive_ind-nat_negative_ind)/sampling_rate*1000; %in ms
            
            %% calculate led response per cluster
            if tagging ==1
                led=unique(spikes.led);%get led value for led-on trials
                led(isnan(led))=[];
                led(led<1e-8)=[];

                filt_spikes_unit_led = filtspikes(spikes, 0, 'assigns', assigns(i),'led',led);

                % find the corresponding led associated with each first spikes
                uniledind=unique(filt_spikes_unit_led.ledind);
                uniledind(abs(uniledind)<1e-8)=[];
                firstspikeIND=[];
                for j=1:length(uniledind)
                   firstspikeIND=[firstspikeIND find(abs(filt_spikes_unit_led.ledind-uniledind(j))<1e-8,1)];
                end
                %abs led time for the first spikes of the cluster
                unit_led_abstime=filt_spikes_unit_led.led_abstime(filt_spikes_unit_led.ledind(firstspikeIND),1);%in sec
                
                %aligning to LED event, each row is one led event
                reponse_led=zeros(length(unit_led_abstime), sample_per_led_period);
                for n=1:length(unit_led_abstime)
                    sample_led_start_ind=int32(round(unit_led_abstime(n)*sampling_rate));%get start sample index
                    reponse_led(n,:)=Chan_data(sample_led_start_ind:(sample_led_start_ind+sample_per_led_period-1))';
                end
                average_response_led=mean(reponse_led,1); %average recording for each relative time point (column by column), in uV
                shifted_average_response_led=circshift(average_response_led,10/1000*sampling_rate,2);

                % align the peaks of the led response and the natural spikes average
                [temp,nat_max_ind]=max(abs(average_reponse_nat));
                if average_reponse_nat(nat_max_ind)>0
                   [temp,led_peak_ind]=max(shifted_average_response_led);
                else
                   [temp,led_peak_ind]=min(shifted_average_response_led);
                end 
                peak_ind_diff=nat_max_ind-led_peak_ind;
                shifted_average_response_led=circshift(shifted_average_response_led,peak_ind_diff);%align shifted_average_response_led and average_reponse_nat
            
                %% calculate jitter, waveform energy and cross correlation

                % jitter: std of the first spike for each led
                first_led_spiketime=filt_spikes_unit_led.led_spiketimes(firstspikeIND)*1000; %in ms
                first_led_spiketime_trunct=first_led_spiketime(first_led_spiketime<T_hist(2) & first_led_spiketime>T_hist(1));
                latency(latency(:,1)==assigns(i),2:3) = [std(first_led_spiketime_trunct) mean(first_led_spiketime_trunct)];
                    
                % 1st column: clusterID; 2nd column: jitter (ms); 3nd columns: mean latency (ms)

                if size(reponse_nat,1)>1
                    %find start and end index of natural trace for coeff calculation
                    [h,p]=ttest(reponse_nat);
                    digtemp=zeros(length(h), digfilter_winsize);
                    digtemp(:,1)=h;
                    for k=1:(digfilter_winsize-1)
                        digtemp(:,k+1)=circshift(h, k*-1);
                    end
                    h_filt=any(digtemp, 2); 
                    left_ind=find(~h_filt(1:nat_max_ind),1,'last')+1;
                    right_ind=(nat_max_ind-1)+find(~h_filt(nat_max_ind:end),1,'first')-1;

                    %waveform energy is the integral of v^2 over time ((uV)^2*ms)
                    waveform_energy_led = sum((shifted_average_response_led(left_ind:right_ind).^2)*(1/sampling_rate*1000));
                    waveform_energy_nat = sum((average_reponse_nat(left_ind:right_ind).^2)*(1/sampling_rate*1000));
                    waveform_energy (waveform_energy(:,1)==assigns(i),2:3) =  [waveform_energy_led waveform_energy_nat];
                    % 1st column: clusterID; 2nd column: waveform energy of the averaged LED response; 3rd column: waveform energy of the average natural spikes.

                    cross_cor = corrcoef(shifted_average_response_led(left_ind:right_ind), average_reponse_nat(left_ind:right_ind));
                    cross_correlation (cross_correlation(:,1)==assigns(i),2) = cross_cor(1,2);
                    % 1st column: clusterID; 2nd column: cross-correlation between averaged LED response and average natural spikes ; 3rd cloumn: cross-correlation lag
                else %no waveform_energy or cross_corr calculated if only 1 nature response spike exist
                    waveform_energy (waveform_energy(:,1)==assigns(i),2:3) =  [nan nan];
                    cross_correlation (cross_correlation(:,1)==assigns(i),2) = nan;
                end
            end
            %% other features of the cluster
            % tuning curve
            [stim_cond_list, tuning_fr_total, tuning_fr_spon] = compute_tuning(spikes, assigns(i), tuning, 0,1);
            tuning_fr_evoked = tuning_fr_total - tuning_fr_spon;
            
            % PSTH

      %% plotting the led response and natrual spike averages per cluster
            %plot in ms
            T=(1:numsamp_led_trunctrace)/sampling_rate*1000;
            figure;
            %averaged led response, align by the led event
            subplot(3,3,1);
            if tagging ==1
                if any(average_response_led) %if there is averaged led response
                    plot(T, average_response_led(1:numsamp_led_trunctrace)); %plot average trace
                else
                    text(1,0.5,'no averaged led response')
                end
                xlim([0 T_led_trunctrace*1000]);
                xlabel('ms');
                ylabel('mV');
                title(['led response' ' (' num2str(size(reponse_led, 1)) ')']);   
            else
                h1=subplot(3,3,1);
                set(h1,'Visible','off')
            end


            
            xlabel('ms');
            %averaged led reponse with raw trace
            subplot(3,3,2)
            if tagging ==1
                if ~isempty(reponse_led) %if there is led response
                    plot(T, reponse_led(1:numsamp_led_trunctrace));
                    hold on
                    plot(T, average_response_led(1:numsamp_led_trunctrace),'linewidth', 5);
                else
                    text(1,0.5,'no led response')
                end
                xlim([0 T_led_trunctrace*1000]);
                ylabel('mV');
                title('led spikes raw trace');
            else
                h2=subplot(3,3,2);
                set(h2,'Visible','off')
            end

            
            %plot first spike histogram
            subplot(3,3,3); 
            if tagging ==1
                if ~isempty(reponse_led) %if there is led response and 
                    binsize=round((max(first_led_spiketime)-min(first_led_spiketime))/binwidth_hist);
                    hist(first_led_spiketime,binsize); %histogram in ms
                else
                    text(1,0.5,'no led response')
                end
                xlim([0 T_led_trunctrace*1000]);
                xlabel('ms');
                ylabel('nSpikes');
                title(['led first spikes (jitter' num2str(latency(latency(:,1)==assigns(i),2)) ')']);
            else
                h3=subplot(3,3,3);
                set(h3,'Visible','off')
            end
            
            
            %plot all spike histogram
            subplot(3,3,6)
            if tagging == 1
                if ~isempty(reponse_led) %if there is led response and 
                    binsize=round((max(first_led_spiketime)-min(first_led_spiketime))/binwidth_hist);
                    if binsize ~= 0 
                        [~,~,~] = make_first_led_spiketimes(spikes, assigns(i), binsize, [0 T_led_disp], 0, 1)
                        hist(first_led_spiketime,binsize); %histogram in ms
                    else
                        text(1,0.5,'max(first_led_spiketime)=min(first_led_spiketime)')
                    end
                else
                    text(1,0.5,'no led response')
                end
                xlim([0 T_led_trunctrace*1000]);
                xlabel('ms');
                ylabel('nSpikes');
                title(['led all spikes (jitter' num2str(latency(latency(:,1)==assigns(i),2)) ')']);
            else
                h6=subplot(3,3,6);
                set(h6,'Visible','off')
            end
            
            %average spontanouse spikes, align by the peak
            subplot(3,3,4)
            t=(1:(numsamp_nat_half_trunctrace*2+1))/sampling_rate*1000;
            plot(t, average_reponse_nat); %plot average trace
            hold on
            plot(nat_positive_ind/sampling_rate*1000,peak_positive,'r*')
            plot(nat_negative_ind/sampling_rate*1000,peak_negative,'r*')
            xlim ([0 t(end)])
            xlabel('ms');
            ylabel('mV');
            title(['natural spikes average' ' (' num2str(size(reponse_nat, 1)) ')']);
            if isempty(reponse_nat)
               text(1,0.5,'no visual driven activity')
            end
            
            %compare average led response and spontaneouse spikes
            %align by either the largest positive or the largest negative
            %peak, depending on the max abs amplitude in average natural
            %response
            subplot(3,3,5)
            if tagging == 1
                plot(t, average_reponse_nat); %natural spikes average plot in default color
                hold on
                plot(t, shifted_average_response_led(1:length(average_reponse_nat)));
                xlabel('ms');
                ylabel('mV');
                title(['compare (' num2str(cross_cor(1,2)) ')']);
            else
                h4=subplot(3,3,5);
                set(h4,'Visible','off')
            end
              
            
            %tunning curve of the cluster
            subplot(3,3,7)
            switch tuning
            case 'sf'
                plot(stim_cond_list, tuning_fr_evoked);

            case 'tf'
                plot(stim_cond_list, tuning_fr_evoked);

            case 'ori'
                if any(tuning_fr_evoked<0)%if there is a negative value, force the smallest value to be 0
                    tuning_fr_evoked = tuning_fr_evoked-min(tuning_fr_evoked);
                    fprintf('tuning_fr_evoked contains negative value, force minimum to be zero\n');
                end
                theta = stim_cond_list/180*pi;
                theta = [theta; theta(1)];
                tuning_fr_evoked = [tuning_fr_evoked; tuning_fr_evoked(1)];
                polar(theta, tuning_fr_evoked);
            end
            title(['tuning: ' tuning]);
            
            %PSTH of the cluster
            subplot(3,3,8)
            filt_unit_psth = psth_bl(filt_spikes_unit_nat, bin_psth);
            title('PSTH')
            
            set(gcf, 'Position',  [100, 0, 1560, 960]) %set display parameter, 
            %[100 0] is the x,y coordinate of the lower right corner, [1560 960] is the length and height of the image
            suptitle(['cluster ' num2str(assigns(i)) ' (' num2str(sum(spikes.assigns==assigns(i))) ')'])
            
            %save imageg
            cd ([spikes.kilosort_path 'cluster_analysis' '\'])
            saveas(gcf, ['cluster ' num2str(assigns(i)) '.jpeg']);        
            saveas(gcf, ['cluster ' num2str(assigns(i)) '.fig']);
        else
            fprintf (['cluster' num2str(assigns(i)) 'is not a good unit\n'])
        end
    else
        fprintf (['cluster' num2str(assigns(i)) 'does not exit\n'])
    end
    

end
if tagging ==1
    spikes.cross_correlation = cross_correlation ;
    spikes.latency = latency;
    spikes.waveform_energy=waveform_energy;
end
spikes.peak_dis= peak_dis;
%% plot parameters
figure
scatter (latency(:,2), latency(:,3))
xlabel('jitter (ms)');
ylabel('mean latency (ms)');
saveas (gcf, 'jitter_mean latency.jpeg')
saveas (gcf, 'jitter_mean latency.fig')

figure
scatter(waveform_energy(:,2), cross_correlation(:,2))
xlabel('led response waveform energy ((uV)^2*ms)');
ylabel('waveform cross-correlation');
saveas (gcf, 'waveform energy_cross-correlation .jpeg')
saveas (gcf, 'waveform energy_cross-correlation .fig')

close all % close all figure window

save([spikes.kilosort_path 'spikes.mat'], 'spikes');%save the spikes construct to folder
end