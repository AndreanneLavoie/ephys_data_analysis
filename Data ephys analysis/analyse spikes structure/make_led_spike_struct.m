function filt_spikes = make_led_spike_struct(spikes, assign)
%INPUT
%spike structure (requires: ) 
%assign (cluster) that want to extract spiketimes from [1x1]

percT_thres = 0.1;% 10% of the led period difference
batchsize=50000;

%filter spiketimes by assign
if ~isempty(assign)
    temp=unique(spikes.led);%get led conditions
    temp(isnan(temp))=[];%remove nan
    temp(abs(temp)<1e-8)=[];%remove 0
    filt_spikes = filtspikes(spikes, 0, 'assigns', assign,'led', temp);
else
    filt_spikes = spikes;
end

%create extended matrix to allow comparision between all possible led
%times and all possible abs spiketimes
length_led = length(filt_spikes.led_abstime(:,1));


%calculate LED period
led_period_temp = filt_spikes.led_abstime(2:(end),1) - filt_spikes.led_abstime(1:(end-1),1);
%remove led period which is very different from standard value (e.g.
%between trails and noise)
led_period_temp(abs(led_period_temp - filt_spikes.vs_params(1,13)/1000)>filt_spikes.vs_params(1,13)/1000*percT_thres) = [];
led_period = mean(led_period_temp);

expttime=[filt_spikes.led_abstime(1,1) filt_spikes.led_abstime(end,1)+led_period];
filt_spikes.ledind=nan(size(filt_spikes.abs_spiketimes)); %the accumulative led number for each individual spike within recording; 
filt_spikes.ledind(filt_spikes.abs_spiketimes < expttime(1) | filt_spikes.abs_spiketimes > expttime(2))= 0; %0 = spikes happen before / after experiment

abs_spiketimes_durexp = filt_spikes.abs_spiketimes(isnan(filt_spikes.ledind)); %extract only spikes that happen within exp
numofspikes_durexp=length(abs_spiketimes_durexp);
numofcycle=floor(numofspikes_durexp/batchsize);
remainder=rem(numofspikes_durexp, batchsize);
ledindtemp=zeros(1, numofspikes_durexp);
for i=1:numofcycle
    temp_abs_spiketimes=abs_spiketimes_durexp((1+(i-1)*batchsize):(i*batchsize));
    abs_spikestimes_matrix = repmat(temp_abs_spiketimes, length_led, 1);%sec
    abs_ledtimes_matrix = repmat(filt_spikes.led_abstime(:,1), 1, batchsize); %sec

    %calculate rel spiketime to ref led
    diff_matrix = abs_ledtimes_matrix-abs_spikestimes_matrix;

    led_bef_spikestemp = diff_matrix < 0; %determin all triggers before each spikes
    led_bef_spikestemp(end+1,:) = 0; %add last row = 0 to make sure there is a 1 0 interface at the end;
    led_ind_bef_spikes = diff(led_bef_spikestemp, [], 1); %find interface (1 change to 0) by taking derivitive to find 
    [r c] = find(led_ind_bef_spikes < 0);  %extract coordinates where -1 located in matrix
    ledindtemp((1+(i-1)*batchsize):(i*batchsize))=r;
end
%process remainder
temp_abs_spiketimes=abs_spiketimes_durexp((1+numofcycle*batchsize):end);
abs_spikestimes_matrix = repmat(temp_abs_spiketimes, length_led, 1);%sec
abs_ledtimes_matrix = repmat(filt_spikes.led_abstime(:,1), 1, remainder); %sec

%calculate rel spiketime to ref led
diff_matrix = abs_ledtimes_matrix-abs_spikestimes_matrix;

led_bef_spikestemp = diff_matrix < 0; %determin all led triggers before each spikes
led_bef_spikestemp(end+1,:) = 0; %add last row = 0 to make sure there is a 1 0 interface at the end;
led_ind_bef_spikes = diff(led_bef_spikestemp, [], 1); %find interface (1 change to 0) by taking derivitive to find 
[r c] = find(led_ind_bef_spikes < 0);  %extract coordinates where -1 located in matrix
ledindtemp((1+numofcycle*batchsize):end)=r;

filt_spikes.ledind(isnan(filt_spikes.ledind)) = ledindtemp; %assign all exp spiketimes a trial ID number corresponding to order of stimuli presented

filt_spikes.led_spiketimes = filt_spikes.abs_spiketimes; 
filt_spikes.led_spiketimes(abs(filt_spikes.ledind-0)<1e-8)= nan;
filt_spikes.led_spiketimes(1,filt_spikes.ledind > 0) = filt_spikes.abs_spiketimes(1, filt_spikes.ledind > 0) - filt_spikes.led_abstime(filt_spikes.ledind(1, filt_spikes.ledind > 0),1)';
temp=filt_spikes.led_spiketimes>led_period;
filt_spikes.led_spiketimes(temp)= nan;
filt_spikes.ledind(temp)=0;

end