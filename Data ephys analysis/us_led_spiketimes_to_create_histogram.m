temp.led_spiketimes;
uniledind=unique(temp.ledind);
uniledind(abs(uniledind)<1e-8)=[];
firstspikeIND=[];
for i=1:length(uniledind)
   firstspikeIND=[firstspikeIND find(abs(temp.ledind-uniledind(i))<1e-8,1)];
end
timetemp=temp.led_spiketimes(firstspikeIND);
m=round((max(timetemp)-min(timetemp))*1000);
hist(timetemp,m)

%get led abs time for filtered spikes
temp.led_abstime(temp.ledind(firstspikeIND),1);