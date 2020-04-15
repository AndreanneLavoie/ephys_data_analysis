%for every assign computer wether it photo-tagged

spikes.cluster_char = spikes.labels;
spikes.est_width = spikes.labels(:,1);

for c=1:length(spikes.labels(:,1))
    
   assign = spikes.labels(c,1);
   
   %%determine if each assign has visual driven activity
   vis_flag = is_visual_driven(spikes, assign);
   
   %%determinr if each assing has led driven activity
   led_flag = shuffle_led_spiketimes_10000(spikes, assign);
   spikes.cluster_char(c, 2) = led_flag;

   %%assing var to spike field
   %spikes.cluster_char(c, 2:3) = [vis_flag led_flag];

   %previous used to test
   spikes.cluster_char(c, 3) = vis_flag;
   
   %calculate the width of hist (estimation); 
   spikes.est_width(c) = calc_jitters(spikes, assign);
end