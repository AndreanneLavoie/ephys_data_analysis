%save single unit information for temporary (presentation) purposes

units = [0 2 19 18 12 18 38 35 47 16 10 216 1 27 29 86 20];
cd (spikes.kilosort_path)
%can add make dir function to create special folder to save files in. Note
%will have to re cd into said folder once created.

for i = 1:length(units)
   assign = units(i);
   
   figure;
   compute_tuning(spikes, assign, 'ori', 1, 1);
   savefig(['unit' num2str(assign) '_tuning']); 
   
   spikes_unit = filtspikes(spikes, 0, 'assigns', assign);
   figure;
   psth_bl(spikes_unit, 10) %10 = bin size in ms
   savefig(['unit' num2str(assign) '_psth_bin_10ms']);
   
   figure; 
   [length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes, assign, 0.1, [0 25], 1, 1);
   savefig(['unit' num2str(assign) '_first_led_hist']);
   
end