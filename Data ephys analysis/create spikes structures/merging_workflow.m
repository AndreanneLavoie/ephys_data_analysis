% MERGING WORKFLOWs

%1 - mark noise and mua in phy

path = ''; 
spikes = create_spikes(path);
assigns = [];

% 2 - merge clusters in phy
 evaluate_merging_cluster(spikes, [assigns], 'sf', 1);
 %save changes in phy
 
%3 - update spike construct after merging in phy 
 [spikes.ref_period_violation, spikes.assigns] = checkRefractoryPeriodViolations_KiloSort_AL(spikes, [], 'all');
 
 
 %repeat step 2 and 3