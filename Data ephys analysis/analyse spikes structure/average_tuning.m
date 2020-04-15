%average tuning curve normalization 



positive_flag = sf_tuning_fr > 0;

selected_units = sum(positive_flag, 2) >= 4.5;

max_sf_tuning_fr = max(sf_tuning_fr,[], 2);

max_sf_tuning_fr = repmat(max_sf_tuning_fr, 1, size(sf_tuning_fr,2));

norm_sf_tuning_fr = sf_tuning_fr ./ max_sf_tuning_fr;

ave_sf_tuning_fr = mean(norm_sf_tuning_fr(selected_units, :));

