%exmaine changing alpha val on num units that are led_driven

spikes.is_led_driven = spikes.labels(spikes.labels(:,2) == 2);
num_units = length(spikes.is_led_driven);

alpha_value = [0.05 0.005 0.0005 0.00005 0.000005]; %%%%
num_alpha_vals = length(alpha_value);
is_led_driven_temp = nan(num_units, num_alpha_vals);

for k = 1:length(alpha_value)
   spikes = is_led_driven(spikes, alpha_value(k));
   
   is_led_driven_temp(:, k) = spikes.is_led_driven(:, 2);
end

% plot
figure;
sum_led_driven_units = sum(is_led_driven_temp);
plot(sum_led_driven_units);

figure;
bar(log(alpha_value), sum_led_driven_units)

xlabel(['log(alpha) ' num2str(flip(alpha_value))])
ylabel('num led-driven units')
