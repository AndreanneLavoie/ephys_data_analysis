%code will plot PSTH and led hist for led of different strenghts (has to be within one
%recording session)

assign = spikes.labels(spikes.labels(:,2) ~= 0);
num_trials_cond = 100;
colour = [[1 0 0];[0 1 0];[0 0 1]; [0 0 0]]; %low = red; med = green; high = blue
%good_spikes = filtspikes(spikes, 0, 'assigns', assigns);

%% plot psth at different led powers
figure;
hold on

[led_cond, index, ~] = unique(spikes.vs_params(:,14));

if length(led_cond)>1 
    trial_per_power = min(diff(index));

    spikes_led = cell(1, length(led_cond));
   
    for i=1: length(led_cond)

       if isempty(assign)
            spikes_led{i} = filtspikes(spikes, 0, 'trials', [index(i):(index(i)+trial_per_power)]);
       else
            spikes_led{i} = filtspikes(spikes, 0, 'trials', [index(i):(index(i)+trial_per_power)], 'assigns', assign);
       end

       psth_bl(spikes_led{i}, [], [], [], colour(i,:))

    end
elseif rem(length(spikes.vs_params(:,1)), num_trials_cond) == 0
    num_led_cond = length(spikes.vs_params(:,1))/num_trials_cond;
    led_cond = 1:num_led_cond;
    spikes_led = cell(1, num_led_cond);
    
    for i = 1:num_led_cond
       if isempty(assign)
            spikes_led{i} = filtspikes(spikes, 0, 'trials', [(num_trials_cond*i - (num_trials_cond -1)):(num_trials_cond*i)]);
       else
            spikes_led{i} = filtspikes(spikes, 0, 'trials', [(num_trials_cond*i - (num_trials_cond -1)):(num_trials_cond*i)], 'assigns', assign);
       end

       psth_bl(spikes_led{i}, [], [], [], colour(i,:))
    end
end

lgd = legend(num2str(led_cond));
lgd.FontSize = 13;
lgd.Title.String = 'Led Power'; 

title(['PSTH led power : ' num2str(assign)]);

%% plot first led spike hist
if ~ isempty(assign)
    figure;
    ax = nan(1, length(led_cond));

    bin_psth = 50; %ms
    bin_LED = 1; %ms
    bin_limit = [0 25]; %ms

    for  i=1: length(led_cond)
       ax(i) = subplot(1, length(led_cond), i); 
       [length_led_1, led_spiketimes_1, led_hist_1] = make_first_led_spiketimes(spikes_led{i}, assign, bin_LED, bin_limit, 1, 1); 
       title(['led power = ' num2str(led_cond(i))]);
    end

    linkaxes(ax,'xy')
end
