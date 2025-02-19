% This script takes trial-structured fluorescence traces and returns 
% inferred neural spiking rate traces (with the same structure).
% It finds the best minimal penalty, lambda (required for noise removal), that fits the data.
% Knowing gamma (the calcium decay in a time bin) is required to run the script.
% Time bins should be longer than the typical rise time of the calcium.

close
clear

load("one_pix_dff_trials.mat")                     %%%%%%%%%%%%%%%%%%%%%%% UPLOAD YOUR OWN

frames=size(whisker_data,1);
trials=size(whisker_data,2);

dt = 0.1; %sec
trial_time = (frames)*dt;

t_trial = dt:dt:trial_time;

% gamma - the percentage of calcium left after a single recording step
% for gcamp6s
gamma_100hz = 0.95; % here, place your own          %%%%%%%%%%%%%%%%%%%%%%% PLACE YOUR OWN

gamma_d = gamma_100hz^2;

t_trace = frames; % should be an even number, 400 and above recomended
odd_indx = 1:2:t_trace-1;
even_indx = 2:2:t_trace;

% for later use (to rebuid the calcium)
Dinv_d = zeros(t_trace/2); 
insert_vec = 1;
for k = 1:t_trace/2
    Dinv_d(k,1:k) = insert_vec;
    insert_vec = [gamma_d^k, insert_vec];
end

% all penalty levels (sort of smoothing)
% If needed, additional lambdas cab be added to the list to refine the search for best minimal lambda.
all_lambda = [1000000 100000 10000 5000 2000 1000 500 300 200 150 100 50 20 10 1 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0001];

% The following finds the right deconvolution penalty (needs to be done once)
% Load data; in the example the data is stored in the "act" parameter.
% It needs to be:
% 1. one long column vector 
% 2. free of inf and nan 
% It can be numbers between 10^0 and 10^2 (runs faster and more accurate)
% the next step fits the data we have to the required structure 

act_matrix = whisker_data*1000; % time x "trials"  %%%%%%%%%%%%%%%%%%%%%%% YOUR DATA GOES HERE 

% divide to odd and even
act_odd = act_matrix(odd_indx,:);
act_even = act_matrix(even_indx,:);

% deconvolve with different penalties
parfor i_lambda = 1:length(all_lambda)
     lambda = all_lambda(i_lambda);
     [r_inferred_convar_odd,beta0_odd,r0_odd] = convar(act_odd,gamma_d,lambda);
     [r_inferred_convar_even,beta0_even,r0_even] = convar(act_even,gamma_d,lambda);

    % rebuild calcium
    c_inferred_odd = Dinv_d*[r0_odd;r_inferred_convar_odd];
    c_inferred_even = Dinv_d*[r0_even;r_inferred_convar_even];

    % genearate expectations
    predicts_c_even = (c_inferred_odd(1:end-1,:)+c_inferred_odd(2:end,:))/2;
    predicts_c_odd = (c_inferred_even(1:end-1,:)+c_inferred_even(2:end,:))/2;

    % shift to mean substructed traces
    predicts_c_even_nodc = predicts_c_even-repmat(mean(predicts_c_even,1),t_trace/2-1,1);
    predicts_c_odd_nodc = predicts_c_odd-repmat(mean(predicts_c_odd,1),t_trace/2-1,1);
    act_even_nodc = act_even-repmat(mean(act_even,1),t_trace/2,1);
    act_odd_nodc = act_odd-repmat(mean(act_odd,1),t_trace/2,1);

    % check deviations from predictions
    err_cf_oe = (predicts_c_even_nodc-act_even_nodc(1:end-1,:)).^2;
    mean_err_cf_oe_per_trial = (2/t_trace)*sum(err_cf_oe,1);
    err_cf_eo = (predicts_c_odd_nodc-act_odd_nodc(2:end,:)).^2;
    mean_err_cf_eo_per_trial = (2/t_trace)*sum(err_cf_eo,1);
    mean_err_cf_per_trial = (mean_err_cf_oe_per_trial+mean_err_cf_eo_per_trial)/2;
    std_err_cf_per_trial = sqrt(((mean_err_cf_oe_per_trial-mean_err_cf_per_trial).^2+(mean_err_cf_eo_per_trial-mean_err_cf_per_trial).^2)/2);
    mean_err_cf(i_lambda) = mean(mean_err_cf_per_trial);
    mean_std_err_cf(i_lambda) = mean(std_err_cf_per_trial);
end

[min_err,lambda_index_at_min] = min(mean_err_cf);
figure(1)
loglog(all_lambda,mean_err_cf)
hold on
loglog(all_lambda,mean_err_cf*0+min_err+mean_std_err_cf(lambda_index_at_min),'k--')
xlabel('\lambda')
ylabel('error')
axis([10^-4 10^4 0 2000])
title('mean c to flu error odd-even')

best_min_lambda = all_lambda(find(mean_err_cf-(min_err+mean_std_err_cf(lambda_index_at_min))<0,1));

% if you have multiple long traces from the same experiment , you can loop over them and average over the scale of lambda to find best lamba
% lambda = 10^mean(log10(all_lambda). remember to clear mean_err_cf and
% clear mean_std_err_cf before each looping


%% now we infer

lambda = best_min_lambda; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% YOUR LAMBDA GOES HERE
gamma_100hz = 0.95;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% YOUR GAMMA GOES HERE 
dt = 0.1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  YOUR TIME BIN GOES HERE 
t = (1:t_trace)*dt;

[rates,beta0,r0] = convar(act_matrix,gamma_100hz,lambda);

%rates = [rates(2,:); rates]; % r(1) has no biological meaning

%% plot results

mean_rates = zeros(size(rates,1),1);
figure(2)

subplot(2,1,1)
plot(t_trial,whisker_data)
xlabel('time (sec, stim onset time 1)')
ylabel('dFF (normalized per trial)')
title('200 trials plotted from a whisker field pixel')
hold on
%plot average
plot(t_trial,mean(whisker_data,2),'k','LineWidth',2)


subplot(2,1,2)
for i = 1:size(rates,2)
    cur_rates = rates(:,i)-mean(rates(:,i));
    plot(t(2:end),cur_rates)
    mean_rates = mean_rates+cur_rates;
    hold on
end
mean_rates = mean_rates/size(rates,2);
xlabel('time (sec, stim onset time 1)')
ylabel('inferred rate w.r.t. mean rate')
title('200 trials plotted from an inferred whisker field pixel')
hold on
%plot average
plot(t(2:end),mean_rates,'k','LineWidth',2)