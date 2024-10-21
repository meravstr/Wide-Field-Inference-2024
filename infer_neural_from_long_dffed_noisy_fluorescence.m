% This script takes a dffed fluorescence trace and 
% REUTRNS AN INFERRED NUERAL SPIKING RATE activity trace.

% It finds the best minimal lambda (a penalty required for noise removal)
% that fits the data.

% To run the script one needs to KNOW GAMMA (the calcium decay in a bin).

% Bins should be longer than the typical rise time of the calcium.
% The longer your recordings are the better
% (50,000 points is great. 10,000 is possible but reducing t_trace to 400
% is recomended).

close all
clear

% Gamma - the percentage of calcium left after a single recording step,
% for gcamp6s in 10hz recordings:
gamma_10hz = 0.95; 
gamma = gamma_10hz^(1/1.5); % PLACE YOUR OWN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_d = gamma^2;

t_trace = 1000; % should be an even number, 400 and above are recomended
odd_indx = 1:2:t_trace-1;
even_indx = 2:2:t_trace;

% For later use (to rebuild the calcium)
Dinv_d = zeros(t_trace/2); 
insert_vec = 1;
for k = 1:t_trace/2
    Dinv_d(k,1:k) = insert_vec;
    insert_vec = [gamma_d^k, insert_vec];
end

% All penalty levels (sort of smoothing):
% *This list can be changed to refine the search for best minimal lambda
all_lambda = [1000000 100000 10000 5000 2000 1000 500 300 200 150 100 50 20 10 1 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0001];

% The following code finds the right inference penalty (needs to be done once)

% Load data; the data is stored in the "act" parameter in the file
% it needs to be:
% 1. one long column vector 2. free of inf and nan 3. it is easier to work
% with numbers between 10^0 and 10^2
load('example_single_pixel_fluorescence_gcamp6s_15hz.mat');
% The next step fits the data we have to the required structure 
act = act'*1000; % YOUR DATA GOES HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_traces = floor(length(act)/t_trace);
act = act(1:n_traces*t_trace);
act_matrix = reshape(act,t_trace,n_traces); % time x "trials"

% Dividing to odd and even
act_odd = act_matrix(odd_indx,:);
act_even = act_matrix(even_indx,:);

% Inferrin with different penalties
for i_lambda = 1:length(all_lambda)
     lambda = all_lambda(i_lambda);
     [r_inferred_convar_odd,beta0_odd,r0_odd] = convar(act_odd,gamma_d,lambda);
     [r_inferred_convar_even,beta0_even,r0_even] = convar(act_even,gamma_d,lambda);

    % Rebuilding the calcium
    c_inferred_odd = Dinv_d*[r0_odd;r_inferred_convar_odd];
    c_inferred_even = Dinv_d*[r0_even;r_inferred_convar_even];

    % Genearating expectations
    predicts_c_even = (c_inferred_odd(1:end-1,:)+c_inferred_odd(2:end,:))/2;
    predicts_c_odd = (c_inferred_even(1:end-1,:)+c_inferred_even(2:end,:))/2;

    % Shifting to mean substructed traces
    predicts_c_even_nodc = predicts_c_even-repmat(mean(predicts_c_even,1),t_trace/2-1,1);
    predicts_c_odd_nodc = predicts_c_odd-repmat(mean(predicts_c_odd,1),t_trace/2-1,1);
    act_even_nodc = act_even-repmat(mean(act_even,1),t_trace/2,1);
    act_odd_nodc = act_odd-repmat(mean(act_odd,1),t_trace/2,1);

    % Checking deviations from predictions
    err_cf_oe = (predicts_c_odd_nodc-act_even_nodc(2:end,:)).^2;
    mean_err_cf_oe_per_trial = (2/t_trace)*sum(err_cf_oe,1);
    err_cf_eo = (predicts_c_even_nodc-act_odd_nodc(2:end,:)).^2;
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

% If you have multiple long traces from the same experiment , you can loop over them and average over the scale of lambda to find best lamba
% lambda = 10^mean(log10(all_lambda). 
% Remember to clear mean_err_cf and mean_std_err_cf before each loop


%% Now we infer

lambda = 0.1; % YOUR BEST LAMBDA GOES HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_trace = 1000; % should be a multiplication of 8 
gamma_10hz = 0.95; 
gamma = gamma_10hz^(1/1.5); % YOUR GAMMA GOES HERE %%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/15; % YOUR time bin GOES HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (1:t_trace)*dt;

n_traces = floor(length(act)/t_trace);
act_matrix_1 = reshape(act,t_trace,n_traces);
act_matrix_2 = reshape(act(t_trace/2+1:end-t_trace/2),t_trace,n_traces-1);
[r_1,beta0_1,r0_1] = convar(act_matrix_1,gamma,lambda);
[r_2,beta0_2,r0_2] = convar(act_matrix_2,gamma,lambda);

r = zeros(length(act),1);
r(1) = 0;
r(2:3*t_trace/4) = r_1(1:3*t_trace/4-1,1);
% Aligning and concatenating
for i_n = 1:n_traces-1
    mean_1r = mean(r_1(t_trace/2+t_trace/8+1:t_trace/2+t_trace/8+t_trace/4,i_n),1);
    means_2l = mean(r_2(t_trace/8+1:t_trace/8+t_trace/4,i_n),1);
    r_2(:,i_n) = r_2(:,i_n) + mean_1r - means_2l;
    mean_1l = mean(r_1(t_trace/8+1:t_trace/8+t_trace/4,i_n+1),1);
    mean_2r = mean(r_2(t_trace/2+t_trace/8+1:t_trace/2+t_trace/8+t_trace/4,i_n),1);
    r_1(:,i_n+1) = r_1(:,i_n+1) + mean_2r - mean_1l;
    % using r_2 to smooth transitions between parts
    r(3*t_trace/4+1+(i_n-1)*t_trace:3*t_trace/4+t_trace/2+(i_n-1)*t_trace) = r_2(t_trace/4:3*t_trace/4-1,i_n);
    r(2+i_n*t_trace:3*t_trace/4+i_n*t_trace) = r_1(1:3*t_trace/4-1,(i_n+1));
end

r = r+max(0,-min(r)); % HERE IS YOUR NEURAL ACTIVITY
r(1) = r(2); % r(1) has no biological meaning

%% Let us plot examples
t_trace_short = 200;
t = (1:t_trace_short)*dt;
for i = 1:8
figure(i+1)
subplot(2,1,1)
plot(t,act(t_trace_short*i+1:t_trace_short*(i+1)),'LineWidth',2, 'Color',[0.7 0.7 0.7])
xlabel('time[sec]')
ylabel('fluorescence')
subplot(2,1,2)
plot(t,r(t_trace_short*i+1:t_trace_short*(i+1)),'LineWidth',2,'Color',[0.3 0.3 1])
xlabel('time[sec]')
ylabel('inferred spiking rate')
box('off') 
end

%% For the example: averaging according to events (and comparing to fluorescence)
% The data includes two types of events, the following finds the mean and
% standard deviation of the fluorescence and inferred spiking rate with
% respect to these events

load('events.mat');

frames_prior = 15;
frames_post = 30;
frames_shift = 1000;
total_frames = 49000;
t = (-frames_prior:1:frames_post)*dt;

act1 = zeros(frames_prior+frames_post+1,length(event1));
r1 = zeros(frames_prior+frames_post+1,length(event1));
for i_trial = 1:length(event1)
    act1(:,i_trial) = act(event1(i_trial)-frames_prior-frames_shift:event1(i_trial)+frames_post-frames_shift);
    r1(:,i_trial) = r(event1(i_trial)-frames_prior-frames_shift:event1(i_trial)+frames_post-frames_shift);
end
act1_mean = mean(act1,2);
act1_std = std(act1,[],2);
r1_mean = mean(r1,2);
r1_std = std(r1,[],2);

figure(10)
subplot(2,1,1)
plot(t,act1_mean,'b')
hold on
fill([t fliplr(t)],[act1_mean+act1_std;flipud(act1_mean-act1_std)],'b','FaceAlpha',0.2,'EdgeColor','none')
subplot(2,1,2)
plot(t,r1_mean,'b')
hold on
fill([t fliplr(t)],[r1_mean+r1_std;flipud(r1_mean-r1_std)],'b','FaceAlpha',0.2,'EdgeColor','none')

act2 = zeros(frames_prior+frames_post+1,length(event2));
r2 = zeros(frames_prior+frames_post+1,length(event2));
for i_trial = 1:length(event2)
    act2(:,i_trial) = act(event2(i_trial)-frames_prior-frames_shift:event2(i_trial)+frames_post-frames_shift);
    r2(:,i_trial) = r(event2(i_trial)-frames_prior-frames_shift:event2(i_trial)+frames_post-frames_shift);
end
act2_mean = mean(act2,2);
act2_std = std(act2,[],2);
r2_mean = mean(r2,2);
r2_std = std(r2,[],2);

figure(10)
subplot(2,1,1)
plot(t,act2_mean,'r')
hold on
fill([t fliplr(t)],[act2_mean+act2_std;flipud(act2_mean-act2_std)],'r','FaceAlpha',0.2,'EdgeColor','none')
ylabel('fluorescence [a.u.]')
xlabel('t[sec]')
legend('event 1 mean','event 1 standard deviation.','event 2 mean','event 2 standard deviation','Location','northwest')
%plot((0:0.03:50)*0,0:0.03:50,'--k')
subplot(2,1,2)
plot(t,r2_mean,'r')
hold on
fill([t fliplr(t)],[r2_mean+r2_std;flipud(r2_mean-r2_std)],'r','FaceAlpha',0.2,'EdgeColor','none')
plot((20:0.03:30)*0,20:0.03:30,'--k')
%legend('event 1 mean','event 1 standard deviation.','event 2 mean','event 2 standard deviation','Location','northwest')
xlabel('t[sec]')
ylabel('neuroal spiking rate [a.u.]')

