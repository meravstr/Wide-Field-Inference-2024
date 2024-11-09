% This script takes a dffed fluorescence trace and 
% returns an inferred neural spiking rate trace.
% It finds the penalty, lambda (required for noise removal), 
% that maintains the fluctuations structure of the data, while removing a significant amount of noise 
% Knowing gamma (the calcium decay in a time bin) is required to run the script.
% Time bins should be longer than the typical rise time of the calcium.

% If the fluorescence trace includes 50,000 points or more, 
% the script is expected to perform well. 
% Inferring a trace with 10,000 points or less is possible, 
% but reducing t_trace to 400 is recommended.

close all
clear

% gamma - the percentage of calcium left after a single recording step
% for gcamp6s
gamma_10hz = 0.95; 
gamma = gamma_10hz^(1/1.5); % PLACE YOUR OWN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_trace = 1000; % should be an even number, 400 and above are recomended

% for later use (to rebuid the calcium)
Dinv = zeros(t_trace); 
insert_vec = 1;
for k = 1:t_trace
    Dinv(k,1:k) = insert_vec;
    insert_vec = [gamma^k, insert_vec];
end

% all penalty levels (smoothing in a way)
% This can be change to refine the search for best minimal lambda, if
% needed
all_lambda = [1000000 100000 10000 5000 2000 1000 500 300 200 150 100 50 20 10 1 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0001];

load('another_example_fluorescence_audio_gcamp6.mat');

% find the right deconvolution penalty (needs to be done once)
% load data, the data is stored in the "act" parameter in the file
% it needs to be:
% 1. one long column vector 2. free of inf and nan 3. it is easier to work
% with numbers between 10^0 and 10^2
% the next step fits the data we have to the required structure 
act = act(1,1001:end)'*1000; % YOUR DATA GOES HERE %%%%%%%%%%%%%%%%%%%%%%%

n_traces = floor(length(act)/t_trace);
act = act(1:n_traces*t_trace);
act_matrix = reshape(act,t_trace,n_traces); % time x "trials"

% measures fluctuations level of the data
fluc_act = mean(diff(act_matrix,[],1).^2,1)';

% deconvolve with different penalties
for i_lambda = 1:length(all_lambda)
     lambda = all_lambda(i_lambda);
     [r_inferred_convar,beta0,r0] = convar(act_matrix,gamma,lambda);
    % rebuild calcium
    c_convar = Dinv*[r0;r_inferred_convar];
    % flactuations
    fluc_c(i_lambda,:) =  mean((c_convar(2:end,:)-c_convar(1:end-1,:)).^2,1);
    
    % analyzing fluc correlations
    curr_coef = corrcoef(fluc_act,fluc_c(i_lambda,:));
    curr_coef = curr_coef(1,2);
    cf_corr_con(i_lambda)  = curr_coef;
end

figure
loglog(all_lambda,cf_corr_con)
[peak_val,peak_indx] = findpeaks(cf_corr_con);
hold on
loglog(all_lambda(peak_indx),peak_val,'o')

% if you have multiple long traces from the same experiment , you can loop over them and average over the scale of lambda to find best lamba
% lambda = 10^mean(log10(all_lambda). remember to clear cf_corr_con
%% now we infer

lambda = all_lambda(peak_indx); % YOUR LAMBDA GOES HERE %%%%%%%%%%%%%%%%%%
t_trace = 1000; % should be a multiplication of 8 
gamma_10hz = 0.95; 
gamma = gamma_10hz^(1/1.5);  % YOUR GAMMA GOES HERE %%%%%%%%%%%%%%%%%%%%%%
dt = 1/15; % YOUR time bin GOES HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (1:t_trace)*dt;

n_traces = floor(length(act)/t_trace);
act_matrix_1 = reshape(act,t_trace,n_traces);
act_matrix_2 = reshape(act(t_trace/2+1:end-t_trace/2),t_trace,n_traces-1);
[r_1,beta0_1,r0_1] = convar(act_matrix_1,gamma,lambda);
[r_2,beta0_2,r0_2] = convar(act_matrix_2,gamma,lambda);

r = zeros(length(act),1);
r(1) = 0;
r(2:3*t_trace/4) = r_1(1:3*t_trace/4-1,1);

%% aligning and concatinating
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

%% lets plot examples
t_trace_short = 200;
t = (1:1:t_trace_short)/15; %15hz
for i = 1:10
figure
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