function [r_final,r1,beta_0,iter] = dynbin_wstop(y,gamma,lambda,r_initial,err_or_iter)
% This function infers the rate r from the fluorescence recordings y, 
% assuming the calcium decays to gamma of its value each time bin if no spikes accrued, 
% and the penalty for an L1 cost on change in spiking rate is lambda 

% The function implements the algorithm dynamically binned to calculate r 
% including calculating beta0 and r at t=1 which is not a rate 
% but rather an artifact of the matrix notation 

% NOTATIONS: INPUTS y is t x n - time by traces matrix ; 
% gamma a number between 0 and 1 (typically close to 1); the calcium decay between two measurement points
% lambda a number, weight penalty of the inferred rate fluctuations
% r_initial - a T-1xn guess for the solutions to r_final (if empty [] a random guess is used)
% err_or_iter - number, this parameter defines the stopping criterion. If
% the number is larger than 1, it is used as the number of iterations to
% perform. If the number is a fraction, a stopping criterion is defined according to this fraction.
% In this case, the iterations would stop once the change in the rate on
% average is smaller than the chosen fraction from the mean magnitude of the rate r (not including its constant shift above zero, if it exists). 
% If empty, a default of 1000 iterations is used.
% OUTPUTS r is t-1 x n matrix ; r0 and beta0 the offsets 1 x n vectors

t = size(y,1);
n = size(y,2);

Dinv = zeros(t);
insert_vec = 1;
for i_t = 1:t
    Dinv(i_t,1:i_t) = insert_vec;
    insert_vec = [gamma^i_t, insert_vec];
end

P = eye(t)-1/t*ones(t);
tildey = P*y;
A = P*Dinv;
% largest step size that ensures converges
s = 0.5*((1-gamma)/(1-gamma^t))^2;

% initializing
if isempty(r_initial)
    r = rand(size(y));
else % finding a good r1 to start with, according to the guess r_initial and y
    r = zeros(size(y));
    a = A(:,1)'*A(:,1);
    b = -2*(A(:,1)'*tildey-A(:,1)'*A(:,2:t)*r_initial);
    c = diag((tildey-A(:,2:t)*r_initial)'*(tildey-A(:,2:t)*r_initial));
    %for j = 1:size(y,2)
            delta = b^2-4*a*c;%b(j)^2-4*a*c(j);
            %if delta>0
            %    j
                r(1,:) = (-b+sqrt(delta))/(2*a);%(-b(j)+sqrt(delta))/(2*a);
            %else 
            %    r(1,j) = y(1,j);
            %end
    %end
    r(2:end,:) = r_initial;
end

if isempty(err_or_iter)
    err_or_iter = 1000;
end

if err_or_iter>1
    for i = 1:ceil(err_or_iter)
        r_old = r;
        Ar = A*r;
        tmAr = (tildey-Ar);
        At_tmAr = A'*tmAr;
        x = r + s*At_tmAr;
        for j = 1:size(y,2)
            r(2:end,j) = fTVdenoise(s*lambda,x(2:end,j));
        end
        r(r<0) = 0;
        r(1,:) = x(1,:);
    end
    
else
    i = 1;
    relative_change = zeros(size(y,2),1);
    test_err = 1; % precentage of error compare to rate magnitude
                  % running, with stopping according to error or iterations
    indx = 1:size(y,2);
    to_update = ones(size(y,2),1);
    while test_err > (err_or_iter/2)
        r_old = r;
        Ar = A*r;
        tmAr = (tildey-Ar);
        At_tmAr = A'*tmAr;
        x = r + s*At_tmAr;
        for j = 1:size(y,2)
            if indx(to_update)
            r(2:end,j) = fTVdenoise(s*lambda,x(2:end,j));
            else 
            r(2:end,j) = r_old(2:end,j);
            end
        end
        r(r<0) = 0;
        for j = indx(to_update)
           relative_change(j) = mean(abs(r_old(2:end,j)-r(2:end,j)))/mean((r(2:end,j)-min(r(2:end,j))));
           if relative_change(j) < (err_or_iter/2)
               to_update(j) = 0;
           end
        end
        r(1,:) = x(1,:);
        i = i+1;
        test_err = max(relative_change);
    end
end
    
    
r_final = r(2:end,:);
r1 = r(1,:);
beta_0 = mean(y-Dinv*r);
iter = i;
end

