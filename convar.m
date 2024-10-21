function [r,beta0,r0] = convar(y,gamma,lambda)
% This function infers the rate r from the fluorescence recordings y, 
% assuming the calcium decays to gamma of its value each time bin if no spikes accrued, 
% and the penalty for a (squared) change in spiking rate is lambda 

% The function uses an analytical solution to calculate r and shift for
% positivity, including calculating beta0 and r at t=1 which is not a rate 
% but rather an artifact of this shift and the matrix notation 
% NOTATIONS: INPUTS y is t x n - time by traces matrix ; 
% gamma a number between 0 and 1 (typically close to 1) the calcium decay between two measurement points; 
% lambda a number, weight penalty of the inferred rate fluctuations 
% OUTPUTS r is t-1 x n matrix ; r0 and beta0 1 x n vectors

t = size(y,1);
n = size(y,2);

Dinv = zeros(t);
insert_vec = 1;
for i_t = 1:t
    Dinv(i_t,1:i_t) = insert_vec;
    insert_vec = [gamma^i_t, insert_vec];
end

P = eye(t)-1/t*ones(t);
ytilde = P*y;
A = P*Dinv;
L = [zeros(t,1) [zeros(1,t-1); [zeros(1,t-1); [-eye(t-2), zeros(t-2,1)] + [zeros(t-2,1), eye(t-2)]]]];
Z = L'*L;

multiplies_r = (A'*A+lambda*Z);
multiplies_r = pinv(multiplies_r);
r_analytic = multiplies_r*A'*ytilde;

d = -min([r_analytic(2:end,:); zeros(1,n)],[],1);
rm_ratio = A'*A*ones(t,1)./(A'*A*[1; zeros(t-1,1)]);
rm = -d*rm_ratio(1);
r = r_analytic+repmat(d,t,1)+[rm;zeros(t-1,n)];

beta0 = 1/t*ones(1,t)*(y-Dinv*r);

% internal verification, for example , if needed, all next fuctions should yeild the
% same number
% a1 = (y(:,19)-ones(t,1)*beta0(1,19)-Dinv*r_anlytic(:,19))'*(y(:,19)-ones(t,1)*beta0(1,19)-Dinv*r_anlytic(:,19));
% a2 = (ytilde(:,19)-A*r_anlytic(:,19))'*(ytilde(:,19)-A*r_anlytic(:,19));
% a3 = (ytilde(:,19)-A*(r_anlytic(:,19)+repmat(d(1,19),t,1)+[r0(1,19);zeros(t-1,1)]))'*(ytilde(:,19)-A*(r_anlytic(:,19)+repmat(d(1,19),t,1)+[r0(1,19);zeros(t-1,1)]));
% a4 = (ytilde(:,19)-A*r(:,19))'*(ytilde(:,19)-A*r(:,19));

r0 = r(1,:);
r = r(2:end,:);

end

