% MH update a_0,a_0 with proposal exp(1)
function [a0_new, b0_new, rho]=MH_a0b0(a0_old,b0_old,sigma,J)

% exponential proposal distribution
a0_p=exprnd(1);
b0_p=exprnd(1);


% a_0 and b_0 are also given EXP(1) prior so that the proposal and the
% prior cancel out in the MH ratio

rho = J*(gammaln(a0_p + b0_p)-gammaln(a0_p)-gammaln(b0_p))...
    + a0_p*sum(log(sigma)) + b0_p*sum(log(1-sigma)) ... 
    -J*(gammaln(a0_old + b0_old)-gammaln(a0_old)-gammaln(b0_old))...
    - a0_old*sum(log(sigma)) - b0_old*sum(log(1-sigma));
    
rho=min(1,rho);


U=unifrnd(0,1);
if U<=rho
    a0_new=a0_p;
    b0_new=b0_p;
else
    a0_new=a0_old;
    b0_new=b0_old;
end