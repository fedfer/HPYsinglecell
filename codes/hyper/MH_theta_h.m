% MH update theta_h 
function [theta_h_new, rho]=MH_theta_h(theta_h_old,theta,J)

% exponential proposal distribution
theta_h_p=exprnd(1);

% theta_h is given EXP(1) prior so that the proposal and the
% prior cancel out in the MH ratio

rho = - J*gammaln(theta_h_p) + theta_h_p*sum(log(theta)) ...
    + J*gammaln(theta_h_old) - theta_h_old*sum(log(theta)) ;
    
rho=min(1,rho);

U=unifrnd(0,1);
if U<=rho
    theta_h_new=theta_h_p;
else
    theta_h_new=theta_h_old;
end