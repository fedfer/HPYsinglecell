% MH update theta and theta_0 with proposal exp(1)
function [theta_new rho]=MH_theta_hyper(theta_old,sigma,Lj,nj,theta_h)

% exponential proposal distribution
theta_p=exprnd(1);

vec1=1:(Lj-1);
rho=exp(sum(log(theta_p+vec1*sigma))-gammaln(theta_p+nj)+gammaln(theta_p+1) ...
    -sum(log(theta_old+vec1*sigma))+gammaln(theta_old+nj)-gammaln(theta_old+1)-...
    theta_p+theta_old+log(exppdf(theta_old,1))-log(exppdf(theta_p,1)));

% add hyperprior contribution, \theta_j | \theta_h \gamma Beta(\theta_h,1)
rho= rho + (theta_h-1)*log(theta_p) - theta_p...
    - (theta_h-1)*log(theta_old) - theta_old;

rho=min(1,rho);


U=unifrnd(0,1);
if U<=rho
    theta_new=theta_p;
else
    theta_new=theta_old;
end