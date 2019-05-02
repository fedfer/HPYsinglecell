% MH per l'aggiornamento di theta e theta_0 con proposal exp(1)
function [theta_new rho]=MH_theta(theta_old,sigma,Lj,nj)

% come proposal scelgo exponenziale

theta_p=exprnd(1);

vec1=1:(Lj-1);
rho=exp(sum(log(theta_p+vec1*sigma))-gammaln(theta_p+nj)+gammaln(theta_p+1) ...
    -sum(log(theta_old+vec1*sigma))+gammaln(theta_old+nj)-gammaln(theta_old+1)-...
    theta_p+theta_old+log(exppdf(theta_old,1))-log(exppdf(theta_p,1)));
rho=min(1,rho);


U=unifrnd(0,1);
if U<=rho
    theta_new=theta_p;
else
    theta_new=theta_old;
end
