% MH for sigma
function [sigma_new acc]=MH_sigma_hyper(sigma_old,theta,q1_star,a0,b0)

k1=length(q1_star);
acc=0;

sigma_p=unifrnd(0,1);

vec1=1:(k1-1);
rho=exp(sum(log(theta+vec1*sigma_p))+sum(gammaln(q1_star-sigma_p))-k1*gammaln(1-sigma_p)...
    -sum(log(theta+vec1*sigma_old))-sum(gammaln(q1_star-sigma_old))+k1*gammaln(1-sigma_old));
% add hyperprior contribution, \sigma_j | a_0,b_0 \sim Beta(a_0,b_0)
rho= rho + (a0-1)*log(sigma_p) + (b0-1)*log(1-sigma_p)...
    - (a0-1)*log(sigma_old) + (b0-1)*log(1-sigma_old);
rho=min(1,rho);


U=unifrnd(0,1);
if U<=rho
    sigma_new=sigma_p;
    acc=1;
else
    sigma_new=sigma_old;
end
