% Metropolis-Hastings for sigma and sigma_0
function [sigma_new]=MH_sigma_independent(sigma_old,theta,n_star)

k=length(n_star);

% Use the uniform distribution as proposal
sigma_p=unifrnd(0,1);

vec=cumsum(ones(1,k-1));
rho=exp(sum(log(theta+vec*sigma_p))+sum(gammaln(n_star-sigma_p))-k*gammaln(1-sigma_p)...
    -sum(log(theta+vec*sigma_old))-sum(gammaln(n_star-sigma_old))+k*gammaln(1-sigma_old));
rho=min(1,rho);


U=unifrnd(0,1);
if U<=rho
    sigma_new=sigma_p;
else
    sigma_new=sigma_old;
end