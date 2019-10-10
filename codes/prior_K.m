% function to compute the probability distribution a priori for K_n (number of clusters)
% output: probability of K_n=1,2,...,n


function P=prior_K(n,N)

a=11;
a0=11;
b=7;
b0=7;
valori_sigma=0.01:0.05:0.99;


LogC=cell(1,length(valori_sigma));
for i=1:length(valori_sigma)
LogC{i}=generalized_factorial(n,n,valori_sigma(i));
i
end


M_P=zeros(N,n);
for j=1:N
    % generates theta, sigma, theta_0 sigma_0 from prior
    theta_0=gamrnd(a0,b0);
    theta=gamrnd(a,b);
    [sigma indice_sigma]=genera_sigma(valori_sigma);
    [sigma_0 indice_sigma0]=genera_sigma(valori_sigma);

    vettore_prodotti=theta_0+sigma_0*(1:n-1);
    for k=1:n

         q=k:n;
         somma=sum(exp(LogC{indice_sigma}(n+1,k+1:n+1)+LogC{indice_sigma0}(k+1:n+1,k+1)'+gammaln(theta/sigma+q)+gammaln(theta_0+1)-...
                 gammaln(theta_0+q)-gammaln(theta/sigma+1)-gammaln(theta+n)-k*log(sigma_0)+...
                 sum(log(vettore_prodotti(1:(k-1))))+gammaln(theta+1)));
        M_P(j,k)=somma/sigma;
     end
    j
end
%  approximate prior
P=sum(M_P)/N;

% plot the distribution 
plot(1:n,P,'-');