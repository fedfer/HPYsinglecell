% sample sigma from a discrete distribution 

function [sigma indice_sigma]=genera_sigma(valori_sigma)

U=unifrnd(0,1);
indice_sigma=find(U<=valori_sigma);
% if sigma is bigger than the last discrete point, I assign to it the biggest value
if isempty(indice_sigma);
    sigma=valori_sigma(end);
    indice_sigma=length(valori_sigma);
    return
end
indice_sigma=indice_sigma(1);
sigma=valori_sigma(indice_sigma);