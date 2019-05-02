% genera i valori di sigma da una uniforme discretizzata

function [sigma indice_sigma]=genera_sigma(valori_sigma)

U=unifrnd(0,1);
indice_sigma=find(U<=valori_sigma);
% se supero il valore della discretizzazione assegno a sigma l'ultimo
% valore disponibile
if isempty(indice_sigma);
    sigma=valori_sigma(end);
    indice_sigma=length(valori_sigma);
    return
end
indice_sigma=indice_sigma(1);
sigma=valori_sigma(indice_sigma);