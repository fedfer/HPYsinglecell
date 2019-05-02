function [Xs ns]=clusterizza(X)

% calcolo i valori distinti in X  e le frequenze

[Xs indici]=unique(sort(X),'legacy');
% Xs ? l'insieme dei valori distinti, invece indici(k)  contiene gli indici
% del'ultima occorrenza di Xs in X. 
ns=diff([0 indici]);