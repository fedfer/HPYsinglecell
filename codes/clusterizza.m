function [Xs ns]=clusterizza(X)

% compute the distinct values and relative frequencies

[Xs indici]=unique(sort(X),'legacy');
% Xs is the set of distict values, index(k) contain the indeces of the last occurrence of Xs in X
ns=diff([0 indici]);