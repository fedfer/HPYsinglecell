function [Xs ns]=cluster_fct(X)

% compute the distinct values and relative frequencies

[Xs indexes]=unique(sort(X),'legacy');
% Xs is the set of distict values, index(k) contain the indeces of the last occurrence of Xs in X
ns=diff([0 indexes]);