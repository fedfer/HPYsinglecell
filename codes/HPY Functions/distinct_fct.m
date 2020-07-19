function [Xs Ts XTs ns qs ls indici_XTs] = distinct_fct(X, T)

% Function that finds the freq of distincts values for dishes X and tables T

% Outputs:
% XTs: matrix with pairs of dish-table of distinct values in XT = [X;T]
% Xs: distinct dishes in restaurant
% Ts= distinct tables in restaurant
% ns= number of uniques in Xs
% qs= number of uniques in Ts
% ls= number of tables in which a certain dish is served


XT_ord=sortrows([X; T]')';
X_ord= XT_ord(1,:);
T_ord= XT_ord(2,:);

% Creat to vectors
% Xs = containing the distinct values in X_ordered (same of  X).
% Thos values are ordered 
% ns = vector of frequencies of Xs in X.

[Xs index_X]=unique(X_ord,'legacy');
ns=diff([0 index_X]);

[Ts index_T]=unique(T_ord);
qs=histc(T_ord,Ts);

[index_T indici_qs]=sort(index_T);
Ts=T_ord(index_T);
qs=qs(indici_qs);

XTs=[X_ord(index_T);Ts];

[Ls indici_L]=unique(XTs(1,:),'legacy');
ls=diff([0 indici_L]);



indici_XTs=0;


end

