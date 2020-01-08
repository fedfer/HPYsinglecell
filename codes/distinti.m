function [Xs Ts XTs ns qs ls indici_XTs] = distinti(X, T)

% Function that finds the freq of distincts values for dishes X and tables T

% Outputs:
% XTs: matrix with pairs of dish-table of distinct values in XT = [X;T]
% Xs: distinct dishes in restaurant
% Ts= distinct tables in restaurant
% ns= number of uniques in Xs
% qs= number of uniques in Ts
% ls= number of tables in which a certain dish is served


XT_ordinati=sortrows([X; T]')';
X_ordinati= XT_ordinati(1,:);
T_ordinati= XT_ordinati(2,:);

% creo due vettori:
% Xs = vettore costituito dai valori distinti contenuti in oX (cio?? in  X).
% Tali valori sono ordinati.
% ns = vettore con le frequenze con cui compaiono i valori di Xs in X.

[Xs indici_X]=unique(X_ordinati,'legacy');
ns=diff([0 indici_X]);

[Ts indici_T]=unique(T_ordinati);
qs=histc(T_ordinati,Ts);

[indici_T indici_qs]=sort(indici_T);
Ts=T_ordinati(indici_T);
qs=qs(indici_qs);

XTs=[X_ordinati(indici_T);Ts];

[Ls indici_L]=unique(XTs(1,:),'legacy');
ls=diff([0 indici_L]);



indici_XTs=0;


end

