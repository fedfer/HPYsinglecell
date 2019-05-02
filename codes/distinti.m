function [Xs Ts XTs ns qs ls indici_XTs] = distinti(X, T)
      
% Funzione che trova la frequenza dei distinti per i piatti X e per i
% tavoli T

% gli output sono:
% XTs= ? la matrice contenete le coppie (piatto tavolo) dei valori distinti
% in XT=[X;T]
% Xs= piatti distinti serviti nel ristorante
% Ts= tavoli distinti nel ristorante
% ns= numerosit? dei distinti Xs
% qs= numerosit? dei distinti Ts
% ls= numerosit? dei tavoli in cui viene servito un certo piatto

% X_ordinati ? il vettore X ordinato in modo crescente
% T_ordinati il corrispondente. Uso la funzione sortrows : lavora sulle
% righe e ordina prima gli elementi in base alla prima colonna e poi a
% parit? di valori sulla prima colonna ordina in base alla seconda colonna
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
% riordino il Ts in base all'ordine degli X : attenzione devo anche
% riordinare qs per questo salvo gli indici
[indici_T indici_qs]=sort(indici_T);
Ts=T_ordinati(indici_T);
qs=qs(indici_qs);

XTs=[X_ordinati(indici_T);Ts];

% per trovare ls conto i distinti in XTs(1,:)
[Ls indici_L]=unique(XTs(1,:),'legacy');
ls=diff([0 indici_L]);


% t=1;
% qs=zeros(1,n);
% XTs=
% ls=ones(1,k); % frequenze dei tavoli distinti
% for i=1:k
%     for j=1:ns(i)
%     Ts(t)=T_locale(j);
%     qs(t)=qs(t)+1;
%     if j<ns(i) && T_locale(j+1) ~= T_locale(j);
%        t=t+1;
%        ls(i)=ls(i)+1;
%     end  
%     end
%     t=t+1;
% end
% L=size(Ts);
% L=L(2);
% qs=qs(1:L);
% Ts=Ts(1:L);
% 
% % creo un vettore che contenga i tavoli distinti accoppiati al piatto
% % distinto
% XT=[X; T];
% XTs=[Xs(1)*ones(1,ls(1));Ts(1:ls(1))];
% cumls=cumsum(ls);
% for i=1:k-1
%     XTs=[XTs [Xs(i+1)*ones(1,ls(i+1));Ts((cumls(i)+1):cumls(i+1))]];
% end

indici_XTs=0;
% creo un vettore indici_XTs che indica le posizione degli elementi di XT,
% cio? del vettore contenete le coppie (piatto, tavolo), in XTs, ossia nel
% vettore dei distinti (piatto, tavolo)
% for i=1:n;
%     for j=1:L;
%         if XT(1:2,i)==XTs(1:2,j);
%             indici_XTs(i)=j;
%         end
%     end
% end

end

