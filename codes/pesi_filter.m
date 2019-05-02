% funzione per il calcolo dei pesi (gi? normalizzati) nel particle filter 
function p=pesi_filter(mjk,m_j_dot,m_dd,m_dot_k,n_j_dot_k,J,nn,medie_mx,bigK,N_iter)
%MANCANO DA METTERE LE PRIOR ANCGE SE NON SERVONO A
%MOLTO!!!!!!!!!!!!XXXXXXXXXXXXXXXXXXX
vec=1:(bigK-1);
% vettore dei logaritmi delle probabilit?
%N_iter ? il numero di elementi nella MC importance sample
logp=zeros(1,N_iter);
for jj=1:N_iter
    alpha=medie_mx(jj,1:J);
    d=medie_mx(jj,(J+1):(2*J));
    gamma=medie_mx(jj,2*J+1);
    nu=medie_mx(jj,2*J+2);
Phi=0;
for j=1:J
    % calcolo l'EPPF relativa al ristorante j esimo (explog)
    vec_j=1:(m_j_dot(j)-1);
    % calcolo il massimo di n_j.k perch? lo uso nel calcolo del coefficiente
    % fattoriale generalizzato: cos? non genero cose che non mi servono
    max_nj=max(n_j_dot_k(j,:));
    max_mj=max(mjk(j,:));
    % genero tutta la matrice dei coefficienti fattoriali:generando il
    % massimo genero pure tutti i coeffiecienti di ordine inferiore
    LogC=generalized_factorial(max_nj,max_mj,d(j));
    % mi calcolo ora il vettore dei coefficienti che mi interssa per il
    % calcolo della likelihood
    %big K ? il numero di unique values
    vec_loggfc=zeros(1,bigK);
    for ii=1:bigK
        vec_loggfc(ii)=LogC(n_j_dot_k(j,ii)+1,mjk(j,ii)+1);
    end
    Phi_j=sum(log(alpha(j)+vec_j*d(j)))-gammaln(alpha(j)+nn(j))+gammaln(alpha(j)+1)+...
        sum(vec_loggfc)-m_j_dot(j)*log(d(j));
    Phi=Phi+Phi_j;
end
logp(jj)=Phi+sum(log(gamma+nu*vec))-gammaln(gamma+m_dd)+gammaln(gamma+1)+...
    sum(gammaln(m_dot_k-nu))-bigK*gammaln(1-nu);
end
% calcolo i pesi: passo agli esponenziali e poi NORMALIZZO... Non faccio il
% ciclo ma lo calcolo in forma matriciale,  
% il ' fa il trasposto
M=ones(N_iter,1)*logp;
%sum(,2) fa la somma sulle colonne
p=1./sum(exp(M-M'),2)';
% for ii=1:N_iter
%     p(ii)=1./sum(exp(logp-logp(ii)));
% end