% funzione per il calcolo dei pesi (gi? normalizzati) nel particle filter
% relativi alle osservazioni nuove che vengono generate dalla normale
% multivariata ad ogni passo. Al passo jj campiono un nuovo valore degli
% iperparametri theta_k e conseguentemente creo un vettore di probabilit?
% p_new in cui metter? p(y|theta_k) e un vettore p in cui metto le probabilit?
%  p(y|m_k) (sono semplicemente un riordinamento di quelle in g che ho
%  calcolato prima
    
function [M_iperparametri_new_unnormalized omega]=pesi_filter_new(mjk,m_j_dot...
    ,m_dd,m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,g,g_cum,h,medie_mx,cov_parametri)
%MANCANO DA METTERE LE PRIOR ANCGE SE NON SERVONO A
%MOLTO!!!!!!!!!!!!XXXXXXXXXXXXXXXXXXX
% vettore in cui metto i pesi p(y|m_k): ? un riordinamento di g
p=zeros(1,N_iter);
% vettore dei logaritmi delle probabilit?
logp_new=zeros(1,N_iter);
% vettore che contiene i nuovi valori dei parametri
M_iperparametri_new_unnormalized=zeros(N_iter,2*J+2);

vec=1:(bigK-1);
for jj=1:N_iter
    U=unifrnd(0,1);
    ind=find(U<g_cum);
    ind=ind(1);
    M_iperparametri_new_unnormalized(jj,:)=mvnrnd(medie_mx(ind,:),h^2*cov_parametri);
    alpha=M_iperparametri_new_unnormalized(jj,1:J);
    d=M_iperparametri_new_unnormalized(jj,(J+1):(2*J));
    gamma=M_iperparametri_new_unnormalized(jj,2*J+1);
    nu=M_iperparametri_new_unnormalized(jj,2*J+2);
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
        vec_loggfc=zeros(1,bigK);
        for ii=1:bigK
            vec_loggfc(ii)=LogC(n_j_dot_k(j,ii)+1,mjk(j,ii)+1);
        end
        Phi_j=sum(log(alpha(j)+vec_j*d(j)))-gammaln(alpha(j)+nn(j))+gammaln(alpha(j)+1)+...
            sum(vec_loggfc)-m_j_dot(j)*log(d(j));
        Phi=Phi+Phi_j;
    end
    %qua calcola i pesi per la base measure P, conta i tavoli come uniche
    %osservazioni
    logp_new(jj)=Phi+sum(log(gamma+nu*vec))-gammaln(gamma+m_dd)+gammaln(gamma+1)+...
        sum(gammaln(m_dot_k-nu))-bigK*gammaln(1-nu);
    p(jj)=g(ind);
end
% calcolo i pesi: passo agli esponenziali e poi NORMALIZZO... Non faccio il
% ciclo ma lo calcolo in forma matriciale (questi sarebbero i pesi
% dell'importance sample, solo che poi ricampiono e quindi non ne tango
% traccia)
M=ones(N_iter,1)*logp_new;
p_new=1./sum(exp(M-M'),2)';
%qua viene calcolato il punto (5) dell'algoritmo generale di West
omega=p_new./p;
omega=omega/sum(omega);
