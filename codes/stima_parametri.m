% funzione che stima i parametri di uscita attraverso i tavoli
% generati con il metodo MCMC
function [mjk, m_dot_k, m_j_dot, m_dd alpha d gamma eta]=stima_parametri(M_l_star,tot_dist,M_parametri,J,N)

% stima di mjk: sarebbero quelli che io ho chiamato l_star
mjk=zeros(J,tot_dist);
for j=1:J
    for i=1:tot_dist
        mjk(j,i)=mean(M_l_star{j}(:,i));
    end
end
mjk=round(mjk);

% stima di m.k
m_dot_k=zeros(1,tot_dist);
for i=1:tot_dist
    m_dot_k(i)=sum(mjk(:,i));
end

% stima di mj.
m_j_dot=zeros(1,J);
for i=1:J
    m_j_dot(i)=sum(mjk(i,:));
end   

% stima di m..
m_dd=sum(m_j_dot);

% M parametri contiene tutti i parametri  da 1 a J alpha , da J+1 a 2*J d
% ed infine in 2J+1 e 2J+2 gamma ed eta
stima_iper=sum(M_parametri)/N;
alpha=stima_iper(1:J);
d=stima_iper((J+1):(2*J));
gamma=stima_iper(2*J+1);
eta=stima_iper(2*J+2);




