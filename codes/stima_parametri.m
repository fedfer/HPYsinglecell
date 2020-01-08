% estimate the parameters using the tables from MCMC


function [mjk, m_dot_k, m_j_dot, m_dd alpha d gamma eta]=stima_parametri(M_l_star,tot_dist,M_parametri,J,N)

% mjk
mjk=zeros(J,tot_dist);
for j=1:J
    for i=1:tot_dist
        mjk(j,i)=mean(M_l_star{j}(:,i));
    end
end
mjk=round(mjk);

%  m.k
m_dot_k=zeros(1,tot_dist);
for i=1:tot_dist
    m_dot_k(i)=sum(mjk(:,i));
end

%  mj.
m_j_dot=zeros(1,J);
for i=1:J
    m_j_dot(i)=sum(mjk(i,:));
end   

% m..
m_dd=sum(m_j_dot);


stima_iper=sum(M_parametri)/N;
alpha=stima_iper(1:J);
d=stima_iper((J+1):(2*J));
gamma=stima_iper(2*J+1);
eta=stima_iper(2*J+2);




