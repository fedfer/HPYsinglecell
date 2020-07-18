% estimate the parameters using the tables from MCMC


function [mjk, m_dot_k, m_j_dot, m_dd alpha d gamma eta]=estimate_parameters(M_l_star,tot_dist,M_parameters,J,N)

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


estimate_hyper=sum(M_parameters)/N;
alpha=estimate_hyper(1:J);
d=estimate_hyper((J+1):(2*J));
gamma=estimate_hyper(2*J+1);
eta=estimate_hyper(2*J+2);




