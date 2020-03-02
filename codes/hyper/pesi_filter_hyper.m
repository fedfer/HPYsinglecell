% update the normalized weights of particle filter
function p=pesi_filter_hyper(mjk,m_j_dot,m_dd,m_dot_k,n_j_dot_k,J,nn,medie_mx,bigK,N_iter)

vec=1:(bigK-1);
% vector of log probabilities 
%N_iter number of elements in sample
logp=zeros(1,N_iter);
for jj=1:N_iter
    alpha=medie_mx(jj,1:J);
    d=medie_mx(jj,(J+1):(2*J));
    gamma=medie_mx(jj,2*J+1);
    nu=medie_mx(jj,2*J+2);
    theta_h=medie_mx(jj,2*J+2+1);
    a0=medie_mx(jj,2*J+2+2);
    b0=medie_mx(jj,2*J+2+3);
Phi=0;
for j=1:J
    % compute the EPPF for jth restaurant (explog)
    vec_j=1:(m_j_dot(j)-1);
    % max of n_j.k for generalized factorial coefficient
    max_nj=max(n_j_dot_k(j,:));
    max_mj=max(mjk(j,:));

    % generate the full matrix
    LogC=generalized_factorial(max_nj,max_mj,d(j));
    % compute for likelihood
    % big K is the number of unique values
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

M=ones(N_iter,1)*logp;

p=1./sum(exp(M-M'),2)';
