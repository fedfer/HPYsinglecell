% function to compute the normalized weights of particle filter

    
function [M_iperparametri_new_unnormalized, omega]=pesi_filter_new_v1(mjk,m_j_dot...
    ,m_dd,m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,g,g_cum,h,medie_mx,cov_parametri,...
    mjk_old,m_j_dot_old,m_dd_old,m_dot_k_old,nj_dot_k_old,nn_old,bigK_old)

% weights of p(y|m_k)
p=zeros(1,N_iter);
% log of probabilities
logp_new=zeros(1,N_iter);

M_iperparametri_new_unnormalized=zeros(N_iter,2*J+2);

vec=1:(bigK-1);
vec_old=1:(bigK_old-1);

for jj=1:N_iter
    U=unifrnd(0,1);
    ind=find(U<g_cum);
    ind=ind(1);
    % step 2 of Algorithm 2
    M_iperparametri_new_unnormalized(jj,:)=mvnrnd(medie_mx(ind,:),h^2*cov_parametri);
    alpha=M_iperparametri_new_unnormalized(jj,1:J);
    d=M_iperparametri_new_unnormalized(jj,(J+1):(2*J));
    gamma=M_iperparametri_new_unnormalized(jj,2*J+1);
    nu=M_iperparametri_new_unnormalized(jj,2*J+2);
    Phi=0;
    for j=1:J
        %  EPPF for restaurant j (explog)
        vec_j=1:(m_j_dot(j)-1);
        vec_j_old=1:(m_j_dot_old(j)-1);
        % find max of n_j.k for generalized factorial coeffs
        max_nj=max(n_j_dot_k(j,:));
        max_mj=max(mjk(j,:));
        max_nj_old=max(nj_dot_k_old(j,:));
        max_mj_old=max(mjk_old(j,:));

        LogC=generalized_factorial(max_nj,max_mj,d(j));
        LogC_old=generalized_factorial(max_nj_old,max_mj_old,d(j));
        % likelihood
        vec_loggfc=zeros(1,bigK);
        vec_loggfc_old=zeros(1,bigK_old);
        for ii=1:bigK
            vec_loggfc(ii)=LogC(n_j_dot_k(j,ii)+1,mjk(j,ii)+1);
        end
        for ii=1:bigK_old
            vec_loggfc_old(ii)=LogC_old(nj_dot_k_old(j,ii)+1,mjk_old(j,ii)+1);
        end
        Phi_j=sum(log(alpha(j)+vec_j*d(j)))-gammaln(alpha(j)+nn(j))+gammaln(alpha(j)+1)+...
            sum(vec_loggfc)-m_j_dot(j)*log(d(j));
        
        Phi_j_old=sum(log(alpha(j)+vec_j_old*d(j)))-gammaln(alpha(j)+nn_old(j))+gammaln(alpha(j)+1)+...
            sum(vec_loggfc_old)-m_j_dot_old(j)*log(d(j));
        
        Phi=Phi + Phi_j - Phi_j_old;
    end

    logp_new(jj)=Phi+sum(log(gamma+nu*vec))-gammaln(gamma+m_dd)+gammaln(gamma+1)+...
        sum(gammaln(m_dot_k-nu))-bigK*gammaln(1-nu) - ...
        (sum(log(gamma+nu*vec_old))-gammaln(gamma+m_dd_old)+gammaln(gamma+1)+...
        sum(gammaln(m_dot_k_old-nu))-bigK_old*gammaln(1-nu));
    
    p(jj)=g(ind);
end

M=ones(N_iter,1)*logp_new;
p_new=1./sum(exp(M-M'),2)';
omega=p_new./p;
omega=omega/sum(omega);