% function to compute the normalized weights of particle filter

    
function [M_iperparametri_new_unnormalized omega]=pesi_filter_new_hyper(mjk,m_j_dot...
    ,m_dd,m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,g,g_cum,h,medie_mx,cov_parametri)

% weights of p(y|m_k)
p=zeros(1,N_iter);
% log of probabilities
logp_new=zeros(1,N_iter);

M_iperparametri_new_unnormalized=zeros(N_iter,2*J+2+3);

vec=1:(bigK-1);
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
    theta_h=M_iperparametri_new_unnormalized(jj,2*J+2+1);
    a0=M_iperparametri_new_unnormalized(jj,2*J+2+2);
    b0=M_iperparametri_new_unnormalized(jj,2*J+2+3);
    Phi=0;
    for j=1:J
        %  l'EPPF for restaurant j (explog)
        vec_j=1:(m_j_dot(j)-1);
        % find max of n_j.k for generalized factorial coeffs
        max_nj=max(n_j_dot_k(j,:));
        max_mj=max(mjk(j,:));

        LogC=generalized_factorial(max_nj,max_mj,d(j));
        % likelihood
        vec_loggfc=zeros(1,bigK);
        for ii=1:bigK
            vec_loggfc(ii)=LogC(n_j_dot_k(j,ii)+1,mjk(j,ii)+1);
        end
        Phi_j=sum(log(alpha(j)+vec_j*d(j)))-gammaln(alpha(j)+nn(j))+gammaln(alpha(j)+1)+...
            sum(vec_loggfc)-m_j_dot(j)*log(d(j));
        Phi=Phi+Phi_j;
    end

    logp_new(jj)=Phi+sum(log(gamma+nu*vec))-gammaln(gamma+m_dd)+gammaln(gamma+1)+...
        sum(gammaln(m_dot_k-nu))-bigK*gammaln(1-nu);
    p(jj)=g(ind);
end

M=ones(N_iter,1)*logp_new;
p_new=1./sum(exp(M-M'),2)';
omega=p_new./p;
omega=omega/sum(omega);