% function to compute the normalized weights of particle filter

    
function [M_hyper_new_unnormalized omega]=weight_filter_new(mjk,m_j_dot...
    ,m_dd,m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,g,g_cum,h,mean_mx,cov_parameters)

% weights of p(y|m_k)
p=zeros(1,N_iter);
% log of probabilities
logp_new=zeros(1,N_iter);

M_hyper_new_unnormalized=zeros(N_iter,2*J+2);

vec=1:(bigK-1);
for jj=1:N_iter
    U=unifrnd(0,1);
    ind=find(U<g_cum);
    ind=ind(1);
    % step 2 of Algorithm 2
    M_hyper_new_unnormalized(jj,:)=mvnrnd(mean_mx(ind,:),h^2*cov_parameters);
    alpha=M_hyper_new_unnormalized(jj,1:J);
    d=M_hyper_new_unnormalized(jj,(J+1):(2*J));
    gamma=M_hyper_new_unnormalized(jj,2*J+1);
    nu=M_hyper_new_unnormalized(jj,2*J+2);
    Phi=0;
    for j=1:J
        %  EPPF for restaurant j (explog)
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
