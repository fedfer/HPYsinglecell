% Update the hyperparameters using particle filter of Liu and West
function [ alpha, d ,gamma ,nu ,M_hyper_new]=Filter_hyperparameters_v1(...
    mjk,m_j_dot,m_dd,m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,M_hyper_old,...
    mjk_old, m_j_dot_old,m_dd_old,m_dot_k_old,...
    nj_dot_k_old,nn_old,bigK_old)


h=1/N_iter;
a=sqrt(1-h^2);
E_parameters=mean(M_hyper_old);
cov_parameters=cov(M_hyper_old);
mean_mx=a*M_hyper_old+(1-a)*ones(N_iter,1)*E_parameters;

% same weight_filter_filter function as before, just need to compute the other
% weights as well
g=weight_filter_(mjk,m_j_dot,m_dd,m_dot_k,n_j_dot_k,J,nn,mean_mx,bigK,N_iter);
g_prec=weight_filter(mjk_old,m_j_dot_old,m_dd_old,m_dot_k_old,...
    nj_dot_k_old,J,nn_old,mean_mx,bigK_old,N_iter);

g = g ./ g_prec;
g_cum=cumsum(g);

M_hyper_new=zeros(N_iter,2*J+2);
% compute the vector of new parameters and weights
% step (2) di pag 11 West
[M_hyper_new_unnormalized, omega]=weight_filter_new_v1(mjk,m_j_dot,m_dd,...
    m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,g,g_cum,h,mean_mx,cov_parameters,...
    mjk_old,m_j_dot_old,m_dd_old,m_dot_k_old,nj_dot_k_old,nn_old,bigK_old);
% resample to have the weights equal to 1/N
omega_cum=cumsum(omega);
for ii=1:N_iter
    U=unifrnd(0,1);
    ind=find(U<omega_cum);
    ind=ind(1);
    M_hyper_new(ii,:)=M_hyper_new_unnormalized(ind,:);
end

% Find the parameters estimates, as the mean
alpha=mean(M_hyper_new(:,1:J));
d=mean(M_hyper_new(:,(J+1):(2*J)));
gamma=mean(M_hyper_new(:,2*J+1));
nu=mean(M_hyper_new(:,2*J+2));
