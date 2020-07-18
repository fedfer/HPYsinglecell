% Update the hyperparameters using particle filter of Liu and West

function [ alpha, d ,gamma ,nu ,M_hyper_new] = Filter_hyperparameters(...
    mjk, m_j_dot, m_dd, m_dot_k, n_j_dot_k, J, nn, bigK, N_iter, M_hyperparameters_old)

h=1/N_iter;
a=sqrt(1-h^2);
E_parameters=mean(M_hyperparameters_old);
cov_parameters=cov(M_hyperparameters_old);
mean_mx=a*M_hyperparameters_old+(1-a)*ones(N_iter,1)*E_parameters;
g=weight_filter(mjk,m_j_dot,m_dd,m_dot_k,n_j_dot_k,J,nn,mean_mx,bigK,N_iter);
g_cum=cumsum(g);
M_hyper_new=zeros(N_iter,2*J+2);
% compute the vector of new parameters and weights
% step (2) di pag 11 West
[M_hyper_new_unnormalized, omega]=weight_filter_new(mjk,m_j_dot,m_dd,...
    m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,g,g_cum,h,mean_mx,cov_parameters);
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
