% Update the hyperparameters using particle filter of Liu and West

function [ alpha, d ,gamma ,nu ,M_iperparametri_new]=Filter_iperparametri(...
    mjk,m_j_dot,m_dd,m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,M_iperparametri_old)

h=1/N_iter;
a=sqrt(1-h^2);
E_parametri=mean(M_iperparametri_old);
cov_parametri=cov(M_iperparametri_old);
medie_mx=a*M_iperparametri_old+(1-a)*ones(N_iter,1)*E_parametri;
g=pesi_filter(mjk,m_j_dot,m_dd,m_dot_k,n_j_dot_k,J,nn,medie_mx,bigK,N_iter);
g_cum=cumsum(g);
M_iperparametri_new=zeros(N_iter,2*J+2);
% compute the vector of new parameters and weights
% step (2) di pag 11 West
[M_iperparametri_new_unnormalized, omega]=pesi_filter_new(mjk,m_j_dot,m_dd,...
    m_dot_k,n_j_dot_k,J,nn,bigK,N_iter,g,g_cum,h,medie_mx,cov_parametri);
% resample to have the weights equal to 1/N
omega_cum=cumsum(omega);
for ii=1:N_iter
    U=unifrnd(0,1);
    ind=find(U<omega_cum);
    ind=ind(1);
    M_iperparametri_new(ii,:)=M_iperparametri_new_unnormalized(ind,:);
end

% Find the parameters estimates, as the mean
alpha=mean(M_iperparametri_new(:,1:J));
d=mean(M_iperparametri_new(:,(J+1):(2*J)));
gamma=mean(M_iperparametri_new(:,2*J+1));
nu=mean(M_iperparametri_new(:,2*J+2));
