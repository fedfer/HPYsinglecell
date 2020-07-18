% partially exchangeable case
% input: vector that contains the observations from each population 
% output: estimates for tables using MCMC


function [M_Tables, M_l_star, M_parameters, Data_star, k_populations]=chinese_restaurant(Data,M0,V0,J,n_init,T,t_burn)

N=T-t_burn;

% generate initial parameters
gamma=1;
alpha=ones(1,J);
nu=0.5;
d=ones(1,J)/2;

Tables=cell(1,J);
Tables_star=cell(1,J);
M_Tables=cell(1,J);
Data_star=cell(1,J);
TablesData_star=cell(1,J);
n_star=cell(1,J);
l_star=cell(1,J);
q_star=cell(1,J);

L=zeros(1,J);
Data_total=[];
% number of species in each pop
k_populations=zeros(1,J);
for j=1:J
    Tables{j}=normrnd(sum(Data{j})/n_init(j),1,1,n_init(j));
    M_Tables{j}=zeros(N,n_init(j));
    Data_star{j}=cluster_fct(Data{j});
    k_populations(j)=length(Data_star{j});
    Data_total=[Data_total Data{j}];
end
% distinct species in all pop
Data_total_star=unique(Data_total);
k=length(Data_total_star);

M_parameters=zeros(N,2*J+2);
M_Tables=cell(1,J);
M_l_star=cell(1,J);
for j=1:J
    M_Tables{j}=zeros(N,n_init(j));
    M_l_star{j}=zeros(N,k);
end

mjk=zeros(J,k);
m_dot_k=zeros(1,k);
m_j_dot=zeros(1,J);

for ii=1:T

    for j=1:J
        Tables{j}=update_tables_n(j,J,Data,Tables,gamma,nu,alpha(j),d(j),n_init,M0,V0);
    end


    % update the parametes using 4 MH
    for j=1:J
        [Data_star{j}, Tables_star{j}, TablesData_star{j}, n_star{j}, q_star{j}, l_star{j}] = distinct_fct(Data{j}, Tables{j});
        L(j)=sum(l_star{j});
        if ii>t_burn

            I=(ismember(Data_total_star,Data_star{j}));
            M_l_star{j}(ii-t_burn,I)=l_star{j};
        else
            I=(ismember(Data_total_star,Data_star{j}));
        end
        mjk(j,I)=l_star{j};
        m_dot_k=sum(mjk);
        m_dd=sum(m_dot_k);
        m_j_dot=sum(mjk');
        
    end
    
[ alpha, d ,gamma ,nu]=MCMC_hyper(mjk,m_j_dot,m_dd,J,n_init,2,0,alpha,d,gamma,nu,k);
    
    if ii>t_burn
        for j=1:J
            M_Tables{j}(ii-t_burn,:)=Tables{j};
        end
        M_parameters(ii-t_burn,:)=[alpha d gamma nu];
    end
end
    


