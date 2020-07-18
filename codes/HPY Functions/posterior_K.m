% partially exchangeable case
% T= total iteration and t_burn= burn in

function [M_Tables, M_l_star, M_parameters, Data_star, k_populations]=posterior_K(Data,M0,V0,J,n_init,T,t_burn)


N=T-t_burn;

% generate initial hyperparameters
theta_0=1;
theta=ones(1,J);
sigma_0=0.5;
sigma=ones(1,J)/2;


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
k_populations=zeros(1,J);
for j=1:J
    Tables{j}=normrnd(sum(Data{j})/n_init(j),1,1,n_init(j));
    M_Tables{j}=zeros(N,n_init(j));
    % cluster the observations and occurencies
    Data_star{j}=clusterizza(Data{j});
    k_populations(j)=length(Data_star{j});
    Data_total=[Data_total Data{j}];
end
% distinct species in all pop
Data_total_star=unique(Data_total);

k=length(Data_total_star);

M_parameters=zeros(N,2*J+2);
M_Tables=cell(1,J);
% numerosity of tables that eat a certani dish 
M_l_star=cell(1,J);
for j=1:J
    M_Tables{j}=zeros(N,n_init(j));
    M_l_star{j}=zeros(N,k);
end


for ii=1:T

    for j=1:J
        Tables{j}=update_tables_n(j,J,Data,Tables,theta_0,sigma_0,theta(j),sigma(j),n_init,M0,V0);
    end

    % update hyperparameters 4 MH
    l_star_franchise=zeros(1,k);
    for j=1:J
        [Data_star{j}, Tables_star{j}, TablesData_star{j}, n_star{j}, q_star{j}, l_star{j}] = distinct_fct(Data{j}, Tables{j});
        L(j)=sum(l_star{j});
        I=(ismember(Data_total_star,Data_star{j}));
        l_star_franchise(I)=l_star_franchise(I)+l_star{j};
        if ii>t_burn
            M_l_star{j}(ii-t_burn,I)=l_star{j};
        end
        
        % update parameters
        theta(j)=MH_theta(theta(j),sigma(j),L(j),n_init(j));
        [sigma(j)]=MH_sigma(sigma(j),theta(j),q_star{j});
    end
    
    theta_0=MH_theta(theta_0,sigma_0,sum(L),k);

	[sigma_0]=MH_sigma(sigma_0,theta_0,l_star_franchise);
    
    if ii>t_burn
        for j=1:J
            M_Tables{j}(ii-t_burn,:)=Tables{j};
        end
        M_parameters(ii-t_burn,:)=[theta sigma theta_0 sigma_0];
    end
end
    


