% partially exchangeable case

% T= total iteration and t_burn= burn in

function [M_Tavoli, M_l_star, M_parametri, Dati_star, k_popolazioni]=posterior_K(Dati,M0,V0,J,n_init,T,t_burn)


N=T-t_burn;

% generate initial hyperparameters
theta_0=1;
theta=ones(1,J);
sigma_0=0.5;
sigma=ones(1,J)/2;


Tavoli=cell(1,J);
Tavoli_star=cell(1,J);
M_Tavoli=cell(1,J);
Dati_star=cell(1,J);
TavoliDati_star=cell(1,J);
n_star=cell(1,J);
l_star=cell(1,J);
q_star=cell(1,J);
L=zeros(1,J);
Dati_totali=[];
k_popolazioni=zeros(1,J);
for j=1:J
    Tavoli{j}=normrnd(sum(Dati{j})/n_init(j),1,1,n_init(j));
    M_Tavoli{j}=zeros(N,n_init(j));
    % cluster the observations and occurencies
    Dati_star{j}=clusterizza(Dati{j});
    k_popolazioni(j)=length(Dati_star{j});
    Dati_totali=[Dati_totali Dati{j}];
end
% distinct species in all pop
Dati_totali_star=unique(Dati_totali);

k=length(Dati_totali_star);

M_parametri=zeros(N,2*J+2);
M_Tavoli=cell(1,J);
% numerosity of tables that eat a certani dish 
M_l_star=cell(1,J);
for j=1:J
    M_Tavoli{j}=zeros(N,n_init(j));
    M_l_star{j}=zeros(N,k);
end


for ii=1:T

    for j=1:J
        Tavoli{j}=update_tables_n(j,J,Dati,Tavoli,theta_0,sigma_0,theta(j),sigma(j),n_init,M0,V0);
    end

    % update hyperparameters 4 MH
    l_star_franchise=zeros(1,k);
    for j=1:J
        [Dati_star{j}, Tavoli_star{j}, TavoliDati_star{j}, n_star{j}, q_star{j}, l_star{j}] = distinct_fct(Dati{j}, Tavoli{j});
        L(j)=sum(l_star{j});
        I=(ismember(Dati_totali_star,Dati_star{j}));
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
            M_Tavoli{j}(ii-t_burn,:)=Tavoli{j};
        end
        M_parametri(ii-t_burn,:)=[theta sigma theta_0 sigma_0];
    end
end
    


