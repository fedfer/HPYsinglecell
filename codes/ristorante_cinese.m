% partially exchangeable case
% input: vector that contains the observations from each population 
% output: estimates for tables using MCMC


function [M_Tavoli, M_l_star, M_parametri, Dati_star, k_popolazioni]=ristorante_cinese(Dati,M0,V0,J,n_init,T,t_burn)

N=T-t_burn;

% generate initial parameters
gamma=1;
alpha=ones(1,J);
nu=0.5;
d=ones(1,J)/2;

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
% number of species in each pop
k_popolazioni=zeros(1,J);
for j=1:J
    Tavoli{j}=normrnd(sum(Dati{j})/n_init(j),1,1,n_init(j));
    M_Tavoli{j}=zeros(N,n_init(j));
    Dati_star{j}=cluster_fct(Dati{j});
    k_popolazioni(j)=length(Dati_star{j});
    Dati_totali=[Dati_totali Dati{j}];
end
% distinct species in all pop
Dati_totali_star=unique(Dati_totali);
k=length(Dati_totali_star);

M_parametri=zeros(N,2*J+2);
M_Tavoli=cell(1,J);
M_l_star=cell(1,J);
for j=1:J
    M_Tavoli{j}=zeros(N,n_init(j));
    M_l_star{j}=zeros(N,k);
end

mjk=zeros(J,k);
m_dot_k=zeros(1,k);
m_j_dot=zeros(1,J);

for ii=1:T

    for j=1:J
        Tavoli{j}=update_tables_n(j,J,Dati,Tavoli,gamma,nu,alpha(j),d(j),n_init,M0,V0);
    end


    % update the parametes using 4 MH
    for j=1:J
        [Dati_star{j}, Tavoli_star{j}, TavoliDati_star{j}, n_star{j}, q_star{j}, l_star{j}] = distinct_fct(Dati{j}, Tavoli{j});
        L(j)=sum(l_star{j});
        if ii>t_burn

            I=(ismember(Dati_totali_star,Dati_star{j}));
            M_l_star{j}(ii-t_burn,I)=l_star{j};
        else
            I=(ismember(Dati_totali_star,Dati_star{j}));
        end
        mjk(j,I)=l_star{j};
        m_dot_k=sum(mjk);
        m_dd=sum(m_dot_k);
        m_j_dot=sum(mjk');
        
    end
    
[ alpha, d ,gamma ,nu]=MCMC_hyper(mjk,m_j_dot,m_dd,J,n_init,2,0,alpha,d,gamma,nu,k);
    
    if ii>t_burn
        for j=1:J
            M_Tavoli{j}(ii-t_burn,:)=Tavoli{j};
        end
        M_parametri(ii-t_burn,:)=[alpha d gamma nu];
    end
end
    


