% caso parzialmente scambiabile: prende in ingresso il cell array Dati che 
% contiene i dati delle varie popolazioni. Voglio applicare un metodo MCMC
% per trovarmi le stime dei parametri riferiti ai tavoli
% T= iterazioni totali e t_burn= numero di burn in

function [M_Tavoli, M_l_star, M_parametri, Dati_star, k_popolazioni]=posterior_K(Dati,M0,V0,J,n_init,T,t_burn)

% M0 e V0 sono media e varianza della misura base che assumiamo gaussiana
% T ? il numero di iterazioni
N=T-t_burn;

% genero gli iperparametri iniziali
theta_0=1;
theta=ones(1,J);
sigma_0=0.5;
sigma=ones(1,J)/2;

% fisso a caso i valori iniziali di T
Tavoli=cell(1,J);
Tavoli_star=cell(1,J);
M_Tavoli=cell(1,J);
Dati_star=cell(1,J);
TavoliDati_star=cell(1,J);
n_star=cell(1,J);
l_star=cell(1,J);
q_star=cell(1,J);
% vettore contenete le numerosit? dei tavoli
L=zeros(1,J);
Dati_totali=[];
% numero di specie in ogni popolazione
k_popolazioni=zeros(1,J);
for j=1:J
    Tavoli{j}=normrnd(sum(Dati{j})/n_init(j),1,1,n_init(j));
    M_Tavoli{j}=zeros(N,n_init(j));
    %la funzione clusterizza calcola i dati distinti e le occorrenze
    Dati_star{j}=clusterizza(Dati{j});
    k_popolazioni(j)=length(Dati_star{j});
    Dati_totali=[Dati_totali Dati{j}];
end
% specie distinte in tutte le popolazioni
Dati_totali_star=unique(Dati_totali);
% numero di specie totali
k=length(Dati_totali_star);
% vettori per memorizzare i valori di uscita: prima menorizzo il vettore
% thea, poi sigma, poi theat_0 ed infine sigma_0
M_parametri=zeros(N,2*J+2);
% ogni elemento del cell array ? una matrice
M_Tavoli=cell(1,J);
% qui inserisco le numerosit? dei tavoli in cui viene mangiato il piatto k
M_l_star=cell(1,J);
for j=1:J
    M_Tavoli{j}=zeros(N,n_init(j));
    M_l_star{j}=zeros(N,k);
end


for ii=1:T
    % aggiorno i tavoli da 1,...,n1 per Xn1 e i tavoli da 1, ..., n2 per
    % Xn2 ecc
    for j=1:J
        Tavoli{j}=aggiorna_tavoli_n(j,J,Dati,Tavoli,theta_0,sigma_0,theta(j),sigma(j),n_init,M0,V0);
    end
    % non devo fare passi di accelerazione perch? non mik interessa il
    % valore delle osservazioni

    % aggiorno gli iperparametri dalle full conditionals secondo 4 MH
    % calcolo \overline{l} per tutti i franchise
    l_star_franchise=zeros(1,k);
    for j=1:J
        [Dati_star{j}, Tavoli_star{j}, TavoliDati_star{j}, n_star{j}, q_star{j}, l_star{j}] = distinti(Dati{j}, Tavoli{j});
        L(j)=sum(l_star{j});
        I=(ismember(Dati_totali_star,Dati_star{j}));
        l_star_franchise(I)=l_star_franchise(I)+l_star{j};
        if ii>t_burn
            % trovo gli indici di dati_star{j} in Dati_totali_star: entrambi i
            % vettori devono essere ordinati!!
            M_l_star{j}(ii-t_burn,I)=l_star{j};
        end
        
        % aggiornamento dei parametri
        theta(j)=MH_theta(theta(j),sigma(j),L(j),n_init(j));
        [sigma(j)]=MH_sigma(sigma(j),theta(j),q_star{j});
    end
    
    theta_0=MH_theta(theta_0,sigma_0,sum(L),k);
    % cell2mat trasforma un cell array in un vettore normale
	[sigma_0]=MH_sigma(sigma_0,theta_0,l_star_franchise);
    
    if ii>t_burn
        for j=1:J
            M_Tavoli{j}(ii-t_burn,:)=Tavoli{j};
        end
        M_parametri(ii-t_burn,:)=[theta sigma theta_0 sigma_0];
    end
end
    


