%clear all;
close all;

% numero di iterazioni
Runs=60;

%% Insieme dei parametri

% numero di iterazioni in MCMC per il numero di tavoli e dei parametri di
% HPY dato il campione iniziale
iter=35000;
burnin=15000;

% numero di iterazioni per l'aggiornamento dei parametri con particle
% filter
N_iter=1000;

% iterazioni in MCMC per gli hyper HPY dopo ogni campionamento
% iter1=200;
% burnin1=100;

% parametro di Good-Turing
C=(1+sqrt(2))*sqrt(3);

% scarico i dati dei microbi dai file .csv con il comando csvread leggo il
% file dalla riga 1 e colonna 1. Le colonne rappresentano le varie
% popolazioni, mentre le righe rappresentano le specie
freq_GB_pop=csvread('GlobalPatterns_data.csv',1,1);
% numero delle specie
k=size(freq_GB_pop,1);

% dal file .sample ricavo le superpopolazioni e quindi gli estremi delle
% colonne da sommare
col_super_GB=[1 4 6 9 11 13 16 19 22 24 27];
freq1=zeros(k,length(col_super_GB)-1);
for ii=1:(length(col_super_GB)-1)
    for jj=col_super_GB(ii):(col_super_GB(ii+1)-1)
    freq1(:,ii)=freq_GB_pop(:,jj)+freq1(:,ii);
    end
end

% le righe della matrice rappresentano le specie quindi i loro labels sono
% numeri 
% numero totale delle specie presenti in tutte le popolazioni
N=k;
labels=1:N;
% numero delle popolazioni
J=size(freq1,2);

%frequenze nella popolazione
freq=cell(1,J);
% labels delle specie nelle varie popolazioni
pop=cell(1,J);
% numero delle specie distinte nella popolazione i esima
NN=zeros(1,J);
for i=1:J
    % trovo le frequenze empiriche normalizzate relative alla popolazione
    % i-esima, togliendo quelle nulle
    freq{i}=freq1(:,i)/sum(freq1(:,i));
    freq{i}=freq{i}(freq1(:,i)~=0);
    % mettto in modo corrispondente i labels
    pop{i}=1:N;
    pop{i}=pop{i}(freq1(:,i)~=0);
    NN(i)=length(pop{i});
end


% ampiezza del campione iniziale
n_init=30*ones(J,1);
% lunghezza del campione addizionale
addsample=300;

% Dati finali
M=zeros(Runs,addsample);
DATAfinal=struct('HPY',M,'uniform',M,'Oracle',M,'GoodTuring',M);

%weights
M=zeros(Runs,J);
WEIGTHS=struct('weight_HPY',M,'weight_uniform',M,'weight_Oracle',M,'weight_GoodTuring',M);
clear M;
    
   %genera i dati iniziali per le varie popolazioni
   % data{j}= vettore che contiene il campione della popolazione j
   data=cell(1,J);
   for i=1:J 
       data{i}=zeros(1,n_init(i));
   end
   for j=1:J
       for i=1:n_init(j)
           data{j}(i)=gendiscr(pop{j},freq{j});
       end
   end
   dati_totali=cell2mat(data);
   % specie dsitinte di tutta la popolazione iniziale
   Kini=unique(dati_totali);
   %numero di specie iniziali distinte
   tot_dist=length(Kini);
   
   %% Inferenza MCMC sul numero di specie e sul numero di parametri
   % questi parametri sono quelli riferiti alle prior degli hyperparametri
   % del PY
   M0=1;
   V0=4;
   bigK=length(Kini);
   
   % aggiornamento con l'algoritmo marginale che ho sviluppato con Antonio
   % ed Igor
   [M_Tavoli M_l_star M_parametri Dati_star k_popolazioni]=posterior_K(data,M0,V0,J,n_init,iter,burnin);
  % save ristorante_cinese  M_Tavoli M_l_star M_parametri Dati_star k_popolazioni Kini tot_dist dati_totali data pop DATAfinal WEIGTHS
   
  % % Aggiornamento con l'algoritmo di Marco: cambia come si aggiornano gli
  % % iperparametri
  % [M_Tavoli, M_l_star, M_parametri, Dati_star, k_popolazioni]=ristorante_cinese(data,M0,V0,J,n_init,iter,burnin);
   
   % ora stimo i parametri che mi servono dopo
   [mjk_ini, m_dot_k_ini, m_j_dot_ini, m_dd_ini alpha d gamma nu]=stima_parametri(M_l_star,tot_dist,M_parametri,J,iter-burnin);
   
      
   % creo un vettore di pesi per la mistura delle normali dei parametri: ne
% prendo solo una parte.
M_parametri_ini=M_parametri(end-N_iter+1:end,:);
   
    %% Algoritmi vari per scegliere da dove campionare la prossima
    % osservazione: Unif, HPY Oracle GT
    
for III=1:Runs
    % Uniform sampling
    supp=1:J;
    masses=ones(1,J)/J;
    KuniUni=Kini;
    bigK=length(Kini);
    
    %HPY
    nn=n_init;
    KuniHPY=Kini;
    mjk=mjk_ini;
    m_dot_k=m_dot_k_ini;
    m_j_dot=m_j_dot_ini;
    m_dd=m_dd_ini;
    M_parametri=M_parametri_ini;
    
    % calcolo n.k
    nj_dot_k=zeros(J,tot_dist);
    for j=1:J
        for w=1:tot_dist
            nj_dot_k(j,w) = sum(data{j}==KuniHPY(w));
        end
    end
    
    
    % oracle
    KuniOracle=Kini;
    missingmass=zeros(1,J);
    for j=1:J
        labelsseen=ismember(pop{j},KuniOracle);
        % ritorna un vettore di 1 e 0: ho 1 se pop{j} è elemento di
        % KuniOracle
        missingmass(j) = 1-sum(freq{j}(labelsseen));
    end
    
    % Good--turing
    
   KuniGT=Kini;
   nGT=n_init;
   % trovo il numero di osservazioni con frequenza 1 nel campione congiunto
   [data_uni freq_uni]=clusterizza(cell2mat(data));
   unGT_tot=data_uni(freq_uni==1);
   clear data_uni freq_uni;
   % trovo le specie con frequenza uno in ogni popolazione
   unGT=cell(1,J);
   for j=1:J
       indici=ismember(data{j},unGT_tot);
       unGT{j}=unique(data{j}(indici));
   end
   clear data_uni freq_uni;
   
   
   % indicatori delle scoperte di nuove specie
   newobsindHPY=zeros(1,addsample);
   newobsindOra=zeros(1,addsample);
   newobsindUni=zeros(1,addsample);
   newobsindGT=zeros(1,addsample);
   
   % indicatori per i pesi
   w_HPY=zeros(1,J);
   w_Ora=zeros(1,J);
   w_Uni=zeros(1,J);
   w_GT=zeros(1,J);
   
   % spezio per contenere i risultati relativi alla stima degli
   % iperparametri
   
   
   %% Inizio degli algoritmi 
   

   
   for i=1:addsample
       
       % scegliamo la popolazione per la strategia Uni
       unif=gendiscr(supp,masses);
       w_Uni(unif)=w_Uni(unif)+1;
       
       
       % scegliamo la popolazone per la strategia HPY
       betadraws=zeros(1,J);
       betazero=betarnd(gamma+nu*bigK,m_dd-bigK*nu);
       for j=1:J
           betadraws(j)=betarnd(betazero*(alpha(j)+(m_j_dot(j)*d(j))), ((1-betazero)*(alpha(j)+m_j_dot(j)*d(j))+ nn(j)- m_j_dot(j)*d(j)));
       end
       [v_max, armchosen]= max(betadraws);
       clear v_max;
       w_HPY(armchosen)=w_HPY(armchosen)+1;
       
       % scegliamo la popolazione per la strategia Oracle
       [v_max armchosenOrac]=max(missingmass);
       clear v_max;
       if length(armchosenOrac)>1
           Oracprob=ones(1,length(armchosenOrac))/armchosenOrac;
           armchosenOrac=gendiscr(armchosenOrac,Oracprob);
       end
       
       w_Ora(armchosenOrac)=w_Ora(armchosenOrac)+1;
       
       % scegliamo la popolazione per la strategia GT
       R_GT=zeros(1,J);
       for j=1:J
           R_GT(j)=(length(unGT{j})/nGT(j))+C*sqrt(log(4*sum(nGT))/nGT(j));
       end
       [v_max armchosenGT]=max(R_GT);
       clear v_max;
       if length(armchosenGT)>1
           GTprob=ones(1,length(armchosenGT))/armchosenGT;
           armchosenGT=gendiscr(armchosenGT,GTprob);
       end
       
       w_GT(armchosenGT)=w_GT(armchosenGT)+1;
       
       % campioniamo un nuovo valore per ogni popolazione scelta
       
       newobservations=zeros(1,J);
       for j=1:J
           newobservations(j) = gendiscr(pop{j},freq{j});
       end
       newunif=newobservations(unif);
       newobsHPY=newobservations(armchosen);
       newobsOrac=newobservations(armchosenOrac);
       newobsGT=newobservations(armchosenGT);
       
       % aggiornamento dell'uniforme
       if sum(KuniUni==newunif)==0
           newobsindUni(i)=1;
           KuniUni=[KuniUni newunif];
       end
       
       % aggiornamento dei tavoli per HPY
       if sum(KuniHPY==newobsHPY)==0
           bigK=bigK+1;
           m_j_dot(armchosen)=m_j_dot(armchosen)+1;
           m_dd=m_dd+1;
           KuniHPY=[KuniHPY newobsHPY];
           newobsindHPY(i)=1;
           m_dot_k=[m_dot_k 1];
           nj_dot_k=[nj_dot_k , zeros(J,1)];
           mjk=[mjk , zeros(J,1)];
           nj_dot_k(armchosen,bigK)=1;
           mjk(armchosen,bigK)=1;
       else
           
           % calcolo la probabilità di sedersi a un tavolo vecchio se
           % l'osservazione è vecchia
           olddistinct=find(KuniHPY==newobsHPY);
           probnewold=[nj_dot_k(armchosen,olddistinct)-d(armchosen)*mjk(armchosen,olddistinct),...
               (alpha(armchosen)+m_j_dot(armchosen)*d(armchosen))*((m_dot_k(olddistinct)-nu)/(gamma+m_dd))];
           probnewold=probnewold/sum(probnewold);
           bern=binornd(1,probnewold(1));
           if bern==0
               % l'osservazione forma un nuovo tavolo
               m_j_dot(armchosen)=m_j_dot(armchosen)+1;
               m_dot_k(olddistinct)= m_dot_k(olddistinct)+1;
               mjk(armchosen,olddistinct)=mjk(armchosen,olddistinct)+1;
               m_dd=m_dd+1;
           end
           nj_dot_k(armchosen,olddistinct)=nj_dot_k(armchosen,olddistinct)+1;
       end
       nn(armchosen)=nn(armchosen)+1;
       
       % aggiorno gli iperparamegtri del HPY utilizzando un MH
     %  [ alpha, d ,gamma ,nu]=MCMC_iperparametri(mjk,m_j_dot,m_dd,J,nn,iter1,burnin1,alpha,d,gamma,nu,bigK);
     % aggiorno gli iperparametri utilizzando il particle filter
             [ alpha, d ,gamma ,nu , M_parametri]=Filter_iperparametri(...
           mjk,m_j_dot,m_dd,m_dot_k,nj_dot_k,J,nn,bigK,N_iter,M_parametri);
       
       % Aggiorno i parametri di Oracle
       
       if sum(KuniOracle==newobsOrac)==0
           KuniOracle=[KuniOracle newobsOrac];
           newobsindOra(i)=1;
           % aggiorno la missin mass
           for j=1:J
               newlabel=find(pop{j}==newobsOrac);
               if isnan(newlabel)==0
                   missingmass(j)=missingmass(j)-freq{j}(newlabel);
               end
           end
       end
      
       
       % aggiorno i parametri di GT
       
       if sum(KuniGT==newobsGT)==0
           KuniGT=[KuniGT newobsGT];
           unGT{armchosenGT}=[unGT{armchosenGT} newobsGT];
           unGT_tot=[unGT_tot , newobsGT];
           KuniGT=[KuniGT newobsGT];
           newobsindGT(i)=1;
       end
       if sum(unGT_tot==newobsGT)
           for j=1:J
               if sum(unGT{j}==newobsGT)
                   indice=find(unGT{j}==newobsGT);
                   unGT{j}(indice)=[];
               end
           end
           indice=find(unGT_tot==newobsGT);
           unGT_tot(indice)=[];
       end      
       nGT(armchosenGT)=nGT(armchosenGT)+1;
   end
      
   % risultati MCMC
   
   distcumHPY=cumsum(newobsindHPY);
   distcumUni=cumsum(newobsindUni);
   distcumOra=cumsum(newobsindOra);
   distcumGT=cumsum(newobsindGT);
   
   DATAfinal.HPY(III,:)=distcumHPY;
   DATAfinal.uniform(III,:)=distcumUni;
   DATAfinal.Oracle(III,:)=distcumOra;
   DATAfinal.GoodTuring(III,:)=distcumGT;
   
   WEIGTHS.weigth_HPY(III,:)=w_HPY;
   WEIGTHS.weight_uniform(III,:)=w_Uni;
   WEIGTHS.weight_Oracle(III,:)=w_Ora;
   WEIGTHS.weight_GoodTuring(III,:)=w_GT;
    
III
    
    
end

   % plotto la scoperta di nuove
   plot(1:addsample,sum(DATAfinal.HPY)/Runs,'r');
   hold on 
   plot(1:addsample,sum(DATAfinal.uniform)/Runs,'g');
   plot(1:addsample,sum(DATAfinal.Oracle)/Runs,'b');
   plot(1:addsample,sum(DATAfinal.GoodTuring)/Runs,'k');
   
   legend('HPY','Uniform','Oracle','UCB','Location','NorthWest');
