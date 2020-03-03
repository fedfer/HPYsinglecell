% Same as main_highpop.m but with PARFOR Loop

clear all;
close all;

% add path for Computing
addpath('/work/sta790/ff31/HPYsinglecell/codes')  
addpath('/work/sta790/ff31/HPYsinglecell/codes/hyper')  
% addpath('/Users/felpo/MATLAB/projects/untitled/codes')  
% addpath('/Users/felpo/MATLAB/projects/untitled/codes/hyper')  

%add the Good-Tulming estimator of Bianca

% numero di iterazioni su cui fare average
% Runs=48;
Runs=5;

% numero delle popolazioni
% J = 100;
J = 10;

% Setting paper 
% numero totale delle specie tra tutte le popolazioni
N=2000;
% parametri per la Zipf
Zipfpar=[1.3; 1.3; 1.3; 1.3; repelem(2,J - 4).'];
% numero delle specie presenti nelle J popolazioni
% NN=250*ones(J,1);
NN=25*ones(J,1);


% Reviewer Answer 6 
% parametri per la Zipf
% Zipfpar=[repelem(2,33).'; repelem(2.1,33).'; repelem(1.9,34).'];

% numero di iterazioni in MCMC per il numero di tavoli e dei parametri di
% HPY dato il campione iniziale
% iter = 35000;
% burnin = 15000;
iter=350;
burnin=150;

% Numero di iterazioni per il particle filter: il numero delle iterazioni
% deve essere inferiore a iter-burnin
% N_iter=1000;
N_iter=100;

% normalizing parameter of GT strategy, in order to give some prob to be selected to 
% populations that have u_Gt = 0
alpha_GT = 0.1;

% ampiezza del campione iniziale
n_init=20*ones(J,1);
n_init_hyper=20*ones(J,1);
% lunghezza del campione addizionale
% addsample=50;
addsample=5;
% units sampled each additional trial
n_inc = 50;
n_inc_hyper = 50;


% Dati finali
M=zeros(Runs,addsample);
%DATAfinal=struct('HPY',M,'uniform',M,'Oracle',M,'GoodTulming',M);
DATAfinal_HPY = M;
DATAfinal_HPY_hyper = M;
DATAfinal_uniform = M;
DATAfinal_Oracle = M;
DATAfinal_GoodTulming = M;


%weights
M=zeros(Runs,J);
%WEIGTHS=struct('weight_HPY',M,'weight_uniform',M,'weight_Oracle',M,'weight_GoodTulming',M);
WEIGTHS_HPY = M;
WEIGTHS_HPY_hyper = M;
WEIGTHS_uniform = M;
WEIGTHS_Oracle = M;
WEIGTHS_GoodTulming = M;
clear M;

% simulare la legge vera in ogni popolazione
labels=1:N;

%frequenze nella popolazione
freq=cell(1,length(NN));
for i=1:J
    freq{i}=zeros(1,NN(i));
end

for i=1:J
    freq1=1:NN(i);
    freq2=sum(freq1.^(-Zipfpar(i)));
    for j=1:NN(i)
        % Zipf
        freq{i}(j)=((1/j)^Zipfpar(i))/freq2;
        
        % Uniform
        % freq{i}(j) = 1/NN(i);
    end
end

% labels delle specie nelle varie popolazioni
pop=cell(1,J);
for i=1:J
    labels_corr=labels;
    for h=1:NN(i)
        pop{i}(h)=gendiscr(labels_corr,ones(1,length(labels_corr))/length(labels_corr));
        ind=find(pop{i}(h)==labels_corr);
        labels_corr(ind)=[];
    end
end

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
Kini_hyper=unique(dati_totali);
%numero di specie iniziali distinte
tot_dist=length(Kini);
tot_dist_hyper=length(Kini);
%% Inferenza MCMC sul numero di specie e sul numero di parametri
% questi parametri sono quelli riferiti alle prior degli hyperparametri
% del PY
M0=1;
V0=4;
bigK=length(Kini);
bigK_hyper=length(Kini);

% aggiornamento con l'algoritmo marginale che ho sviluppato con Antonio
% ed Igor
[M_Tavoli_hyper,M_l_star_hyper,M_parametri_hyper,Dati_star_hyper,k_popolazioni_hyper]=posterior_K_hyper(data,M0,V0,J,n_init,iter,burnin);
[M_Tavoli,M_l_star,M_parametri,Dati_star,k_popolazioni]=posterior_K(data,M0,V0,J,n_init,iter,burnin);


% ora stimo i parametri che mi servono dopo
% hyper 
[mjk_ini_hyper, m_dot_k_ini_hyper, m_j_dot_ini_hyper, m_dd_ini_hyper,alpha_hyper, ...
    d_hyper,gamma_hyper,nu_hyper,theta_h,a0,b0]=stima_parametri_hyper(M_l_star_hyper,tot_dist_hyper, ...
    M_parametri_hyper,J,iter-burnin);
gamma_hyper = repmat(gamma_hyper,Runs,1);
nu_hyper = repmat(nu_hyper,Runs,1);
alpha_hyper = repmat(alpha_hyper,Runs,1);
d_hyper = repmat(d_hyper,Runs,1);
theta_h = repmat(theta_h,Runs,1);
a0 = repmat(a0,Runs,1);
b0 = repmat(b0,Runs,1);

% HPT-TS without hyperprior
[mjk_ini, m_dot_k_ini, m_j_dot_ini, m_dd_ini,alpha,d,gamma,nu]=stima_parametri(M_l_star,tot_dist,M_parametri,J,iter-burnin);
gamma = repmat(gamma,Runs,1);
nu = repmat(nu,Runs,1);
alpha = repmat(alpha,Runs,1);
d = repmat(d,Runs,1);
theta_h = repmat(theta_h,Runs,1);
% creo un vettore di pesi per la mistura delle normali dei parametri: ne
% prendo solo una parte.
M_parametri_ini=M_parametri(end-N_iter+1:end,:);
M_parametri_ini_hyper=M_parametri_hyper(end-N_iter+1:end,:);


% storage 
% M_parametri_storage = repmat(M_parametri_ini,Runs,addsample,1);


%% Algoritmi vari per scegliere da dove campionare la prossima
% osservazione: Unif, HPY Oracle GT
disp('start parfor loop')

parfor III=1:Runs
%for III=1:Runs
    % Uniform sampling
    supp=1:J;
    masses=ones(1,J)/J;
    KuniUni=Kini;
    bigK=length(Kini);
    bigK_hyper=length(Kini_hyper);
    
    % HPY
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
    
    
    % HPY _hyper
    nn_hyper=n_init_hyper;
    KuniHPY_hyper=Kini_hyper;
    mjk_hyper=mjk_ini_hyper;
    m_dot_k_hyper=m_dot_k_ini_hyper;
    m_j_dot_hyper=m_j_dot_ini_hyper;
    m_dd_hyper=m_dd_ini_hyper;
    M_parametri_hyper=M_parametri_ini_hyper;
    
    % calcolo n.k
    nj_dot_k_hyper=zeros(J,tot_dist_hyper);
    for j=1:J
        for w=1:tot_dist_hyper
            nj_dot_k_hyper(j,w) = sum(data{j}==KuniHPY_hyper(w));
        end
    end
    
    
    % Oracle
    KuniOracle=Kini;
    missingmass=zeros(1,J);
    for j=1:J
        labelsseen=ismember(pop{j},KuniOracle);
        % ritorna un vettore di 1 e 0: ho 1 se pop{j} ? elemento di
        % KuniOracle
        missingmass(j) = 1-sum(freq{j}(labelsseen));
    end
    
    % Good-Tulming
    GT_data=data;
    KuniGT=Kini;
    
    % indicatori delle scoperte di nuove specie
    newobsindHPY=zeros(1,addsample);
    newobsindHPY_hyper=zeros(1,addsample);
    newobsindOra=zeros(1,addsample);
    newobsindUni=zeros(1,addsample);
    newobsindGT=zeros(1,addsample);
    
    % indicatori per i pesi
    w_HPY=zeros(1,J);
    w_HPY_hyper=zeros(1,J);
    w_Ora=zeros(1,J);
    w_Uni=zeros(1,J);
    w_GT=zeros(1,J);
    
    % spezio per contenere i risultati relativi alla stima degli
    % iperparametri
    
    
    %% Inizio degli algoritmi
    
    
    
    for i=1:addsample
        
        alpha_temp = alpha(III,:);
        d_temp = d(III,:);
        alpha_temp_hyper = alpha_hyper(III,:);
        d_temp_hyper = d_hyper(III,:);
        
        % scegliamo la popolazione per la stratqegia Uni
        unif=gendiscr(supp,masses);
        w_Uni(unif)=w_Uni(unif)+1;
        
        
        % scegliamo la popolazone per la strategia HPY
        betadraws=zeros(1,J);
        K_j_post = zeros(1,J);
        betazero=betarnd(gamma(III,1)+nu(III,1)*bigK,m_dd-bigK*nu(III,1));
        for j=1:J
            betadraws(j)=betarnd(betazero*(alpha_temp(j)+(m_j_dot(j))), ((1-betazero)*(alpha_temp(j)+m_j_dot(j))+ nn(j)- m_j_dot(j)));
            K_j_post(j) = E_Kjl_simplified(betadraws(j),nu(III,1),gamma(III,1),bigK,...
                alpha_temp(j),d_temp(j),betazero,m_j_dot(j),n_inc);
        end
  
        [v_max, armchosen]= max(K_j_post);
        w_HPY(armchosen)=w_HPY(armchosen)+1;
        
        
        % scegliamo la popolazone per la strategia HPY hyper
        betadraws_hyper=zeros(1,J);
        K_j_post_hyper = zeros(1,J);
        betazero_hyper=betarnd(gamma_hyper(III,1)+nu_hyper(III,1)*bigK_hyper,m_dd_hyper-bigK_hyper*nu_hyper(III,1));
        for j=1:J
            betadraws_hyper(j)=betarnd(betazero_hyper*(alpha_temp_hyper(j)+(m_j_dot_hyper(j))),...
                ((1-betazero_hyper)*(alpha_temp_hyper(j)+m_j_dot_hyper(j))+ nn_hyper(j)- m_j_dot_hyper(j)));
            K_j_post_hyper(j) = E_Kjl_simplified(betadraws_hyper(j),nu_hyper(III,1),gamma_hyper(III,1),bigK_hyper,...
                alpha_temp_hyper(j),d_temp_hyper(j),betazero_hyper,m_j_dot_hyper(j),n_inc_hyper);
        end
  
        [v_max_hyper, armchosen_hyper]= max(K_j_post_hyper);
        w_HPY_hyper(armchosen_hyper)=w_HPY_hyper(armchosen_hyper)+1;
        
        
        
        % scegliamo la popolazione per la strategia Oracle
        [v_max armchosenOrac]=max(missingmass);
        %clear v_max;
        if length(armchosenOrac)>1
            Oracprob=ones(1,length(armchosenOrac))/armchosenOrac;
            armchosenOrac=gendiscr(armchosenOrac,Oracprob);
        end
        w_Ora(armchosenOrac)=w_Ora(armchosenOrac)+1;
        
        % scegliamo la popolazione per la strategia GT

        u_GT=zeros(1,J);
        for j=1:J
            t = size(data{j},2)/n_inc;
            f = makeFinger(data{j});            
            u_GT(j) = U_GT(f,t) + alpha_GT;
        end
        u_GT = u_GT./sum(u_GT);
        armchosenGT = gendiscr(1:J,u_GT);
        w_GT(armchosenGT)=w_GT(armchosenGT)+1;
        
        % campioniamo un nuovo valore per ogni popolazione scelta
        newobservations=zeros(n_inc,J);
        for j=1:J
            for jj=1:n_inc
                newobservations(jj,j) = gendiscr(pop{j},freq{j});
            end
        end
        newunif=newobservations(:,unif);
        newobsHPY=newobservations(:,armchosen);
        newobsHPY_hyper=newobservations(:,armchosen_hyper);
        newobsOrac=newobservations(:,armchosenOrac);
        newobsGT=newobservations(:,armchosenGT);
        
        % aggiornamento dell'uniforme
        newobsindUni(i) = 0;
        for jj=1:n_inc
        if sum(KuniUni==newunif(jj))==0
            newobsindUni(i)=newobsindUni(i)+1;
            KuniUni=[KuniUni newunif(jj)];
        end
        end
        
        % UPDATE HPY
        newobsindHPY(i) = 0;
        for jj=1:n_inc
        if sum(KuniHPY==newobsHPY(jj))==0
            bigK=bigK+1;
            m_j_dot(armchosen)=m_j_dot(armchosen)+1;
            m_dd=m_dd+1;
            KuniHPY=[KuniHPY newobsHPY(jj)];
            newobsindHPY(i)=newobsindHPY(i)+1;
            m_dot_k=[m_dot_k 1];
            nj_dot_k=[nj_dot_k , zeros(J,1)];
            mjk=[mjk , zeros(J,1)];
            nj_dot_k(armchosen,bigK)=1;
            mjk(armchosen,bigK)=1;
        else
            
            % calcolo la probabilit? di sedersi a un tavolo vecchio se
            % l'osservazione ? vecchia
            olddistinct=find(KuniHPY==newobsHPY(jj));
            probnewold=[nj_dot_k(armchosen,olddistinct)-d_temp(armchosen)*mjk(armchosen,olddistinct),...
                (alpha_temp(armchosen)+m_j_dot(armchosen)*d_temp(armchosen))*((m_dot_k(olddistinct)-nu(III,1))/(gamma(III,1)+m_dd))];
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
        end
        [alpha(III,:), d(III,:) ,gamma(III,1) ,nu(III,1), M_parametri]=Filter_iperparametri(...
            mjk,m_j_dot,m_dd,m_dot_k,nj_dot_k,J,nn,bigK,N_iter,M_parametri);
        
        
        % UPDATE HPY HYPER _hyper
        newobsindHPY_hyper(i) = 0;
        for jj=1:n_inc_hyper
        if sum(KuniHPY_hyper==newobsHPY_hyper(jj))==0
            bigK_hyper=bigK_hyper+1;
            m_j_dot_hyper(armchosen_hyper)=m_j_dot_hyper(armchosen_hyper)+1;
            m_dd_hyper=m_dd_hyper+1;
            KuniHPY_hyper=[KuniHPY_hyper newobsHPY_hyper(jj)];
            newobsindHPY_hyper(i)=newobsindHPY_hyper(i)+1;
            m_dot_k_hyper=[m_dot_k_hyper 1];
            nj_dot_k_hyper=[nj_dot_k_hyper , zeros(J,1)];
            mjk_hyper=[mjk_hyper , zeros(J,1)];
            nj_dot_k_hyper(armchosen_hyper,bigK_hyper)=1;
            mjk(armchosen_hyper,bigK_hyper)=1;
        else

            olddistinct_hyper=find(KuniHPY_hyper==newobsHPY_hyper(jj));
            probnewold_hyper=[nj_dot_k_hyper(armchosen_hyper,olddistinct_hyper)- ...
                d_temp_hyper(armchosen_hyper)*mjk_hyper(armchosen_hyper,olddistinct_hyper),...
                (alpha_temp_hyper(armchosen_hyper)+ ...
                m_j_dot_hyper(armchosen_hyper)*d_temp_hyper(armchosen_hyper))*((m_dot_k_hyper(olddistinct_hyper)-...
                nu_hyper(III,1))/(gamma_hyper(III,1)+m_dd_hyper))];
            probnewold_hyper=probnewold_hyper/sum(probnewold_hyper);
            bern_hyper=binornd(1,probnewold_hyper(1));
            if bern_hyper==0
                % l'osservazione forma un nuovo tavolo
                m_j_dot_hyper(armchosen_hyper)=m_j_dot_hyper(armchosen_hyper)+1;
                m_dot_k_hyper(olddistinct_hyper)= m_dot_k_hyper(olddistinct_hyper)+1;
                mjk_hyper(armchosen_hyper,olddistinct_hyper)=mjk_hyper(armchosen_hyper,olddistinct_hyper)+1;
                m_dd_hyper=m_dd_hyper+1;
            end
            nj_dot_k_hyper(armchosen_hyper,olddistinct_hyper)=nj_dot_k_hyper(armchosen_hyper,olddistinct_hyper)+1;
        end
        
        nn_hyper(armchosen_hyper)=nn_hyper(armchosen_hyper)+1;
        end
        [alpha_hyper(III,:), d_hyper(III,:) ,gamma_hyper(III,1) ,nu_hyper(III,1), ...
            theta_h(III,1),a0(III,1),b0(III,1), M_parametri_hyper]=Filter_iperparametri_hyper(...
            mjk_hyper,m_j_dot_hyper,m_dd_hyper,m_dot_k_hyper,nj_dot_k_hyper,J,nn_hyper,...
            bigK_hyper,N_iter,M_parametri_hyper);
        
        
        % Update Oracle
        newobsindOra(i) = 0;
        for jj=1:n_inc
        if sum(KuniOracle==newobsOrac(jj))==0
            KuniOracle=[KuniOracle newobsOrac(jj)];
            newobsindOra(i)=newobsindOra(i)+1;
            % aggiorno la missin mass
            for j=1:J
                newlabel=find(pop{j}==newobsOrac(jj));
                if isnan(newlabel)==0
                    missingmass(j)=missingmass(j)-freq{j}(newlabel);
                end
            end
        end
        end
        
        % aggiorno i parametri di GT
        newobsindGT(i)=0;
        for jj=1:n_inc
        if sum(KuniGT==newobsGT(jj))==0
            KuniGT=[KuniGT newobsGT(jj)];
            newobsindGT(i)=newobsindGT(i)+1;
        end
        end
        GT_data{armchosenGT} = [GT_data{armchosenGT} newobsGT'];
        
    end
    
    % risultati MCMC
    
    distcumHPY=cumsum(newobsindHPY);
    distcumHPY_hyper=cumsum(newobsindHPY_hyper);
    distcumUni=cumsum(newobsindUni);
    distcumOra=cumsum(newobsindOra);
    distcumGT=cumsum(newobsindGT);
    
    DATAfinal_HPY(III,:)=distcumHPY;
    DATAfinal_HPY_hyper(III,:)=distcumHPY_hyper;
    DATAfinal_uniform(III,:)=distcumUni;
    DATAfinal_Oracle(III,:)=distcumOra;
    DATAfinal_GoodTulming(III,:)=distcumGT;
    
    WEIGTHS_HPY(III,:)=w_HPY;
    WEIGTHS_HPY_hyper(III,:)=w_HPY_hyper;
    WEIGTHS_uniform(III,:)=w_Uni;
    WEIGTHS_Oracle(III,:)=w_Ora;
    WEIGTHS_GoodTulming(III,:)=w_GT;
    
    III
    
    % Storage just for one iteration of algorithm 
    % M_parametri_storage(III,:,i) = M_parametri;
    
    
end

%to save workspace
%M_par_save = M_parametri_storage(1,:,:);
%clear M_parametri_storage; 
save('J100_hyper.mat');

%save weights and add sample units
%xlswrite('HPY.xls',transpose(sum(DATAfinal.HPY)/Runs));
%xlswrite('unif.xls',transpose(sum(DATAfinal.uniform)/Runs));
%xlswrite('GT.xls',transpose(sum(DATAfinal.GoodTulming)/Runs));
%xlswrite('Oracle.xls',transpose(sum(DATAfinal.Oracle)/Runs));

%per caricare dei vecchi valori
% unif= csvread('unif.csv');