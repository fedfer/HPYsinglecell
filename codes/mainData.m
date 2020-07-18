clear all;
close all;

addpath('/work/sta790/ff31/HPYsinglecell/codes')  
addpath('/work/sta790/ff31/HPYsinglecell/data')  

% total number of iterations 
Runs=100;

% Load data
fetal=csvread('fetal_dict2.txt');
adult=csvread('adult_dict2.txt');
embryo=csvread('embryo_dict2.txt');
newborn=csvread('newborn_dict2.txt');

J = 4;

% number of MCMC iterations
iter=35000;
burnin=15000;

% number of iterations for particle filter
% smaller than iter-burnin
N_iter=1000;

% normalizing parameter of GT strategy, in order to give some prob to be selected to 
% populations that have u_Gt = 0
alpha_GT = 0.1;

% initial sample
n_init=50*ones(J,1);
% additional samples
addsample=20;
% units sampled each additional trial
n_inc = 25;


% Storage 
M=zeros(Runs,addsample);
DATAfinal=struct('HPY',M,'uniform',M,'Oracle',M,'GoodTulming',M);

%weights
M=zeros(Runs,J);
WEIGTHS=struct('weight_HPY',M,'weight_uniform',M,'weight_Oracle',M,'weight_GoodTulming',M);
clear M;


pop=cell(1,J);
freq=cell(1,J);
   
pop{1} = fetal(:,1); freq{1} = fetal(:,2);
pop{2} = adult(:,1); freq{2} = adult(:,2);
pop{3} = embryo(:,1); freq{3} = embryo(:,2);
pop{4} = newborn(:,1); freq{4} = newborn(:,2);

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


Kini=unique(dati_totali); % distinct species accross population
tot_dist=length(Kini); % total distincts

%% Algorithm to choose from which sample the next oservation for the competitors
%  Unif, HPY Oracle GT
M0=1;
V0=4;
bigK=length(Kini);

[M_Tavoli, M_l_star, M_parametri, Dati_star, k_popolazioni]=posterior_K(data,M0,V0,J,n_init,iter,burnin);


[mjk_ini, m_dot_k_ini, m_j_dot_ini, m_dd_ini, alpha, d, gamma, nu]=estimate_parameters(M_l_star,tot_dist,M_parametri,J,iter-burnin);
gamma = repmat(gamma,Runs,1);
nu = repmat(nu,Runs,1);
alpha = repmat(alpha,Runs,1);
d = repmat(d,Runs,1);


M_parametri_ini=M_parametri(end-N_iter+1:end,:);


for III=1:Runs
    % Uniform sampling
    supp=1:J;
    masses=ones(1,J)/J;
    KuniUni=Kini;
    bigK=length(Kini);
    
    % HPY
    nn=n_init;
    KuniHPY=Kini;
    mjk=mjk_ini;
    m_dot_k=m_dot_k_ini;
    m_j_dot=m_j_dot_ini;
    m_dd=m_dd_ini;
    M_parametri=M_parametri_ini;
    
    %  n.k
    nj_dot_k=zeros(J,tot_dist);
    for j=1:J
        for w=1:tot_dist
            nj_dot_k(j,w) = sum(data{j}==KuniHPY(w));
        end
    end
    
    
    % Oracle
    KuniOracle=Kini;
    missingmass=zeros(1,J);
    for j=1:J
        labelsseen=ismember(pop{j},KuniOracle);

        missingmass(j) = 1-sum(freq{j}(labelsseen));
    end
    
    % Good-Tulming
    GT_data=data;
    KuniGT=Kini;
    
    newobsindHPY=zeros(1,addsample);
    newobsindOra=zeros(1,addsample);
    newobsindUni=zeros(1,addsample);
    newobsindGT=zeros(1,addsample);
    
    w_HPY=zeros(1,J);
    w_Ora=zeros(1,J);
    w_Uni=zeros(1,J);
    w_GT=zeros(1,J);
    
    
    
    
    
    
    for i=1:addsample
        
        % Uniform
        unif=gendiscr(supp,masses);
        w_Uni(unif)=w_Uni(unif)+1;
        
        
        % HPY
        betadraws=zeros(1,J);
        K_j_post = zeros(1,J);
        betazero=betarnd(gamma(III,1)+nu(III,1)*bigK,m_dd-bigK*nu(III,1));

        for j=1:J
            betadraws(j)=betarnd(betazero*(alpha(j)+(m_j_dot(j))), ((1-betazero)*(alpha(j)+m_j_dot(j))+ nn(j)- m_j_dot(j)));
            K_j_post(j) = E_Kjl_simplified(betadraws(j),nu(III,1),gamma(III,1),bigK,...
                alpha(j),d(j),betazero,m_j_dot(j),n_inc);
        end
  
        [v_max, armchosen]= max(K_j_post);
        clear v_max;
        w_HPY(armchosen)=w_HPY(armchosen)+1;
        
        % Oracle
        [v_max armchosenOrac]=max(missingmass);
        clear v_max;
        if length(armchosenOrac)>1
            Oracprob=ones(1,length(armchosenOrac))/armchosenOrac;
            armchosenOrac=gendiscr(armchosenOrac,Oracprob);
        end
        w_Ora(armchosenOrac)=w_Ora(armchosenOrac)+1;
        
        % GT
        u_GT=zeros(1,J);
        for j=1:J
            t = size(data{j},2)/n_inc;
            f = makeFinger(data{j});            
            u_GT(j) = U_GT(f,t) + alpha_GT;
        end
        u_GT = u_GT./sum(u_GT);
        armchosenGT = gendiscr(1:J,u_GT);
        w_GT(armchosenGT)=w_GT(armchosenGT)+1;
        
        % sample a unit for each population
        newobservations=zeros(n_inc,J);
        for j=1:J
            for jj=1:n_inc
                newobservations(jj,j) = gendiscr(pop{j},freq{j});
            end
        end
        newunif=newobservations(:,unif);
        newobsHPY=newobservations(:,armchosen);
        newobsOrac=newobservations(:,armchosenOrac);
        newobsGT=newobservations(:,armchosenGT);
        
        % update Uniform
        newobsindUni(i) = 0;
        for jj=1:n_inc
        if sum(KuniUni==newunif(jj))==0
            newobsindUni(i)=newobsindUni(i)+1;
            KuniUni=[KuniUni newunif(jj)];
        end
        end
        
        % save old parameters for new filtering algorithm
        mjk_old = mjk;m_j_dot_old = m_j_dot;m_dd_old = m_dd; 
        m_dot_k_old = m_dot_k;nj_dot_k_old = nj_dot_k; nn_old = nn;
        bigK_old = bigK;
        
        % Update HPY
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
  
            olddistinct=find(KuniHPY==newobsHPY(jj));
            probnewold=[nj_dot_k(armchosen,olddistinct)-d(armchosen)*mjk(armchosen,olddistinct),...
                (alpha(armchosen)+m_j_dot(armchosen)*d(armchosen))*((m_dot_k(olddistinct)-nu(III,1))/(gamma(III,1)+m_dd))];
            probnewold=probnewold/sum(probnewold);
            bern=binornd(1,probnewold(1));
            if bern==0
                m_j_dot(armchosen)=m_j_dot(armchosen)+1;
                m_dot_k(olddistinct)= m_dot_k(olddistinct)+1;
                mjk(armchosen,olddistinct)=mjk(armchosen,olddistinct)+1;
                m_dd=m_dd+1;
            end
            nj_dot_k(armchosen,olddistinct)=nj_dot_k(armchosen,olddistinct)+1;
        end
        
        nn(armchosen)=nn(armchosen)+1;
        end
        
        [alpha(III,:), d(III,:) ,gamma(III,1) ,nu(III,1) , M_parametri]=Filter_iperparametri_v1(...
            mjk,m_j_dot,m_dd,m_dot_k,nj_dot_k,J,nn,bigK,N_iter,M_parametri,...
            mjk_old, m_j_dot_old,m_dd_old,m_dot_k_old,nj_dot_k_old,nn_old,...
            bigK_old);
        
    
        
        % update Oracle
        newobsindOra(i) = 0;
        for jj=1:n_inc
        if sum(KuniOracle==newobsOrac(jj))==0
            KuniOracle=[KuniOracle newobsOrac(jj)];
            newobsindOra(i)=newobsindOra(i)+1;
            % missin mass
            for j=1:J
                newlabel=find(pop{j}==newobsOrac(jj));
                if isnan(newlabel)==0
                    missingmass(j)=missingmass(j)-freq{j}(newlabel);
                end
            end
        end
        end
        
        % update GT
        newobsindGT(i)=0;
        for jj=1:n_inc
        if sum(KuniGT==newobsGT(jj))==0
            KuniGT=[KuniGT newobsGT(jj)];
            newobsindGT(i)=newobsindGT(i)+1;
        end
        end
        GT_data{armchosenGT} = [GT_data{armchosenGT} newobsGT'];
        
    end
    
    % results MCMC
    distcumHPY=cumsum(newobsindHPY);
    distcumUni=cumsum(newobsindUni);
    distcumOra=cumsum(newobsindOra);
    distcumGT=cumsum(newobsindGT);
    
    DATAfinal.HPY(III,:)=distcumHPY;
    DATAfinal.uniform(III,:)=distcumUni;
    DATAfinal.Oracle(III,:)=distcumOra;
    DATAfinal.GoodTulming(III,:)=distcumGT;
    
    WEIGTHS.weigth_HPY(III,:)=w_HPY;
    WEIGTHS.weight_uniform(III,:)=w_Uni;
    WEIGTHS.weight_Oracle(III,:)=w_Ora;
    WEIGTHS.weight_GoodTulming(III,:)=w_GT;
    
    III
    
    
end


