clear all;
close all;

% number of iterations
Runs=50;


%% Parameters
% total number of populations
J = 5;

% total number of species
N=3000;

% number of species in each pop
NN=250*ones(J,1);

% parameters of la Zipf
Zipfpar=[1.3; 1.3; 1.3; 1.3; repelem(2,J - 4).'];

% initial sample
n_init=30*ones(J,1);
% additional samples
addsample=300;

% number of MCMC iterations
iter=35000;
burnin=15000;

% number of iterations for particle filter
% smaller than iter-burnin
N_iter=1000;

% parameters Good-Turing
C=(1+sqrt(2))*sqrt(3);

% Final storage
M=zeros(Runs,addsample);
DATAfinal=struct('HPY',M,'uniform',M,'Oracle',M,'GoodTuring',M);

% Weights
M=zeros(Runs,J);
WEIGTHS=struct('weight_HPY',M,'weight_uniform',M,'weight_Oracle',M,'weight_GoodTuring',M);
clear M;


% sample the species in each population 
labels=1:N;

freq=cell(1,length(NN));
for i=1:J
    freq{i}=zeros(1,NN(i));
end

for i=1:J
    freq1=1:NN(i);
    freq2=sum(freq1.^(-Zipfpar(i)));
    for j=1:NN(i)
        freq{i}(j)=((1/j)^Zipfpar(i))/freq2;
    end
end

% labels in each population 
pop=cell(1,J);
for i=1:J
    labels_corr=labels;
    for h=1:NN(i)
        pop{i}(h)=gendiscr(labels_corr,ones(1,length(labels_corr))/length(labels_corr));
        ind=find(pop{i}(h)==labels_corr);
        labels_corr(ind)=[];
    end
end

% sample data for each population 
% data{j} = vector that contains the values for population j
data=cell(1,J);
for i=1:J
    data{i}=zeros(1,n_init(i));
end
for j=1:J
    for i=1:n_init(j)
        data{j}(i)=gendiscr(pop{j},freq{j});
    end
end
data_total=cell2mat(data);
Kini=unique(data_total); % distinct species accross population
tot_dist=length(Kini); % total distincts


% Inference on HPY hyperparameters
M0=1;
V0=4;
bigK=length(Kini);

% update with marginal algorithm 
[M_Tables M_l_star M_parameters Data_star k_pop]=posterior_K(data,M0,V0,J,n_init,iter,burnin);

% estimates the parameters that we need for bandits
[mjk_ini, m_dot_k_ini, m_j_dot_ini, m_dd_ini alpha d gamma nu]=estimate_parameters(M_l_star,tot_dist,M_parameters,J,iter-burnin);

M_parameter_ini=M_parameters(end-N_iter+1:end,:);

%% Algorithm to choose from which sample the next oservation for the competitors
%  Unif, HPY Oracle GT

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
    M_parameters=M_parameter_ini;
    
    % compute n.k
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
    
    % Good--Turing
    KuniGT=Kini;
    nGT=n_init;
  
    [data_uni freq_uni]=cluster_fct(cell2mat(data));
    unGT_tot=data_uni(freq_uni==1);
    clear data_uni freq_uni;

    unGT=cell(1,J);
    for j=1:J
        indici=ismember(data{j},unGT_tot);
        unGT{j}=unique(data{j}(indici));
    end
    clear data_uni freq_uni;
    
    
    % vector for new species
    newobsindHPY=zeros(1,addsample);
    newobsindOra=zeros(1,addsample);
    newobsindUni=zeros(1,addsample);
    newobsindGT=zeros(1,addsample);
    
    % vector for weights
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
        betazero=betarnd(gamma+nu*bigK,m_dd-bigK*nu);
        for j=1:J
            betadraws(j)=betarnd(betazero*(alpha(j)+(m_j_dot(j)*d(j))), ((1-betazero)*(alpha(j)+m_j_dot(j)*d(j))+ nn(j)- m_j_dot(j)*d(j)));
        end
        [v_max, armchosen]= max(betadraws);
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
        
        
        % Sample a unit for each population
        newobservations=zeros(1,J);
        for j=1:J
            newobservations(j) = gendiscr(pop{j},freq{j});
        end
        newunif=newobservations(unif);
        newobsHPY=newobservations(armchosen);
        newobsOrac=newobservations(armchosenOrac);
        newobsGT=newobservations(armchosenGT);
        
        % update Uniform 
        if sum(KuniUni==newunif)==0
            newobsindUni(i)=1;
            KuniUni=[KuniUni newunif];
        end
        
        % update HPY hyperparameters
        % TP=2 New table and new dish 
        % TP=1 New table and old dish
        % TP=0 Old table and old dish
        TP=0;
        olddistinct=find(KuniHPY==newobsHPY);
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
            TP=2;
        else
            

            olddistinct=find(KuniHPY==newobsHPY);
            probnewold=[nj_dot_k(armchosen,olddistinct)-d(armchosen)*mjk(armchosen,olddistinct),...
                (alpha(armchosen)+m_j_dot(armchosen)*d(armchosen))*((m_dot_k(olddistinct)-nu)/(gamma+m_dd))];
            probnewold=probnewold/sum(probnewold);
            bern=binornd(1,probnewold(1));
            if bern==0

                m_j_dot(armchosen)=m_j_dot(armchosen)+1;
                m_dot_k(olddistinct)= m_dot_k(olddistinct)+1;
                mjk(armchosen,olddistinct)=mjk(armchosen,olddistinct)+1;
                m_dd=m_dd+1;
                TP = 1;
            end
            nj_dot_k(armchosen,olddistinct)=nj_dot_k(armchosen,olddistinct)+1;
        end
        nn(armchosen)=nn(armchosen)+1;
        

        % Particle filter for hyperparameters
        [ alpha, d ,gamma ,nu , M_parameters]=Filter_hyperparameters(...
            mjk,m_j_dot,m_dd,m_dot_k,nj_dot_k,J,nn,bigK,N_iter,M_parameters,...
            armchosen,olddistinct,TP);
        
        
        % update Oracle
        if sum(KuniOracle==newobsOrac)==0
            KuniOracle=[KuniOracle newobsOrac];
            newobsindOra(i)=1;
            %  missin mass
            for j=1:J
                newlabel=find(pop{j}==newobsOrac);
                if isnan(newlabel)==0
                    missingmass(j)=missingmass(j)-freq{j}(newlabel);
                end
            end
        end
        
        
        % update GT
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
                    index=find(unGT{j}==newobsGT);
                    unGT{j}(index)=[];
                end
            end
            index=find(unGT_tot==newobsGT);
            unGT_tot(index)=[];
        end
        nGT(armchosenGT)=nGT(armchosenGT)+1;
        
    end
    
    % results MCMC
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

