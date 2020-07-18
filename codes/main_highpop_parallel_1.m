% Same as main_highpop.m but designed for using PARFOR Loop

clear all;
close all;

% add path for Computing
addpath('/work/sta790/ff31/HPYsinglecell/codes/hyper')  

% number of iterations
Runs=100;

% total number of populations
J=100;

% Setting paper 
% number of species in each pop
N=20000;

% parameters of la Zipf
Zipfpar=[1.3; 1.3; 1.3; 1.3; repelem(2,J - 4).'];

% number of species in each pop
NN=2500*ones(J,1);

% Other setting
% Zipfpar=[repelem(2,33).'; repelem(2.1,33).'; repelem(1.9,34).'];

% number of MCMC iterations
iter=35000;
burnin=15000;

% number of iterations for particle filter
% smaller than iter-burnin
N_iter=1000;

% normalizing parameter of GT strategy, in order to give some prob to be selected to 
% populations that have u_Gt = 0
alpha_GT = 0.1;

n_init=20*ones(J,1); % initial sample
addsample=100; % extra samples at each step
n_inc = 50; % units sampled each additional trial

% Final storage
M=zeros(Runs,addsample);
DATAfinal_HPY = M;
DATAfinal_uniform = M;
DATAfinal_Oracle = M;
DATAfinal_GoodTulming = M;


% Weights
M=zeros(Runs,J);
WEIGTHS_HPY = M;
WEIGTHS_uniform = M;
WEIGTHS_Oracle = M;
WEIGTHS_GoodTulming = M;
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
        % Zipf
        freq{i}(j)=((1/j)^Zipfpar(i))/freq2;
        
        % Uniform
        % freq{i}(j) = 1/NN(i);
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

dati_totali=cell2mat(data);
Kini=unique(dati_totali); % distinct species accross population
tot_dist=length(Kini); % total distincts

% Inference on HPY hyperparameters
M0=1;
V0=4;
bigK=length(Kini);

% update with marginal algorithm 
[M_Tavoli M_l_star M_parametri Dati_star k_popolazioni]=posterior_K(data,M0,V0,J,n_init,iter,burnin);

% estimates the parameters that we need for bandits
[mjk_ini, m_dot_k_ini, m_j_dot_ini, m_dd_ini alpha d gamma nu]=stima_parametri(M_l_star,tot_dist,M_parametri,J,iter-burnin);
gamma = repmat(gamma,Runs,1);
nu = repmat(nu,Runs,1);
alpha = repmat(alpha,Runs,1);
d = repmat(d,Runs,1);

M_parametri_ini=M_parametri(end-N_iter+1:end,:);


disp('finish first estimation')

% save workspace
save('MCMC_estimation_paper.mat');

