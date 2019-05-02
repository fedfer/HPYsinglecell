% funzione che calcola la discovery probability e la media delle nuove
% specie osservate nelle prossime m osservazioni per
% porzioni  del campione di 20, 40, 60, 80 , 100%
% considero il campione di ampiezza successiva m solo per il primo dataset,
% perché voglio calcolare la discovery probability solo per il primo
% dataset

function [specie_future media media_quantili P_posterior Valori_predittiva media_pr_nuova media_pr_nuova_quantili...
    Probabilita_nuova M   media_nuove_X1_non_X2  media_nuove_X1_non_X2_quantili ...
    media_distinte_nuove_X1_non_X2  media_distinte_nuove_X1_non_X2_quantili   ...
        media_vecchie_X1_non_X2  media_vecchie_X1_non_X2_quantili media_vecchie_condivise  ...
        media_vecchie_condivise_quantili  media_nuove  media_nuove_quantili]=prediction(M_X1,M_X2,M_T1,M_T2,M_parametri,N,n1,m1,n2)

M=round(m1*[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
X1=M_X1(1,1:n1);
X1_star=clusterizza(X1);
k1=length(X1_star);
X2=M_X2(1,1:n2);
X2_star=clusterizza(X2);
k2=length(X2_star);
k=length(unique([X1_star X2_star]));

% voglio calcolare le specie future che sono state campionate nella
% previsione ad m passi presenti nel secondo dataset ma non nel primo
specie_distinte_X2=ismember(X2_star,X1_star); % ritorna un 1 quando X2_star è comune anche ad X1_star
specie_distinte_X2=X2_star(specie_distinte_X2==0); % specie presenti in X2 ma non in X1 all'inizio
logical_species=ismember(X1_star,X2_star);
specie_distinte_X1=X1_star(logical_species==0);
specie_comuni=X1_star(logical_species==1);

% calcolo la media della probabilità di osservare un nuovo gene, su
% entrambi i campioni, cioè non deve mai essere stato osservato, in nessuna
% libreria. Anche la discovery probability su tutto
Probabilita_nuova=zeros(N,length(M));
specie_future=zeros(N,length(M));
osservazioni_nuove_X1_non_X2=zeros(N,length(M));
osservazioni_vecchie_X1_non_X2=zeros(N,length(M));
osservazioni_vecchie_condivise=zeros(N,length(M));
osservazioni_nuove=zeros(N,length(M));
specie_nuove_X1_non_X2=zeros(N,length(M));
clear nuove_specie;
media=zeros(1,length(M));
media_quantili=zeros(length(M),2);
P_posterior=cell(1,length(M));
Valori_predittiva=cell(1,length(M));
media_pr_nuova=zeros(1,length(M));
media_pr_nuova_quantili=zeros(length(M),2);
media_nuove_X1_non_X2=zeros(1,length(M));
media_nuove_X1_non_X2_quantili=zeros(length(M),2);
media_distinte_nuove_X1_non_X2=zeros(1,length(M));
media_distinte_nuove_X1_non_X2_quantili=zeros(length(M),2);
media_vecchie_X1_non_X2=zeros(1,length(M));
media_vecchie_X1_non_X2_quantili=zeros(length(M),2);
media_vecchie_condivise=zeros(1,length(M));
media_vecchie_condivise_quantili=zeros(length(M),2);
media_nuove=zeros(1,length(M));
media_nuove_quantili=zeros(length(M),2);
for i=1:length(M)
    % a posteriori per n=M(i)
    m=M(i);
    M_Xnm1=M_X1(:,1:n1+m);
    M_Tnm1=M_T1(:,1:n1+m);
    M_Xn2=M_X2(:,1:n2);
    M_Tn2=M_T2(:,1:n2);
    for j=1:N
        theta_0=M_parametri(j,1);
        sigma_0=M_parametri(j,2);
        theta=M_parametri(j,3);
        sigma=M_parametri(j,4);
        Xnm1=M_Xnm1(j,:);
        Tnm1=M_Tnm1(j,:);
        Xn2=M_Xn2(j,:);
        Tn2=M_Tn2(j,:);
        [X_nm1_star T_nm1_star XT_nm1_star n1_star q1_star l1_star] = distinti(Xnm1, Tnm1);
        L1=sum(l1_star);
        kj=length(unique([Xnm1 Xn2]));
        kj1=length(X_nm1_star);
        L2=length(unique(Tn2));
        Probabilita_nuova(j,i)=(theta_0+kj*sigma_0)*(theta+L1*sigma)/...
            ((theta_0+L1+L2)*(theta+n1+m));
        specie_future(j,i)=kj-k;
        % queste sono le osservazioni campionate dalle distinte proprie di X1
        % : posso anche contarle più volte
        osservazioni_nuove_X1_non_X2(j,i)=sum(ismember(Xnm1,specie_distinte_X2));
        % qui come prima ma prendo le distinte
        Xnm1_star=clusterizza(Xnm1);
        specie_nuove_X1_non_X2(j,i)=sum(ismember(Xnm1_star,specie_distinte_X2));
        osservazioni_vecchie_X1_non_X2(j,i)=sum(ismember(Xnm1,specie_distinte_X1));
        osservazioni_vecchie_condivise(j,i)=sum(ismember(Xnm1,specie_comuni));
        osservazioni_nuove(j,i)=sum(1-ismember(Xnm1,[X1_star X2_star]));
        j
    end
    % la funzione clusterizza prende in ingresso dei vettori riga
    [X_local n_local]=clusterizza(specie_future(:,i)');
    Valori_predittiva{i}=X_local';
    P_posterior{i}=n_local';
    clear X_local n_local;
    % normalizzo le frequneze
    sum(P_posterior{i})
    P_posterior{i}=P_posterior{i}/sum(P_posterior{i});
    media(i)=sum(Valori_predittiva{i}.*P_posterior{i});
    media_quantili(i,:)= quantile(Valori_predittiva{i},[0.025 0.975]);
    media_pr_nuova(i)=sum(Probabilita_nuova(:,i))/N;
    media_quantili(i,:)= quantile(Valori_predittiva{i},[0.025 0.975]);
    media_pr_nuova_quantili(i,:)=quantile(Probabilita_nuova(:,i),[0.025 0.975]);
    % calcolo la media delle m specie campionate nella previsione ad m
    % passi di X1 già presenti nel campione di base di X2 ma non in quello
    % di X1
    media_nuove_X1_non_X2(i)=sum(osservazioni_nuove_X1_non_X2(:,i))/N;
    media_nuove_X1_non_X2_quantili(i,:)=quantile(osservazioni_nuove_X1_non_X2(:,i),[0.025 0.975]);
    media_distinte_nuove_X1_non_X2(i)=sum(specie_nuove_X1_non_X2(:,i))/N;
    media_distinte_nuove_X1_non_X2_quantili(i,:)=quantile(specie_nuove_X1_non_X2(:,i),[0.025 0.975]);
    
    media_vecchie_X1_non_X2(i)=sum(osservazioni_vecchie_X1_non_X2(:,i))/N;
    media_vecchie_X1_non_X2_quantili(i,:)=quantile(osservazioni_vecchie_X1_non_X2(:,i),[0.025 0.975]);
    media_vecchie_condivise(i)=sum(osservazioni_vecchie_condivise(:,i))/N;
    media_vecchie_condivise_quantili(i,:)=quantile(osservazioni_vecchie_condivise(:,i),[0.025 0.975]);
    media_nuove(i)=sum(osservazioni_nuove(:,i))/N;
    media_nuove_quantili(i,:)=quantile(osservazioni_nuove(:,i),[0.025 0.975]);
end