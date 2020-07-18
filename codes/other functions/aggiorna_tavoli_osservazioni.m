% funzione che aggiorna  sia il tavolo che l'osservazione associata

function  [X1_new T1_new]=aggiorna_tavoli_osservazioni(X1_old,T1_old,X2,T2,theta_0,sigma_0,theta,sigma,n1,m1,M0,V0)


for i=(n1+1):(n1+m1)
    X1_i=X1_old;
    T1_i=T1_old;
    X1_i(i)=[];
    T1_i(i)=[];
    X_i=[X1_i X2];
    T_i=[T1_i T2];
    % cluster the franchise
    [X_star T_star XT_star n_star q_star l_star] = distinti(X_i, T_i);

    % function ismember: if the table of the entire franchise are in the first restaurant
    % compute the values for first restaurant

    logical_T=ismember(T_star,T1_i);
    indici1=(logical_T==1);
    XT1_star=XT_star(:,indici1);
    q1_star=q_star(indici1);
    kj=length(X_star);
    L=length(q_star);
    L1=length(q1_star);
    Pr=zeros(1,L1+kj+1);
    % probability of new table and new dish
    Pr(1)=(theta_0+kj*sigma_0)*(theta+L1*sigma)/((theta_0+L)*(theta+n1+m1-1));
    % probability of old table and old dish
    Pr(2:kj+1)=(theta+L1*sigma)*(l_star-sigma_0)/((theta_0+L)*(theta+n1+m1-1));
    % probability of old table and old dish
    Pr(kj+2:L1+kj+1)=(q1_star-sigma)/(theta+n1+m1-1);
    
    Pr_cum=cumsum(Pr);
    U=unifrnd(0,1);
    indice=find(U<Pr_cum);
    indice=indice(1);
    
    if indice==1
        X1_old(i)=normrnd(M0,sqrt(V0));
        T1_old(i)=normrnd(M0,sqrt(V0));
    elseif indice>1 && indice<=kj+1
        X1_old(i)=X_star(indice-1);
        T1_old(i)=normrnd(M0,sqrt(V0));
    else
        X1_old(i)=XT1_star(1,indice-kj-1);
        T1_old(i)=XT1_star(2,indice-kj-1);
    end
end

X1_new=X1_old;
T1_new=T1_old;
    