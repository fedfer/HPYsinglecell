% function to update the tables only and not the observations
% j: index of population in which update the tables
function [Tj_new]=aggiorna_tavoli_n(j,J,Dati,Tavoli,theta_0,sigma_0,theta,sigma,n_init,M0,V0)

Xj=Dati{j};
Tj_old=Tavoli{j};

for i=1:n_init(j)
    x=Xj(i);
    Tj_i=Tj_old;
    Tj_i(i)=[];
    Xj_i=Xj;
    Xj_i(i)=[];
    % Find the tables that eat the same dish seperately for all J restaurants
    T_piatto_x=cell(1,J);
    Tavoli_totali_j=[];
    Tavoli_totali_j_piatto_x=[];
    for h=1:J
        if h==j
            T_piatto_x{h}=Tj_i(Xj_i==x);
        else
            T_piatto_x{h}=Tavoli{h}(Dati{h}==x);
            Tavoli_totali_j=[Tavoli_totali_j Tavoli{h}];
            Tavoli_totali_j_piatto_x=[Tavoli_totali_j_piatto_x T_piatto_x{h}];
        end
    end
    njx=length(T_piatto_x{j});

    % the table can be the same as one of which the dish x is eaten, or new
    % number of tables in first restaurant
    Lj=length(unique(Tj_i));
    % total number of tables without i-th restaurant
    L=Lj+length(unique(Tavoli_totali_j));
    if njx>0
        [Tj_star qj_star]=clusterizza(T_piatto_x{j});
        Lj_x=length(Tj_star);
        L_j_x=length(unique(Tavoli_totali_j_piatto_x));
        Pr=zeros(1,Lj_x+1);
        Pr(1)=(Lj_x+L_j_x-sigma_0)*(theta+Lj*sigma)/(theta_0+L);
        Pr(2:Lj_x+1)=(qj_star-sigma);
        Pr=Pr/sum(Pr);
    else
        Pr=1;
    end
    
    Pr_cum=cumsum(Pr);
    U=unifrnd(0,1);
    indice=find(U<Pr_cum);
    indice=indice(1);
    
    if indice==1
        Tj_old(i)=normrnd(M0,sqrt(V0));
    else
        Tj_old(i)=Tj_star(indice-1);
    end
end

Tj_new=Tj_old;
    
    