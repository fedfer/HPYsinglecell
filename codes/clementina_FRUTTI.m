% Ho due dataset
clear all;
close all;

a0=300;
b0=5;
a=300;
b=5;
var=10;
M0=1;
V0=4;

data_set=cell(1,2);
data_set{1}='FRUTTO1_FlavFr1.txt';
data_set{2}='FRUTTO2_RindPdig24.txt';
n=zeros(1,2);
k=zeros(1,2);
% geni contiene le stringhe di geni, il primo cell array si
% riferisce al primo dataset, il secondo si riferisce al secondo
% dataset. clustr contiene i relativi cluster rifeiti a ciascun gene di A
% double : converte una stringa nella codifica ASCII
% char: converte un codice ASCII in un carattere

cluster=cell(1,2);

for j=1:2

    fid=fopen(data_set{j},'r');
    % nella prima riga del dataset ho messo il numero delle osservazioni totale
    % e il numero dei geni distinti ricavati dal sito dell'ncbi. Nel resto del
    % file ho il numero di geni della libreira e vicino il nome del gene
    n(j)=fscanf(fid,'%f',1);
    k(j)=fscanf(fid,'%f',1);
    cluster{j}=zeros(1,k(j));
    % celle array di comodo per le sttringhe di geni
    STR=cell(1,k(j));

    % fid=fopen('filenam.txt','**'): apre il file identificandolo con fid nella
    % modalità indicata da **, lettura o scrittura di solito
    % fscanf: legge da file in un certo formato, per le stringhe meglio usare i
    % comandi successivi
    % fgets(fid): legge da file una riga del file e nella stringa è incluso
    % anche il termine della linea (ossia il comando di invio)
    % fgetl(fid): legge una linea del file e la trasforma in una stringa, ma
    % non include il comando di invio -> MI SERVE
    % nei due precedenti comandi, se si incontra l'eof ritorna -1: quindi mi
    % fermo a leggere il file quando la funzione fgetl ritorna -1

    for i=1:k(j)
        cluster{j}(i)=fscanf(fid,'%f',1);
        STR{i}= fscanf(fid,'%s',1);
    end
    fclose(fid);
    if j==1
        geni_FRUTTO1=STR;
    else
        geni_FRUTTO2=STR;
    end
    clear STR;
end


cluster_FRUTTO1=cluster{1};
cluster_FRUTTO2=cluster{2};
n1=n(1);
n2=n(2);
k1=k(1);
k2=k(2);
clear k n cluster;

% trasformo i geni in numeri
% creo un vetore indici che mi dice quli indici del secondo campione
% contengono geni uguali a quelli del primo campione
indici=[];
indici1=[];
X1_star=zeros(1,k1);
X2_star=zeros(1,k2);
for i=1:k1
    indice_corrente=[];
    for j=1:k2
        if length(geni_FRUTTO1{i})~=length(geni_FRUTTO2{j})
            continue;
        end
        if geni_FRUTTO2{j}==geni_FRUTTO1{i}
            indice_corrente=j;
            indici=[indici indice_corrente];
            indici1=[indici1 i];
        end
    end
    nome_gene=normrnd(M0,sqrt(V0));
    X1_star(i)=nome_gene;
    X2_star(indice_corrente)=nome_gene;
end
indici_rimasti=1:k2;
indici_rimasti(indici)=[];
X2_star(indici_rimasti)=normrnd(M0,sqrt(V0),1,length(indici_rimasti));

Xn1=[];
for i=1:k1
    Xn1=[Xn1 X1_star(i)*ones(1,cluster_FRUTTO1(i))];
end
Xn2=[];
for i=1:k2
    Xn2=[Xn2 X2_star(i)*ones(1,cluster_FRUTTO2(i))];
end  

% Informazioni sulla libreria

[X1_star n1_star]=clusterizza(Xn1);
[X2_star n2_star]=clusterizza(Xn2);
indici1=ismember(X1_star,X2_star);
indici2=ismember(X2_star,X1_star);
fprintf('I geni comuni sono %i : \n',sum(indici1));
fprintf('Nella prima libreria %i EST si riferiscono a tali geni \n',sum(n1_star(indici1==1)));
fprintf('Nella seconda libreria %i EST si riferiscono a tali geni \n',sum(n2_star(indici2==1)));


% genero altre 2000 osservazioni
m1=2000;
m2=2000;

[nuove_specie_totali X1_nuove_specie X2_nuove_specie M_X1 M_X2 M_T1 M_T2 M_parametri N]=posterior_K(m1,m2,Xn1,Xn2,a0,b0,a,b,M0,V0,var);

save clementina_FRUTTI.mat nuove_specie_totali X1_nuove_specie X2_nuove_specie M_X1 M_X2 M_T1 M_T2 M_parametri N n1 m1 n2 m2
