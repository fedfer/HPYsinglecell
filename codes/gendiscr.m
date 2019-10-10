% function to sample from a discrete distribution 

function y=gendiscr(supp,masses)

U=unifrnd(0,1);
Pcum=cumsum(masses);
indice=find(U<=Pcum);
indice=indice(1);
y=supp(indice);