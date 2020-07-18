% function to sample from a discrete distribution 

function y=gendiscr(supp,masses)

U=unifrnd(0,1);
Pcum=cumsum(masses);
index=find(U<=Pcum);
index=index(1);
y=supp(index);