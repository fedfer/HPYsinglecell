function y=ponchamarrel(n,sigma)

if n==0
    y=1;
else
    y=exp(gammaln(sigma+n)-gammaln(sigma));
end