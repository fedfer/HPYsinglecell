% Function to compute the generalized factorial coefficients 

function c=loggfc(alpha,n)

if n<=0
    c=[];
    return
end
logalpha=log(alpha);
c=logalpha;
if n==1
    return;
end
for m=1:(n-1)
    x=zeros(1,m+1);
    a=log(m-alpha*(2:(m)))+c(2:(m));
    b=logalpha+c(1:(m-1));
    themax=max([a ; b]);
    x(1)=log(m-alpha)+c(1);
    x(2:(m))=themax+log(exp(a-themax)+exp(b-themax));
    x(m+1)=logalpha+c(m);
    c=x;
end