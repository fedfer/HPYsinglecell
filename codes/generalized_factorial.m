% Function to generate a matrix with logs of generalized factorial coefficients, using
% the recurrency relation. 
% Instead of C(n,l), we will have the log of C(n-1,l-1,s) since we wanto to have 
% the coefficients with n=0 and l = 0


function LogC=generalized_factorial(N,L,s)

LogC=zeros(N+1,L+1);

LogC(1,1)=0;

for n=1:N
    LogC(n+1,n+1)=n*log(s);
end
for n=2:N
    LogC(n+1,2)=log(n-1-s)+LogC(n,2);
end
for n=2:N
    for l=3:n
 
        x=log(s+exp(log(n-1-s*(l-1))+LogC(n,l)-LogC(n,l-1)));
        if x==inf

            x=log(n-1-s*(l-1))+LogC(n,l)-LogC(n,l-1);
            disp('Approssimazione di x');
        end
        LogC(n+1,l)= LogC(n,l-1)+x;
    end
end
 
