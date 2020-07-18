% Compute the C factorial
function C=C_factorial(n,l,sigma)


if n==0 && l==0
    C=1;
    return;
elseif n<l
    C=0;
    return;
else
    log_vec=cumsum(log(1:l-1));
    log_vec=[log_vec(end:-1:1) 0];
    t_vec=1:l;
    segni=(-1).^t_vec;
    exp_vec=exp(gammaln(n-sigma*t_vec)-gammaln(-t_vec*sigma) ...
        -cumsum(log(t_vec))-log_vec);
    C=sum(segni.*exp_vec);
end
    