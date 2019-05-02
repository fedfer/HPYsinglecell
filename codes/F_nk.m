%function that computes of F(n,l,sigma,theta), which is the probability of
%observing k distinct values in (Y_1,....,Y_n) given the hyperparameters
%I take as input the logarighm of C(n,k,sigma) since I compute it beforehands 
%in order to save space
function [val]=F_nk(n,k,sigma,theta,log_C)
    val_log = 0;
    for r=1:(k-1)
        val_log = val_log+log((theta+r*sigma));
    end
    val_log = val_log - (k*log(sigma)+log(pochhammer(theta+1,n-1))) + log_C;
    val = exp(val_log);
end