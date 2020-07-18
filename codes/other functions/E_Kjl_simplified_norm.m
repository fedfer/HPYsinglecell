function [val]=E_Kjl_simplified_norm(p_j,alpha,gamma,big_K,theta_j,sigma_j,b_0,m_j_dot,l)
    % same as E_Kjl_simplified but using the normal approximation
    % N(np,np(1-p) to the binomial distribution
    C2 = generalized_factorial(l,l,sigma_j);
    val = 0;
    for i=0:l
        val_curr = 0;
        for m_hat=1:i
                val_curr=val_curr+...
                    F_nk(i,m_hat,sigma_j,(theta_j+m_j_dot*sigma_j)*b_0,C2(i+1,m_hat+1))*...
                    pochhammer(gamma+(big_K+1)*alpha,m_hat)/...
                    pochhammer(gamma+big_K*alpha,m_hat);
        end
        
        %sum = log(nchoosek(l,i)) + i*log(p_j)+(l-i)*log(1-p_j)+log(val_curr);
        %val = val + exp(sum);
        %val = val + nchoosek(l,i)*(p_j^i)*((1-p_j)^(l-i))*val_curr;
        val = val + normpdf(i,l*p_j,sqrt(l*(1-p_j)*p_j))*val_curr;
    end
    val = ((gamma+big_K*alpha)/alpha)*(val-1);
end