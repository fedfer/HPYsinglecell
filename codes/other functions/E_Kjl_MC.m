function [val]=E_Kjl_MC(p_j,alpha,gamma,big_K,theta_j,sigma_j,b_0,m_j_dot,l,S)
    C2 = generalized_factorial(l,l,sigma_j);
    probs = zeros(l,1);
    
    for m_hat = 0:l
        probs(l,1) = F_nk(i,m_hat,sigma_j,(theta_j+m_j_dot*sigma_j)*b_0,C2(i+1,m_hat+1));
    end
    
    for i=0:l
        val_curr = 0;
        for m_hat=0:i
                val_curr=val_curr+...
                    F_nk(i,m_hat,sigma_j,(theta_j+m_j_dot*sigma_j)*b_0,C2(i+1,m_hat+1))*...
                    pochhammer(gamma+(big_K+1)*alpha,m_hat)/...
                    pochhammer(gamma+big_K*alpha,m_hat);
        end
        val = val + nchoosek(l,i)*(p_j^i)*((1-p_j)^(l-i))*val_curr;
    end
    val = ((gamma+big_K*alpha)/alpha)*(val-1);
end