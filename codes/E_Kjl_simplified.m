function [val]=E_Kjl_simplified(p_j,nu,gamma,big_K,theta_j,sigma_j,b_0,m_j_dot,l)
    warning('off','all');
    C2 = generalized_factorial(l,l,sigma_j);
    val = 0;
    for i=1:l
        %s = F_nk(i,1,sigma_j,(theta_j+m_j_dot*sigma_j)*b_0,C2(i+1,1+1))*...
        %            pochhammer(gamma+(big_K+1)*nu,1)/...
        %            pochhammer(gamma+big_K*nu,1);
        s = F_nk(i,1,sigma_j,(theta_j+m_j_dot)*b_0,C2(i+1,1+1))*...
                    pochhammer(gamma+(big_K+1)*nu,1)/...
                    pochhammer(gamma+big_K*nu,1);
        val_curr = s;
        for m_hat=2:i
                %s = ((theta_j+m_j_dot*sigma_j)*b_0+(m_hat-1)*sigma_j)*...
                %    exp(C2(i+1,m_hat+1))*s*(gamma+(big_K+1)*nu+m_hat-1)/...
                %    ((gamma+big_K*nu+m_hat-1)*(exp(C2(i+1,m_hat))*sigma_j));
                s = ((theta_j+m_j_dot)*b_0+(m_hat-1)*sigma_j)*...
                    exp(C2(i+1,m_hat+1))*s*(gamma+(big_K+1)*nu+m_hat-1)/...
                    ((gamma+big_K*nu+m_hat-1)*(exp(C2(i+1,m_hat))*sigma_j));
                val_curr=val_curr+s;
        end
        sum_curr = log(nchoosek(l,i)) + i*log(p_j)+(l-i)*log(1-p_j)+log(val_curr - 1);
        val = val + exp(sum_curr);
        %val = val + nchoosek(l,i)*(p_j^i)*((1-p_j)^(l-i))*val_curr;
    end
    val = ((gamma+big_K*nu)/nu)*(val);
end