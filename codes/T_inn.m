%function to compute the tensor with dimensions (l,l,l) 
%for the indeces k,i,m_hat in order to
%compute the inner summation of E(K_j^l)
function [T]=T_inn(alpha,gamma,big_K,theta_j,sigma_j,b_0,m_j_dot,l)
    T = zeros(l,l,l);
    C1 = generalized_factorial(l,l,alpha);
    C2 = generalized_factorial(l,l,sigma_j);
    for k=1:l
        for i=k:l
            for m_hat=k:i
                T(k,i,m_hat)=F_nk(m_hat,k,alpha,gamma+big_K*alpha,C1(m_hat+1,k+1))*...
                    F_nk(i,m_hat,sigma_j,(theta_j+m_j_dot*sigma_j)*b_0,C2(i+1,m_hat+1));
            end
        end
        k
    end
end