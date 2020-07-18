% function to compute E(K_j^l | Y_n,b_0,p_k)
function [val]=E_Kjl(p_j,l,T)
    val_mat = zeros(l,l);
    for k=1:l
        for i=k:l
            val_mat(k,i) = k*nchoosek(l,i)*(p_j^i)*((1-p_j)^(l-i))*...
                sum(squeeze(T(k,i,k:i))) ;
        end
        k
    end
    val = sum(sum(val_mat));
end