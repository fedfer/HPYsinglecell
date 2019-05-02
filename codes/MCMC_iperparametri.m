% aggiorno gli iperparametri dalle full conditionals secondo 4 MH

function  [ alpha, d ,gamma ,nu]=MCMC_iperparametri(mjk,m_j_dot,m_dd,J,nn,iter1,burnin1,alpha,d,gamma,nu,bigK)

M_alpha=zeros(iter1-burnin1,J);
M_d=zeros(iter1-burnin1,J);
M_nu=zeros(iter1-burnin1,1);
M_gamma=zeros(iter1-burnin1,1);

for y=1:iter1
    for j=1:J
        % aggiornamento dei parametri
        % aggiorno alpha
        alphaprop=exprnd(1);
        accept=gammaln((alphaprop/d(j))+m_j_dot(j))-gammaln(alphaprop+nn(j))-...
            gammaln((alpha(j)/d(j))+m_j_dot(j))+gammaln(alpha(j)+nn(j))+...
            gammaln(alpha(j)/d(j))-gammaln(alphaprop/d(j))+gammaln(alphaprop)-gammaln(alpha(j));
        accepttt=exp(accept);
        acceptprob=min([accepttt,1]);
        bern_accept=binornd(1,acceptprob);
        if (bern_accept==1)
            alpha(j)=alphaprop;
        end
        % aggiorno d
        dprop=unifrnd(0,1);
        loggfc_dprop=loggfc(dprop,nn(j));
        loggfc_d=loggfc(d(j),nn(j));
        accept=gammaln((alpha(j)/dprop)+m_j_dot(j))+loggfc_dprop(m_j_dot(j))-...
            gammaln((alpha(j)/d(j))+m_j_dot(j))-loggfc_d(m_j_dot(j))+...
            gammaln(alpha(j)/d(j))-gammaln(alpha(j)/dprop)+log((1-dprop)/(1-d(j)));
        accepttt=exp(accept);
        acceptprob=min([accepttt,1]);
        U=unifrnd(0,1);
        if U<=acceptprob
            d(j)=dprop;
        end
    end
    
    % aggiorno nu
    nuprop=unifrnd(0,1);
    loggfc_nuprop=loggfc(nuprop,m_dd);
    loggfc_nu=loggfc(nu,m_dd);
    accept=gammaln((gamma/nuprop)+bigK)+loggfc_nuprop(bigK)-gammaln(gamma/nuprop)...
        +log(1-nuprop)-gammaln((gamma/nu)+bigK)-loggfc_nu(bigK)...
        +gammaln(gamma/nu)+log(1-nu);
    accepttt=exp(accept);
    acceptprob=min([accepttt,1]);
    U=unifrnd(0,1);
    if U<=acceptprob
        nu=nuprop;
    end
    % aggiorna gamma
    gammaprop=exprnd(1);
    accept=gammaln((gammaprop/nu)+bigK)+gammaln(gammaprop)-gammaln(gammaprop/nu)-...
        gammaln(gammaprop+m_dd)-gammaln((gamma/nu)+bigK)-gammaln(gamma)+gammaln(gamma/nu)+gammaln(gamma+m_dd);
    accepttt=exp(accept);
    acceptprob=min([accepttt,1]);
    U=unifrnd(0,1);
    if U<=acceptprob
        gamma=gammaprop;
    end
    
    if y > burnin1
        M_alpha(y,:)=alpha;
        M_d(y,:)=d;
        M_nu(y)=nu;
        M_gamma(y)=gamma;
        
    end
end

N=iter1-burnin1;

alpha=sum(M_alpha)/N;
d=sum(M_d)/N;
gamma=sum(M_gamma)/N;
nu=sum(M_nu)/N;