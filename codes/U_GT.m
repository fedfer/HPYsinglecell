function u = U_GT(f,t)

% Input:  ratio of future samples and past samples, t
%         fingerprint,f
% Output: GT estimator, if t \ge 1, then we compute the smoothed version
% with the poisson distribution
n = size(f,1);
T = zeros(n,1);
i = 1:n;
T(i, 1) = (-t).^i;
if t < 1
    u = - T'*f;
else 
    r = (1/2*t)*log((n*(t+1)^2)/(t-1));
    P = ones(1,n)-poisscdf(1:n,r);
    u = - T'.*P*f;
end
    u = max(u,0);
end