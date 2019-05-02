function f=makeFinger(v)

% Input:  vector of integers, v
% Output: vector of fingerprints, f where f(i) = |{j: |{k:v(k)=j}|=i }|
%         i.e. f(i) is the number of elements that occur exactly i times 
%         in the vector v

%la prima fct va bene
h1 = hist(v,min(v):max(v));
%histcounts2(h1(:,1),h1(:,2))
%rows and columns start from zero 

f=hist(h1,0:max(h1));
f=f(2:end); f=f(:);