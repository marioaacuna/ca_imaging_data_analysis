function c = dcorr(X, Y)
%DCORR computes the distance correlation between two random variables
% X and Y.
% Reference: http://en.wikipedia.org/wiki/Distance_correlation 
% 

X = X(:);
Y = Y(:);
n2 = numel(X) ^ 2;

a = pdist2mex(X', X', 'euc',[],[],[]);
b = pdist2mex(Y', Y', 'euc',[],[],[]);

A = a -bsxfun(@plus,mean(a,1),mean(a,2)) +mean(a(:));
B = b -bsxfun(@plus,mean(b,1),mean(b,2)) +mean(b(:));

dCovXY = sum(sum(A.*B)) /n2;
dVarX  = sum(sum(A.^2)) /n2;
dVarY  = sum(sum(B.^2)) /n2;
c = sqrt(dCovXY / sqrt(dVarX * dVarY));
if c<0, c=0; end
if c>1, c=1; end
if ~isfinite(c), c=0; end
if ~isreal(c), c=real(c); end
