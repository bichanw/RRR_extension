function [w0,urrr,vrrr] = svd_RRR(X, Y, rnk, lambda)

% ridge regularzation
if nargin < 4
    lambda = 0;
end

XX = X'*X + lambda * eye(size(X,2));
if rcond(XX) < 1e-10
    wridge = pinv(XX)*(X'*Y);  % least squares estimate
else
    wridge = (XX)\(X'*Y);  % least squares estimate
end
[~,~,vrrr] = svd(Y'*X*wridge);  % perform SVD of relevant matrix


vrrr = vrrr(:,1:rnk);      % get column vectors
urrr= wridge * vrrr;  % get row vectors
w0 = urrr*vrrr';  % construct full RRR estimate
vrrr = vrrr';
