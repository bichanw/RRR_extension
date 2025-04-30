function [w0,urrr,vrrr] = svd_RRR_noniso(X, Y, rnk, C)

wls = (X'*X)\(X'*Y);  % least squares estimate

% calculate covariance of Y if there's no input
if nargin < 4
	res_wls = Y - X*wls; % residual for least square estimate
	% covariance for residuals
	C = res_wls'*res_wls/(size(X,1)-1);
end

% tmp = C^(-1/2);
C_sqrt = C^(1/2);
C_sqrt_inv = inv(C_sqrt);
[~,~,vrrr] = svd(C_sqrt_inv*Y'*X*wls*C_sqrt_inv);  % perform SVD of relevant matrix
vrrr = C_sqrt*vrrr(:,1:rnk);      % get column vectors before accounting for covaraince
urrr = (X'*X)\(X'*Y*inv(C)*vrrr);  % get row vectors
w0 = urrr*vrrr';  % construct full RRR estimate
