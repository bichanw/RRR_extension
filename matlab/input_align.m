function [aa,p,aa_rand] = input_align(X, W , r, C)
% calculate how much the communication weights W align with the
% principal components of the input X
% Input: 
%   X: input matrix (stimuli) of size (n_samples, n_input_neurons)
%   W: communication weights of size (n_input_neurons, n_communication_dims)
%   r: rank of the communication weights (optional, default is rank(W))
% Output:
%   aa: alignment index (0-1), where 1 is maximally aligned


% activity axes
if nargin < 4
    C = cov(X);
end
[Upca,Spcavec] = svd(C,"vector");
% Spcavec = diag(Spca);

% Compute SVD of weights
[Uw,Sw,Vw]= svd(W,"vector");
Swvec = [Sw; zeros(size(W,1)-size(W,2),1)];

% compute communication aligment index
amax = Spcavec'*(Swvec.^2);  % maximal value
amin = Spcavec'*(flipud(Swvec.^2)); % minimal value
araw = trace(W'*C*W);  % test statistic

aa = (araw-amin)/(amax-amin); 



if nargout > 1
    % ---------- run statistical test --------------
    % generate null distribuition

    [nx,ny] = size(W);
    if nargin < 3
        r = rank(W);
    end

    nperm = 1000;
    % generate random orthogonal matrix
    % and calculate alignment index
    nperm = 1000;
    aa_rand = zeros(nperm,1);
    upd = textprogressbar(nperm);
    for ii = 1:nperm
        [U,~] = qr(randn(size(W)));
        W_rand = U * diag(Swvec) * U';
        aa_rand(ii) = input_align(X, W_rand);
        upd(ii);
    end

    p = sum(aa_rand >= aa) / nperm;  % p-value

end


return

