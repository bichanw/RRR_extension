function [X,Y,U,V,ops] = simu_RRR(ops)
    % simplest case: Y = XB
    if nargin == 0
        ops = struct;
    end
    [T,ops] = getOr(ops,'T',100);
    [nx,ops] = getOr(ops,'nx',10);
    [ny,ops] = getOr(ops,'ny',5);
    [rnk,ops] = getOr(ops,'rnk',1);
    [signse,ops] = getOr(ops,'signse',0.1);
    [magnitude,ops] = getOr(ops,'magnitude',1);
    [thetas,ops] = getOr(ops,'thetas',[]);
    [Sigma,ops] = getOr(ops,'Sigma',[]);
    % input U, V if not wish to generate a new pair
    [U,ops] = getOr(ops,'U',[]);
    [V,ops] = getOr(ops,'V',[]);
    

    % generate input
    X = randn(T,nx);
    X = X - mean(X,1); % add an extra step to zero data?

    % generate communication matrix
    if isempty(U)
        U = randn(nx,rnk);
        V = randn(ny,rnk);
        if ~isempty(thetas)
            [U,~,V] = svd(U*V','econ');
            U = U(:,1:rnk);
            V = V(:,1:rnk)./sqrt(thetas);
            % [~,S,~] = svd(U*V','vector')
        end
        V = V * magnitude;
    end
    B = U*V';

    % generate Y
    if isempty(Sigma)
        E = signse*randn(T,ny);
        E = E - mean(E,1);
        Y = X * B + E;
    else
        Y = X * B + mvnrnd(zeros(1,ny),Sigma,T);
    end
    
end

