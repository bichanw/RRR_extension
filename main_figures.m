% code to generate figure 2



% fig 2a: advantage of regularization 

    % amount of data, estimation error w/ or w/o regularization
    n_sim = 50;  
    nx = 50;ny = 50;
    NSamples = round(50*2.^(1:.5:4));
    % Lambdas = 0:0.1:2;
    Lambdas = [0:10 20 100 1e3];
    err_fun = @(A, B) mean((A - B).^2, 'all') / mean(B.^2, 'all') ;

    U = [];V = [];
    err = nan(n_sim,numel(Lambdas),numel(NSamples));
    for isim = 1:n_sim
        [X,Y,U,V,ops] = simu_RRR(struct('rnk',50,'nx',nx,'ny',ny,'signse',50,'T',max(NSamples),'U',U,'V',V));
        w_true = U*V';
        for jj = 1:numel(Lambdas)
            w_hat = arrayfun(@(n_sample) svd_RRR(X(1:n_sample,:), Y(1:n_sample,:), ops.rnk,Lambdas(jj)),NSamples,'uni',0);
            err(isim,jj,:) = cellfun(@(w_) err_fun(w_,w_true), w_hat);
        end
    end


    % identify best lambda
    tmp = squeeze(mean(err,1));
    [~,I] = min(tmp);
    % plotting 
    close all;
    ax = subplot(1,1,1,'NextPlot','add');
    % plot best & no lambda
    plot_multiple_lines(permute(err(:,[1 I(3)],:),[1 3 2]),ax(1),'to_plt',{'err'},'x',NSamples);
    text_legend(ax(1),{'w/o reg','w/ reg'},[],2);ylabel(ax(1),'estimation error');xlabel(ax(1),'# samples');
    % plot all lambdas
    set(findall(gcf,'-property','FontSize'),'FontSize',10);
    set(ax,'FontSize',9,'XLim',[min(NSamples) max(NSamples)]);
    set(gcf,'Position',[0 0 3 2.2]*72);
    saveas(gcf,'fig2a','png');
    return



% fig 2b: advantage of non-isotropic noise
    % define error function
    err_fun = @(A, B) mean((A - B).^2, 'all') / mean(B.^2, 'all');

    % simulate with RRR
    % maximum variance 
    M_sigma = linspace(1,500,6);

    % run simulation
    n_sim = 50;
    ny = 50; nx=50;
    err = nan(n_sim,numel(M_sigma),2);
    for ii = 1:numel(M_sigma)
        % generate ground truth covariance
        Sigma = diag(logspace(0,log10(M_sigma(ii)),ny)); % variance on a log scale
        % Sigma = diag(linspace(1,M_sigma(ii),ny)); % variance on a linear scale
        % Sigma = diag([ones(1,ny-1) M_sigma(ii)]); % variance all ones but the last 

        for isim = 1:n_sim
            [X,Y,U,V,ops] = simu_RRR(struct('Sigma',Sigma,'T',100,'nx',nx,'ny',ny,'signse',50));

            % estimate, cross validation
            w0_noniso = svd_RRR_noniso(X, Y, ops.rnk, Sigma);
            w0 = svd_RRR(X, Y, ops.rnk);
            err(isim,ii,:) = cellfun(@(x) err_fun(x,U*V'),{w0,w0_noniso});
        end
    end

    % plot
    close all;
    Colors = colormap('lines');
    ax = subplot(1,1,1,'NextPlot','add');
    for ii = 1:2
        plot_multiple_lines(err(:,:,ii),ax,'x',M_sigma,'base_color',Colors(ii,:),'to_plt',{'err'});
    end
    xlabel('max(\sigma^2)','Interpreter','tex');
    ylabel('estimation error');
    t = text_legend(ax,{'isotropic','non-iso'},[],1);
    set(findall(gcf,'-property','FontSize'),'FontSize',10);
    set(ax,'FontSize',9,'XLim',[min(M_sigma) max(M_sigma)]);
    set(gcf,'Position',[0 0 3 2.2]*72);
    saveas(gcf,'fig2b','png');

    return



