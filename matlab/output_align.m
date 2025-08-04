function outputalignmentidx = output_align(X,Y,W,if_plot)
% OUTPUT_ALIGN computes the alignment index of the output population
%   outputalignmentidx = output_align(X,Y,W)
%   X: input population (n_samples x n_input_features)
%   Y: output population (n_samples x n_output_features)
%   W: communication weights (n_input_features x n_output_features)
% Output:
%   outputalignmentidx: raw alignment index 

% have code here for random generation

if nargin < 4
    if_plot = false; % default is not to plot
end

% ----- calculate raw alignment index -----
% Do PCA on population
[upop,spopvec] = svd(cov(Y),'vector');
spopcum = cumsum(spopvec);
spopvec_nrm = spopvec/sum(spopvec);
spopcum_nrm = cumsum(spopvec_nrm);

% Project communication covariance onto PCs
cov_predicted = cov(X*W);
scomvec = diag(upop'*cov_predicted*upop);
scomvec_nrm = scomvec./sum(scomvec);
scomcum_nrm = cumsum(scomvec)/sum(scomvec);

% compute alignment index
muscom = mean(scomcum_nrm);
muspop = mean(spopcum_nrm);

alignment_raw = muscom-muspop;

% compute communication fraction
commfrac = sum(scomvec)/sum(spopvec);


% ----- normalize alignment index -----
% compute max possible alignment
totcom = sum(scomvec);
ii = find(spopcum>totcom+1e-10,1);  % find how many dimensions we'd need for maximally aligned communication
scommax = spopvec;
scommax(ii:end) = 0;
scommax(ii) = totcom-sum(scommax); % fix up final bin
scommax_cum = cumsum(scommax)/sum(scommax); % compute cumulative
a_max = mean(scommax_cum)-muspop; % max possible alignment score

% compute min possible alignment (by flipping order of pop eigenvalues)
spopvec_rev = flipud(spopvec);
spopcum_rev = cumsum(spopvec_rev);
ii = find(spopcum_rev>totcom+1e-10,1);  % find how many dimensions we'd need for minimally aligned communication
scommin = spopvec_rev;
scommin(ii:end) = 0;
scommin(ii) = totcom-sum(scommin); % fix up final bin
scommin = flipud(scommin); % flip back to normal ordering
scommin_cum = cumsum(scommin)/sum(scommin); % compute cumulative
a_min = mean(scommin_cum)-muspop; % min possible alignment score

% rescale alignment scores above and below zero
if alignment_raw>0
    outputalignmentidx = alignment_raw/a_max;
else
    outputalignmentidx = alignment_raw/abs(a_min);
end
% fprintf('alignment index=%.2f\n',outputalignmentidx);


if if_plot
% ====================================
% make plots
close all;
clrs = get(gca,'colororder');
nd = size(upop,1);  % number of output neurons

subplot(221);
plot(1:nd,spopvec_nrm,'o-',1:nd,scomvec/sum(spopvec),'*-');
ylabel('normalized variance'); box off;
set(gca,'ylim',[0 spopvec_nrm(1)*1.1]);
xlabel('PC dimension');
title('fraction of variance in PC dimensions');
legend('PC','communication','location', 'north');

subplot(223);
plot(1:nd,spopcum_nrm,'-o',1:nd,scomcum_nrm,'-*'); box off;
xlabel('PC dimension');
ylabel('normalized cumulative variance');
title('fraction of variance');
hold on;
plot([1 nd],muspop*[1 1], '--','color',clrs(1,:));
plot([1 nd], muscom*[1 1],'--','color',clrs(2,:));
hold off;
set(gca,'ylim',[0 1]);
legend('PCs', 'communiction', 'mean(PCs)', 'mean(communication)',...
    'location','south');

subplot(222);
plot(commfrac,outputalignmentidx,'k*'); 
hold on; 
plot([0 1],[0 0],'k--', [0.5 0.5], [-1 1], 'k--','linewidth',1); hold off;
%plot(outputalignment*[1 1],[0 1],'k', outputalignment,1,'k*');
set(gca,'ylim',[-1.0 1.01], 'xtick',-1:.5:1,'xlim',[0 1]);
axis square;
ylabel('output alignment index');
xlabel('communication fraction');
box off;
ef;
end