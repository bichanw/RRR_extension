function t = text_legend(ax,Legends,Colors,iCorner)
% function t = text_legend(ax,Legends,Colors,iCorner)
% write legend as texts of different colors at the bottom right of the figure


% set default color scheme
if nargin < 3 || isempty(Colors)
    Colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840];
    Colors = Colors(1:numel(Legends),:);
end

% default to bottom right corner
if nargin < 4
    iCorner = 4;
end

switch iCorner
case 1
t = arrayfun(@(ii) text(ax,ax.XLim(1), ax.YLim(2) - diff(ax.YLim) / 10 * (ii-1),...
                        Legends{end+1-ii},'Color',Colors(end+1-ii,:),'FontSize',ax.FontSize-1,'HorizontalAlignment','left','VerticalAlignment','top'), 1:numel(Legends));
case 2
t = arrayfun(@(ii) text(ax,ax.XLim(2), ax.YLim(2) - diff(ax.YLim) / 10 * (ii-1),...
                        Legends{end+1-ii},'Color',Colors(end+1-ii,:),'FontSize',ax.FontSize-1,'HorizontalAlignment','right','VerticalAlignment','top'), 1:numel(Legends));
case 3
t = arrayfun(@(ii) text(ax,ax.XLim(1), ax.YLim(1) + diff(ax.YLim) / 10 * (ii-1),...
                        Legends{end+1-ii},'Color',Colors(end+1-ii,:),'FontSize',ax.FontSize-1,'HorizontalAlignment','left','VerticalAlignment','bottom'), 1:numel(Legends));
case 4
t = arrayfun(@(ii) text(ax,ax.XLim(2), ax.YLim(1) + diff(ax.YLim) / 10 * (ii-1),...
                        Legends{end+1-ii},'Color',Colors(end+1-ii,:),'FontSize',ax.FontSize-1,'HorizontalAlignment','right','VerticalAlignment','bottom'), 1:numel(Legends));

end


end