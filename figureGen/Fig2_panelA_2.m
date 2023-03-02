function Fig2a_2()

load('Fig2a2_table.mat')
ROIStimSp = Fig2a2_table.ControlProb_sp_syncStim;
ROIStimSpControlTemporal = Fig2a2_table.ROIStimSpControlTemporal;

figName = sprintf('Fig2_panela_TFR_roi_comparison_stim_vs_Sham');
f0 = figure('Name', figName,'NumberTitle','off');
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',8);
set(gcf,'DefaultAxesFontName','arial');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 4 6]); % this size is the maximal to fit on an A4 paper when printing to PDF
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
colormap('jet');

p4 = [0.2 0.1 0.7 0.8];
ii_p = 4;
eval(sprintf('axes(''position'',p%d)',ii_p))
[map,num,typ] = brewermap(565,'Spectral')
for ii = 1:length(ROIStimSpControlTemporal)
    plot([1+rand(1)*0.2,2+rand(1)*0.2], 100*[ROIStimSpControlTemporal(ii),ROIStimSp(ii)],'.r-','linewidth',0.01,'color',map(ii,:))
    hold all
end
axis([0.5 2.5 0 90])
set(gca,'ytick',[0 45 90],'xtick',[1,2],'xticklabel',[],'ticklength',[0.03 0.02])
box off

% calc stats
[p1,~] = signrank(ROIStimSp,ROIStimSpControlTemporal,'alpha',0.05)
