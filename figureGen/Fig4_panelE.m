% Load data saved by analyzeHealthyPtsPerformance__plots.m
% XLS table generated by writeResults2XLS.m
folderName = 'E:\Data_p\ClosedLoopDataset\BEHAV_LINK';
visualPAL_fig_data = load(fullfile(folderName,'dPrimeChange_dataset'),'visualPAL_fig_data');
visualPAL_fig_data = visualPAL_fig_data.visualPAL_fig_data;

PtCellIdx = visualPAL_fig_data.PtCellIdx;
couple_mat = visualPAL_fig_data.couple_mat;
couple_mat_PRUNED = visualPAL_fig_data.couple_mat_PRUNED ;
dprime_change_sleep = visualPAL_fig_data.dprime_change_sleep;
dprime_change_stim = visualPAL_fig_data.dprime_change_stim;
dprime_diff = visualPAL_fig_data.dprime_diff;
N_rec_diff = visualPAL_fig_data.N_rec_diff;
N_falsePos_diff = visualPAL_fig_data.N_falsePos_diff;

% 3 subgroups based on STIM targeting and anatomical area
% subGroup{1} - frontal_well_targeted_stim_pts
% subGroup{2} - posterior_well_targeted_stim_pts
% subGroup{3} - posterior_neg_stim_pts
[subGroup, cmap] = defineSubGroups_stimAnatomy();
global G_targeted; G_targeted = 1;
global G_targeted_diff_anatomy; G_targeted_diff_anatomy = 2;
global G_mixed; G_mixed = 3;
stim_color = cmap(1,:);
SHAM_NIGHT_COLOR = cmap(4,:);

color_i = [];
for ii = 1:size(couple_mat_PRUNED,1)
    for sg_i = 1:3
        if ismember(couple_mat_PRUNED(ii,3),subGroup{sg_i})
            color_i(ii) = sg_i;
        end
    end
end


%%  Preparing behavioral Panel for figure 4
% Embedding all pts
filename = 'couplingAnalysis_forFig';
% mm = matfile('E:\Data_p\ClosedLoopDataset\couplingAnalysis\frontalStim_stimBothHemispheres\figures\couplingAnalysis_forFig.mat');
folderName = 'E:\Data_p\ClosedLoopDataset\couplingAnalysis\frontalStim_stimBothHemispheres_biPolarRipples_highFreqSpindles\';
% folderName = 'E:\Data_p\ClosedLoopDataset\couplingAnalysis\frontalStim_stimBothHemispheres_biPolarRipples_highFreqSpindles_90sec\'
mm = matfile(fullfile(folderName,filename) );
ac = AnalyzeCoupling;
ac.minRipplesPerCond = 20;
ac.minSW_TH = 0;

couplingAnalysisSummary = mm.couplingAnalysisSummary;
pt_num_AC = couplingAnalysisSummary.pt_list;
all_pts = unique(pt_num_AC);

mm = matfile(fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\','rippleInfo'));
rippleInfo = mm.rippleInfo;
hipRippleCh = zeros(1,length(couplingAnalysisSummary.channelCouples));
validRippleCh = zeros(1,length(couplingAnalysisSummary.channelCouples));
minRipplesPreSleep = 0;
neuralIndex_J = []; neuralIndex_otherHem_J = [];
for ii = 1:length(couplingAnalysisSummary.channelCouples)
    
    ind = find(couplingAnalysisSummary.channelCouples(ii,1) == rippleInfo.chNum_v  &...
        pt_num_AC(ii) == rippleInfo.ptNum_v );
    if isempty(ind)
        validRippleCh(ii) = false;
        continue
    else
        if ismember(ind,rippleInfo.artifactInd) | rippleInfo.nRip(ind) < minRipplesPreSleep
            validRippleCh(ii) = false;
        else
            validRippleCh(ii) = true;
        end
    end
    
    chArea = rippleInfo.area{ind};
    if chArea.isHip
        hipRippleCh(ii) = true;
    end
end
hipRippleCh = logical(hipRippleCh);
validRippleCh = logical(validRippleCh) ;


% 3 subgroups based on STIM targeting and anatomical area
[subGroup, cmap] = defineSubGroups_stimAnatomy();
AC_color_i = [];
for ii_p = 1:length(pt_num_AC)
    for sg_i = 1:3
        if ismember(pt_num_AC(ii_p),subGroup{sg_i})
            AC_color_i(ii_p) = sg_i;
        end
    end
end
AC_color_i(~validRippleCh) = 0;


%% Prep panel for Main Fig. 4e
figName = ('Figure4_behavPanels_e');
newA4figure(figName)

% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);

p_sz_x = 0.14;
p_sz_y = 0.18;
p_buf = 0.3;

p1 = [0.1,0.1,p_sz_x,p_sz_y];
p2 = [0.232,0.22,p_sz_x*.8,p_sz_y*.8];

PRE_IND = 6; SHORT_TERM_IND = 3  ;

for ii_c = 1:2
    
    clear neuralI_perPt str_t
    if ii_c == 1
        BASELINE = (couplingAnalysisSummary.EventRateCondsRipSWAll(validRippleCh,[PRE_IND]));
        STIM_EFFECT = (couplingAnalysisSummary.EventRateCondsRipSWAll(validRippleCh,SHORT_TERM_IND));
        str_t{1} = 'rip-sw coupling change';
    elseif ii_c == 2
        BASELINE = (couplingAnalysisSummary.EventRateCondsRipSWSpAll(validRippleCh,[PRE_IND]));
        STIM_EFFECT = (couplingAnalysisSummary.EventRateCondsRipSWSpAll(validRippleCh,SHORT_TERM_IND));
        str_t{1} = 'rip-sw-sp coupling change';
    end
    eval(sprintf('axes(''position'',p%d)',ii_c))
    
    % Single-event coupling
    DIFF_elec = (STIM_EFFECT-BASELINE)';
    nSW_per_conditionAll = couplingAnalysisSummary.nSW_per_conditionAll(validRippleCh,[SHORT_TERM_IND,PRE_IND])';
    nRip_per_conditionAll = couplingAnalysisSummary.nRip_per_conditionAll(validRippleCh,[SHORT_TERM_IND,PRE_IND])';
    % DIFF_elec(logical(sum(nSW_per_conditionAll < ac.minSW_TH ) | sum(nRip_per_conditionAll < ac.minRipplesPerCond ))) = NaN;
    DIFF_elec(logical((STIM_EFFECT == 0 ) & (BASELINE == 0))) = NaN;
    
    pt_num_AC = couplingAnalysisSummary.pt_list(validRippleCh);
    otherHemCouples = couplingAnalysisSummary.otherHemCouples(validRippleCh);
    neuralIndex = []; behavIndex = []; neuralIndex_otherHem = []; behavIndex_otherHem = []; pt_vec = []; pt_vec_oh = [];
    stimType_plot = []; stimType_plot2 = [];
    N_pairs_perPt = []; neuralI_sem_perPt = [];
    
    neuralI_perPt = [];
    
    for ii_p = 1:length(dprime_change_sleep)
        if ismember( couple_mat_PRUNED(ii_p,3),subGroup{2})
            neuralI_perPt(ii_p) = NaN;
            N_pairs_perPt(ii_p) = NaN;
            neuralI_sem_perPt(ii_p) = NaN;
            continue
        end
        indices  = find(couple_mat_PRUNED(ii_p,3) == pt_num_AC  & ~otherHemCouples );
        indices2 = find(couple_mat_PRUNED(ii_p,3) == pt_num_AC & otherHemCouples   );
        
        if (~isempty([indices,indices2]))
            allInd = [indices,indices2];
            neuralI_perPt(ii_p) = nanmedian(DIFF_elec(allInd));
            neuralI_sem_perPt(ii_p) = nanstd(DIFF_elec(allInd))/sqrt(length(allInd));
            if neuralI_perPt(ii_p) == 0
                disp(ii_p)
            end
        else % no contacts
            neuralI_perPt(ii_p) = NaN;
            N_pairs_perPt(ii_p) = NaN;
            neuralI_sem_perPt(ii_p) = NaN;
        end
        a = DIFF_elec(indices);
        
        ind = find(a == 0 | isnan(a));
        a(ind) = [];
        if ~isempty(a)
            b = visualPAL_fig_data.dprime_diff(ii_p) *ones(1,length(a));
            pp = couple_mat_PRUNED(ii_p,3) * ones(1,length(a));
            color_ii = color_i(ii_p) * ones(1,length(a));
            
            neuralIndex = [neuralIndex,a];
            behavIndex = [behavIndex,b];
            pt_vec = [pt_vec, pp];
            stimType_plot = [stimType_plot, color_ii];
        end
        c = DIFF_elec(indices2);
        ind = find(c == 0 | isnan(c));
        c(ind) = [];
        if ~isempty(c)
            d = visualPAL_fig_data.dprime_diff(ii_p) *ones(1,length(c));
            pp = couple_mat_PRUNED(ii_p,3) * ones(1,length(c));
            neuralIndex_otherHem = [neuralIndex_otherHem,c];
            behavIndex_otherHem = [behavIndex_otherHem,d];
            color_ii = color_i(ii_p)*ones(1,length(c));
            stimType_plot2 = [stimType_plot2, color_ii];
            pt_vec_oh = [pt_vec_oh, pp];
        end
        
    end
    sameHem = ones(length(neuralIndex),1);
    subj = pt_vec(:);
    behavI = behavIndex(:);
    neuralI = neuralIndex(:);
    stimType = stimType_plot(:);
    behav_neural_fig4_table1 = table(subj, behavI, neuralI, sameHem, stimType);
    subj = pt_vec_oh(:);
    behavI = behavIndex_otherHem(:);
    neuralI = neuralIndex_otherHem(:);
    stimType = stimType_plot2(:);
    sameHem = zeros(length(neuralIndex_otherHem),1);
    behav_neural_fig4_table2 = table(subj, behavI, neuralI, sameHem, stimType);
    behav_neural_fig4_table = [behav_neural_fig4_table1; behav_neural_fig4_table2];
    N_validpairs = height(behav_neural_fig4_table);
    
    rows = find(behav_neural_fig4_table.stimType == 3);
    
    % correaltion based on median values
    medV_non_nan = neuralI_perPt;
    dPrime_no_nan = visualPAL_fig_data.dprime_diff;
    rmvInd = isnan(medV_non_nan);
    medV_non_nan(rmvInd) = [];
    dPrime_no_nan(rmvInd) = [];
    [rho2,pval2] = corr(medV_non_nan',dPrime_no_nan','Type','Spearman');
    [r,pval_t] = corr(medV_non_nan',dPrime_no_nan');
    
    str_t{2} = sprintf('spearman/Pearson correlation:behav-coupling change');
    str_t{3} = sprintf('rho/r - %2.2e/%2.2e, P-%2.2e/%2.2e',r,rho2, pval2,pval_t);
    str_t{4} = sprintf('N = %d pairs/%d pts',N_validpairs, length(medV_non_nan));
    
    % removing mixed-phase pts
    ind = logical(color_i == 3);
    medV_non_nan_sync = neuralI_perPt ;
    dPrime_no_nan_sync = visualPAL_fig_data.dprime_diff;
    rmvInd = isnan(medV_non_nan_sync) | ind;
    medV_non_nan_sync(rmvInd) = [];
    dPrime_no_nan_sync(rmvInd) = [];
    [rho2_sync,pval2_sync] = corr(medV_non_nan_sync',dPrime_no_nan_sync','Type','Spearman');
    [r_sync,pval_t_sync] = corr(medV_non_nan_sync',dPrime_no_nan_sync');
    
    if ii_c == 1
        range = [-.6,.6,-0.4 0.4];
    else
        range = [-.8,.8,-0.4 0.4];
    end
    
    errorbar(neuralI_perPt,visualPAL_fig_data.dprime_diff,neuralI_sem_perPt,'horizontal','.k')
    hold all
    
    if SAVE_TABLE
        outputSrcFolder = 'C:\Users\mgeva\Documents\GitHub\closedLoop-pub\figureGen'
        a = neuralI_perPt;
        b = visualPAL_fig_data.dprime_diff;
        c = neuralI_sem_perPt;
        varNames  = {'neural_score_mean','behavioral_score','neural_score_sem'};
        data_table = table(a(:),b(:),c(:),'VariableNames',varNames);
        
        if ii_c == 1
            save(fullfile(outputFigureFolder2,'data_table'),'Fig4eMain_table');
        else
            save(fullfile(outputFigureFolder2,'data_table'),'Fig4eInset_table');
        end
    end
    
    for ii_ac = 1:length(neuralI_perPt)
        plot(neuralI_perPt(ii_ac), visualPAL_fig_data.dprime_diff(ii_ac),'.','markersize',14,...
            'color',cmap(color_i(ii_ac),:), 'linewidth',1.7)
        hold all
    end
    
    set(gca,'xlim',range(1:2),'ylim',range(3:4),...
        'xtick',[range(1),0,range(2)],'ytick',[-0.4:0.4:0.4],'TickDir','out')
    
    line(get(gca,'xlim'),zeros(1,2),'color','k')
    line(zeros(1,2),get(gca,'ylim'),'color','k')
    
    inds = ~isnan(neuralI_perPt);
    P = polyfit(neuralI_perPt(inds), visualPAL_fig_data.dprime_diff(inds),1);
    x0 = -1; x1 = 1;
    xi = linspace(x0,x1);
    yi = P(1)*xi+P(2);
    hold on
    plot(xi,yi,'k','linewidth',.7)
    XLIM = get(gca,'xlim');
    YLIM = get(gca,'ylim');
    text(XLIM(1),YLIM(2)+diff(YLIM),str_t)
    clear str_t
    
    neuralInd_all{ii_c} = neuralI_perPt;
    
end

aa = gcf;
set(aa,'renderer','zbuffer');
outputFigureFolder = 'E:\Data_p\ClosedLoopDataset\manuscriptFigures';
res = 400;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(aa.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!

outputFigureFolder = 'E:\Dropbox\Nir_Lab\closedLoopRevision\Figures\links_fig4\';
res = 400;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(aa.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!


