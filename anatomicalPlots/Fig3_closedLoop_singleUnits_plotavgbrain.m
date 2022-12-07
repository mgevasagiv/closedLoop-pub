global globalFsDir;
globalFsDir='E:\Data_p\FreeSurferWinMirror';

patients = {'485','486',  '487',  '488', '489', '490', '496', '497', '498', '499', '505',   '510',  '515','520','538','541','544','545'};
PgroupAvgCoords=[];
PgroupLabels=[];
PgroupIsLeft=[];
SgroupAvgCoords=[];
SgroupLabels=[];
SgroupIsLeft=[];

[subGroup, cmap]  = defineSubGroups_stimAnatomy();

%% Load SW/spindle indices
summaryFile = 'E:\Data_p\ClosedLoopDataset\SWTriggeredSpikeHistResults\SWTriggeredSpikeHist_allSes_AllPts_stableUnits_summaryForPlot.mat';
mm = matfile(summaryFile);
ALL_UNITS = mm.ALL_UNITS;
runData = mm.runData;

% Create appropriate matrices to feed into display function (contact name
% and location):
missingContacts = [];
probe_elec_coord = cell(1,1);
contactID = [];     
for iContact = 1:length(ALL_UNITS.pt_list)
    disp(iContact)
    currPt = num2str(ALL_UNITS.pt_list(iContact));
    if ismember(currPt, subGroup{3})
        continue
    end
    
    % Load new data if pt num changed 
    if iContact == 1 || ALL_UNITS.pt_list(iContact) ~= ALL_UNITS.pt_list(iContact-1)
        disp('switching pt')
        disp(currPt)
        subPath = fullfile(globalFsDir,char(currPt));
        elecReconPath=fullfile(subPath,'elec_recon');
        filename = fullfile(elecReconPath, sprintf('%sPostimpLoc.txt',char(currPt)));
        [elec_name, elec_n, x, y, z, Hem, D] = textread(filename,'%s %d %f %f %f %s %s', 200);
        
        % pt index in runData (needed for MACRO file link)
        for iiP = 1:length(runData)
            if strcmp(runData(iiP).patientName(2:end),currPt)
                iPatient = iiP;
            end
        end
    end
    
    % Get the full contact name for each contact
    try
        % load macro montage for area name
        macroMontage = load(runData(iPatient).macroMontageFileName);
        macroMontage = macroMontage.MacroMontage;
        unitArea = ALL_UNITS.area_list(iContact);
        ChNum = nan;
        for jj = 1:length(macroMontage)
            if strcmpi(macroMontage(jj).Area,unitArea) % first channel of this area
                ChNum = jj;
                break
            end
        end
        if isnan(ChNum); error('channel not found'); end
        ChLabel = [macroMontage(jj).Area,'1'];
    catch
        disp('probe area-name doesn''t match options')
    end
    
    contact_ind = [];
    for ii = 1:length(elec_name)
        if strcmpi(ChLabel(1:end-1),elec_name{ii}) && ...
                strcmpi(ChLabel(end),num2str(elec_n(ii)))
            contact_ind = ii;
            break
        end
    end
    if ~isempty(contact_ind); disp(elec_name{contact_ind}); end
    if isempty(contact_ind); warning('contact missing in mloc file');
        missingContacts = [missingContacts, iContact];
        continue
    else
        elec_coord_pt_space = [x(contact_ind), y(contact_ind), z(contact_ind)];
        cfg=[];
        cfg.plotEm = 0;
        cfg.isSubdural=0; % 0 indicates that an electrode is a depth electrode
        cfg.elecCoord = elec_coord_pt_space;
        cfg.elecNames{1,1} = ChLabel;
        cfg.isLeft = strcmpi(ChLabel(1),'L');
        
        [avgCoords, ELEC_NAMES, isLeft]=sub2AvgBrain(currPt,cfg);
        
        PgroupAvgCoords=[PgroupAvgCoords; avgCoords];
        PgroupLabels=[PgroupLabels, ELEC_NAMES];
        PgroupIsLeft=[PgroupIsLeft; isLeft];
        contactID = [contactID, iContact];
        clear avgCoords ELEC_NAMES isLeft
    end
end

% Display effect for targeted stimulation patients

if isempty(missingContacts)
    anatomyInfo.patients =  ALL_UNITS.pt_list;
    anatomyInfo.isProbeHem =  ALL_UNITS.isProbeHem;
    anatomyInfo.FiringRates =  ALL_UNITS.FiringRates;
    anatomyInfo.mu_e = ALL_UNITS.mu_e;
    anatomyInfo.r_e = ALL_UNITS.r_e;
    anatomyInfo.pVal_mu = ALL_UNITS.pVal_mu;
end

anatomyInfo.PgroupAvgCoords = PgroupAvgCoords;
anatomyInfo.PgroupLabels = PgroupLabels;
anatomyInfo.PgroupIsLeft = PgroupIsLeft;

save(fullfile(globalFsDir,'anatomyPlotsData','spikingUnitsPhaseLocking'),'anatomyInfo')

%% Plot on MNI brain

anatomyInfo = load(fullfile(globalFsDir,'anatomyPlotsData','spikingUnitsPhaseLocking'),'anatomyInfo');
anatomyInfo = anatomyInfo.anatomyInfo;
pt_vec = anatomyInfo.patients;
r_e = anatomyInfo.r_e';
BASELINE = 1;
STIM = 3;
pval_th = 0.05;
DIFF =  (r_e(:,STIM)-r_e(:,BASELINE))./(r_e(:,STIM)+r_e(:,BASELINE));
ids = (anatomyInfo.pVal_mu(:,BASELINE) > pval_th) | (anatomyInfo.pVal_mu(:,STIM) > pval_th);
DIFF(ids) = nan;
DIFF((r_e(:,STIM)==0) | isnan(r_e(:,STIM)) | (r_e(:,BASELINE)==0) | isnan(r_e(:,BASELINE)) ) = nan;



frontalVec = logical(zeros(1,length(pt_vec)));
for ii = 1:length(pt_vec)
    area = classifyArea(anatomyInfo.PgroupLabels{ii}(5:end-1));
    if area.isFrontal
        frontalVec(ii) = true;
    end
end
[p,h,stats] = ranksum(DIFF(frontalVec),DIFF(~frontalVec));

%%


% Reorganize the plotted electrodes
% First remove NAN entries
DIFF_ORIG = DIFF;
rmv_ind = find(isnan(DIFF_ORIG));
PgroupAvgCoords = anatomyInfo.PgroupAvgCoords;
PgroupIsLeft = anatomyInfo.PgroupIsLeft;
PgroupLabels = anatomyInfo.PgroupLabels;
PgroupAvgCoords(rmv_ind,:) = [];
PgroupIsLeft(rmv_ind) = [];
PgroupLabels(rmv_ind) = [];
DIFF(rmv_ind) = [];

% Now merge values in each area 
anatomyInfo_plot.PgroupAvgCoords = PgroupAvgCoords;
anatomyInfo_plot.PgroupIsLeft = PgroupIsLeft;
anatomyInfo_plot.PgroupLabels = PgroupLabels;
clear median_DIFF markerSizeCnt DIFF_tmp_vec
rmv_ind = []; markerSizeCnt(1) = 1; DIFF_tmp_vec = DIFF(1);
cnt = 1;
median_DIFF(1) = DIFF(1); global_cnt = 1;
for ii = 2:length(PgroupLabels) 
    if PgroupAvgCoords(ii) == PgroupAvgCoords(ii-1)
        % merge entries
        DIFF_tmp_vec(end+1) = DIFF(ii);
        cnt = cnt + 1;
        rmv_ind(end+1) = ii;
        markerSizeCnt(ii) = cnt;
        markerSizeCnt(ii-1) = 0;
    else
        global_cnt = global_cnt + 1;
        if length(DIFF_tmp_vec)>1
            median_DIFF(global_cnt) = nanmedian(DIFF_tmp_vec);
        else
           median_DIFF(global_cnt) =  DIFF(ii);
           DIFF_tmp_vec = DIFF(ii);
        end
        markerSizeCnt(ii) = 1;
        cnt = 1;
    end
end
markersize = markerSizeCnt(markerSizeCnt~=0);
anatomyInfo_plot.PgroupAvgCoords(rmv_ind,:) = [];
anatomyInfo_plot.PgroupIsLeft(rmv_ind) = [];
anatomyInfo_plot.PgroupLabels(rmv_ind) = [];

N = 64;
% map = colormap('jet');
map = brewermap(N,'Spectral');
map = flipud(map);
values = linspace(-.25,.25,N);
[colorBin,~] = discretize(median_DIFF,values);
colorBin(find(median_DIFF < values(1))) = 1;
colorBin(find(median_DIFF > values(end))) = N;
elecColors = map(colorBin,:);
Nelec = size(elecColors,1);

%%
figName = sprintf('stimEffects_locations_singleUnits');
f0 = figure('Name', figName,'NumberTitle','off');
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','arial');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 19 24.7]); % this size is the maximal to fit on an A4 paper when printing to PDF
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
colormap('jet');

ax1 = axes('position',[0.1,0.1,.2,.2],'units','centimeters');
ax2 = axes('position',[0.4,0.1,.2,.2],'units','centimeters');
ax3 = axes('position',[0.1,0.3,.2,.2],'units','centimeters');
ax4 = axes('position',[0.4,0.3,.2,.2],'units','centimeters');
ax5 = axes('position',[0.1,0.5,.2,.2],'units','centimeters');
ax6 = axes('position',[0.4,0.5,.2,.2],'units','centimeters');
ax7 = axes('position',[0.1,0.7,.2,.2],'units','centimeters');
ax8 = axes('position',[0.4,0.7,.2,.2],'units','centimeters');


cfg=[]; 
% cfg.surfType='inflated';
cfg.ignoreDepthElec='n';
cfg.elecShape='sphere';
cfg.opaqueness=0.4;
cfg.elecSize = 0.7*markersize';
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.edgeBlack = 'n';
cfg.elecCoord = [anatomyInfo_plot.PgroupAvgCoords,anatomyInfo_plot.PgroupIsLeft];  
cfg.elecNames = anatomyInfo_plot.PgroupLabels;
cfg.pullOut = 0;
cfg.axis = ax7;
cfg.view='li';
cfgOut=plotPialSurf('fsaverage',cfg);
%view(ax7,[-94,27])

cfg.view='ri';
cfg.axis = ax8;
cfgOut=plotPialSurf('fsaverage',cfg);
%view(ax8,[94,27])

cfg.axis = ax5;
cfg.view='l';
cfgOut=plotPialSurf('fsaverage',cfg);
%view(ax5,[-94,27])

cfg.view='r';
cfg.axis = ax6;
cfgOut=plotPialSurf('fsaverage',cfg);
%view(ax6,[94,27])

outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');

a = gcf;
set(f0,'renderer','zbuffer');
% res = 400;
% eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!

figure
figName = 'rendererColorbar_singleU';
imagesc(1*ones(1,64),1:64,values)
colormap(map)
cc = colorbar;
cc.TickDirection = 'out';
cc.Ticks = [-.25,0,.25];
title(sprintf('axis limits- [%f,%f]',values(1),values(end)))
a = gcf;
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!

