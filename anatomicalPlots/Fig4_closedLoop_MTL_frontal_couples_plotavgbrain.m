global globalFsDir;
globalFsDir='E:\Data_p\FreeSurferWinMirror';
% Preparing a figure with final couple view
% 1) removing couples that didn't pass ripple criteria
% 2) removing 489 that had too much IEDs to calculate TFRs

%%
outputFolder = 'E:\Data_p\ClosedLoopDataset\couplingAnalysis\frontalStim_stimBothHemispheres_biPolarRipples';
mm = matfile(fullfile(outputFolder,'coupleAreaStr.mat'));

patients_all = mm.pt_list_str;
%probeChan
probeChans_all = mm.MTL_area ;
%stim channels
stimChannel_all = mm.frontal_area;

cnt = 0; stimChannel = []; probeChans = [];
clear patients probeChans stimChannel
for ii = 1:length(patients_all)
        if strcmp(patients_all{ii},'505')
            continue
        end
    if isempty(probeChans_all{ii}) % this means it's an invalid couple (bad ripple channel)
        continue
    else
        cnt= cnt + 1;
        patients{cnt} = patients_all{ii};
        probeChans{cnt} = probeChans_all{ii};
        stimChannel{cnt} = stimChannel_all{ii};
        
        
%         if strcmpi('505',patients{cnt}) & strcmpi(stimChannel{cnt}, 'LOF1')
%             stimChannel{cnt} = 'LOF-AC1'; % pt 505
%         end


        if strcmpi('538',patients{cnt}) & strcmpi(stimChannel{cnt},'RHSG6')
            stimChannel{cnt} = 'RHSG5'; % pt 538 - use closest macro contact
        end
        
        if strcmpi(stimChannel{cnt},'RAC5')
            stimChannel{cnt} = 'ROPR5'; %
        end
        
        if strcmpi(stimChannel{cnt},'RAI2')
            stimChannel{cnt} = 'ROPR-AI2';
        end
        
    end
end

%stimChannel{12} = 'RSTG4';
%stimChannel{13} = 'RSTG4';

%stimChannel{30} = 'LOF-AC1';
%stimChannel{21} = 'RSTG3';


PgroupAvgCoords=[];
PgroupLabels=[];
PgroupIsLeft=[];
SgroupAvgCoords=[];
SgroupLabels=[];
SgroupIsLeft=[];
for ptIndex=1:length(patients)
    fprintf('Working on Participant %s\n',patients{ptIndex});
    
    
    if strcmpi(patients{ptIndex}, '505')
        continue % probe and channel are opposite hemispheres, cannot view
    end
    
    PATIENTS_INFO{ptIndex}.name = patients{ptIndex};
    
    elecNames_s = patients{ptIndex};
    
    
    subPath = fullfile(globalFsDir,char(PATIENTS_INFO{ptIndex}.name));
    elecReconPath=fullfile(subPath,'elec_recon');
    filename = fullfile(elecReconPath, sprintf('%sPostimpLoc.txt',patients{ptIndex}));
    [elec_name, elec_n, x, y, z, Hem, D] = textread(filename,'%s %d %f %f %f %s %s', 200);
    
    probe_ind = [];
    for ii = 1:length(elec_name)
        if strcmpi(probeChans{ptIndex}(1:end-1),elec_name{ii}) && ...
                strcmpi(probeChans{ptIndex}(end),num2str(elec_n(ii)))
            probe_ind = ii;
        end
    end
    if isempty(probe_ind); warning('probe ind missing'); end
    probe_elec_coord{ptIndex} = [x(probe_ind), y(probe_ind), z(probe_ind)];
    % Export probe x-y-z to an average brain coordnate system
    cfg=[];
    cfg.plotEm = 0;
    cfg.isSubdural=0; % 0 indicates that an electrode is a depth electrode
    cfg.elecCoord = probe_elec_coord{ptIndex};
    cfg.elecNames{1,1} = probeChans{ptIndex};
    cfg.isLeft = strcmpi(Hem{probe_ind},'L');
    
    [avgCoords, ELEC_NAMES, isLeft]=sub2AvgBrain(PATIENTS_INFO{ptIndex}.name,cfg);
    PgroupAvgCoords=[PgroupAvgCoords; avgCoords];
    PgroupLabels=[PgroupLabels, ELEC_NAMES];
    PgroupIsLeft=[PgroupIsLeft; isLeft];
    clear avgCoords ELEC_NAMES isLeft
    
    stim_ind = [];
    for ii = 1:length(elec_name)
        if strcmpi(stimChannel{ptIndex}(1:end-1),elec_name{ii}) && ...
                strcmpi(stimChannel{ptIndex}(end),num2str(elec_n(ii)))
            stim_ind = ii;
        end
    end
    if isempty(stim_ind); warning('stim ind missing'); end
    stim_elec_coord{ptIndex} = [x(stim_ind), y(stim_ind), z(stim_ind)];
    
    % Export stim x-y-z to an average brain coordnate system
    cfg=[];
    cfg.plotEm = 0;
    cfg.isSubdural=0;  % 0 indicates that an electrode is a depth electrode
    cfg.elecCoord = stim_elec_coord{ptIndex};
    cfg.elecNames{1,1} = stimChannel{ptIndex};
    cfg.isLeft = strcmpi(Hem{stim_ind},'L');
    
    [avgCoords, ELEC_NAMES, isLeft]=sub2AvgBrain(PATIENTS_INFO{ptIndex}.name,cfg);
    SgroupAvgCoords=[SgroupAvgCoords; avgCoords];
    SgroupLabels=[SgroupLabels, ELEC_NAMES];
    SgroupIsLeft=[SgroupIsLeft; isLeft];
    clear avgCoords ELEC_NAMES isLeft
    
end

anatomyInfo.patients = patients;
anatomyInfo.probeChans = probeChans;
anatomyInfo.stimChans = stimChannel;
anatomyInfo.SgroupAvgCoords = SgroupAvgCoords;
anatomyInfo.SgroupLabels = SgroupLabels;
anatomyInfo.SgroupIsLeft = SgroupIsLeft;
anatomyInfo.PgroupAvgCoords = PgroupAvgCoords;
anatomyInfo.PgroupLabels = PgroupLabels;
anatomyInfo.PgroupIsLeft = PgroupIsLeft;

save(fullfile(globalFsDir,'anatomyPlotsData','average_CouplingDB'),'anatomyInfo')

%% Plotting stim and probe on one brain (Main Fig 4)
anatomyInfo = load(fullfile(globalFsDir,'anatomyPlotsData','average_CouplingDB'),'anatomyInfo');
anatomyInfo = anatomyInfo.anatomyInfo;

close all

figName = sprintf('MTL_cortical_couples_locations');
f0 = figure('Name', figName,'NumberTitle','off');
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',8);
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


frontalElecColor = [0,0,0];
MTLElecColor = [0.1098,0.0980,0.8431];
elecColors = [repmat(MTLElecColor,length(anatomyInfo.PgroupLabels),1);...
    repmat(frontalElecColor,length(anatomyInfo.SgroupLabels),1)];

[subGroup, cmap] = defineSubGroups_stimAnatomy();
cmap(3,:) = [0.2,0.2,0.2];
for p_ct = 1:length(anatomyInfo.patients)
    pairs{p_ct,1} = anatomyInfo.PgroupLabels{p_ct};
    pairs{p_ct,2} = anatomyInfo.SgroupLabels{p_ct};
    
    for ii_g = 1:3
        if ismember(str2num(anatomyInfo.PgroupLabels{p_ct}(1:3)),subGroup{ii_g})
            pairs{p_ct,3} = cmap(ii_g,:);
        end
    end
    pairs{p_ct,4} = anatomyInfo.PgroupLabels{p_ct}(5);
end

cfg=[];
cfg.surfType = 'pial';
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 2*ones(1,length(elecColors));
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.elecCoord = [[anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];[anatomyInfo.SgroupAvgCoords,anatomyInfo.SgroupIsLeft]];
cfg.elecNames = [anatomyInfo.PgroupLabels,anatomyInfo.SgroupLabels];
% for ii = 1:length(anatomyInfo.PgroupLabels)
%     cfg.elecNames{ii} = num2str(ii);
%     cfg.elecNames{ii+length(anatomyInfo.PgroupLabels)} = num2str(ii);
% end
cfg.pairs = pairs;
cfg.lineWidth = 0.4;
cfg.edgeBlack = 'n';
cfg.axis = ax7;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='ri';
cfg.axis = ax8;
cfgOut=plotPialSurf('fsaverage',cfg);



outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');

a = gcf;
set(a,'renderer','zbuffer');
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!

%% Prep a brain - only with frontal electrodes, for extended data figure 8 (revision)

figName = sprintf('corticalContacts_couples_locations');
f0 = figure('Name', figName,'NumberTitle','off');
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',8);
set(gcf,'DefaultAxesFontName','arial');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 19 24.7]); % this size is the maximal to fit on an A4 paper when printing to PDF
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
colormap('jet');

ax7 = axes('position',[0.1,0.7,.2,.2],'units','centimeters');
ax8 = axes('position',[0.4,0.7,.2,.2],'units','centimeters');


frontalElecColor = [0,0,0];

cfg=[];
cfg.surfType = 'pial';
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 2*ones(1,length(anatomyInfo.SgroupLabels));
cfg.elecColors = repmat(frontalElecColor,length(anatomyInfo.SgroupLabels),1);
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.elecCoord = [anatomyInfo.SgroupAvgCoords,anatomyInfo.SgroupIsLeft];
cfg.elecNames = anatomyInfo.SgroupLabels;
% for ii = 1:length(anatomyInfo.PgroupLabels)
%     cfg.elecNames{ii} = num2str(ii);
%     cfg.elecNames{ii+length(anatomyInfo.PgroupLabels)} = num2str(ii);
% end
cfg.pairs = [];
cfg.lineWidth = 0.4;
cfg.edgeBlack = 'n';
cfg.axis = ax7;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='ri';
cfg.axis = ax8;
cfgOut=plotPialSurf('fsaverage',cfg);

outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');

a = gcf;
set(a,'renderer','zbuffer');
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!




%% Other options

elecColors = repmat([0,0,1],length(anatomyInfo.PgroupLabels),1);
% brain
cfg=[];
cfg.view='l';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 2;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
cfg.elecNames = anatomyInfo.PgroupLabels;
cfg.axis = ax1;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='r';
cfg.axis = ax2;
cfgOut=plotPialSurf('fsaverage',cfg);


elecColors = repmat([1,0,0],length(anatomyInfo.PgroupLabels),1);
% brain
cfg=[];
cfg.view='l';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 2;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.elecCoord = [anatomyInfo.SgroupAvgCoords,anatomyInfo.SgroupIsLeft];
cfg.elecNames = anatomyInfo.SgroupLabelnss;
cfg.axis = ax3;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='r';
cfg.axis = ax4;
cfgOut=plotPialSurf('fsaverage',cfg);

frontalElecColor = [0,0,0];
MTLElecColor = [0.1098,0.0980,0.8431];
elecColors = [repmat(MTLElecColor,length(anatomyInfo.PgroupLabels),1);repmat(frontalElecColor,length(anatomyInfo.SgroupLabels),1)];
cfg=[];
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 2;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.elecCoord = [[anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];[anatomyInfo.SgroupAvgCoords,anatomyInfo.SgroupIsLeft]];
cfg.elecNames = [anatomyInfo.PgroupLabels,anatomyInfo.SgroupLabels];
cfg.axis = ax5;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='ri';
cfg.axis = ax6;
cfgOut=plotPialSurf('fsaverage',cfg);

[subGroup, cmap] = defineSubGroups_stimAnatomy();
cmap(3,:) = [0.2,0.2,0.2];
for p_ct = 1:length(anatomyInfo.patients)
    pairs{p_ct,1} = anatomyInfo.PgroupLabels{p_ct};
    pairs{p_ct,2} = anatomyInfo.SgroupLabels{p_ct};
    
    for ii_g = 1:3
        if ismember(str2num(anatomyInfo.PgroupLabels{p_ct}(1:3)),subGroup{ii_g})
            pairs{p_ct,3} = cmap(ii_g,:);
        end
    end
    pairs{p_ct,4} = anatomyInfo.PgroupLabels{p_ct}(5);
end


cfg=[];
cfg.surfType = 'pial';
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 2;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.elecCoord = [[anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];[anatomyInfo.SgroupAvgCoords,anatomyInfo.SgroupIsLeft]];
cfg.elecNames = [anatomyInfo.PgroupLabels,anatomyInfo.SgroupLabels];
% for ii = 1:length(anatomyInfo.PgroupLabels)
%     cfg.elecNames{ii} = num2str(ii);
%     cfg.elecNames{ii+length(anatomyInfo.PgroupLabels)} = num2str(ii);
% end
cfg.pairs = pairs;
cfg.lineWidth = 0.4;
cfg.edgeBlack = 'n';
cfg.axis = ax7;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='ri';
cfg.axis = ax8;
cfgOut=plotPialSurf('fsaverage',cfg);

coupleMap = brewermap(length(anatomyInfo.PgroupLabels),'Spectral');
cfg=[];
cfg.surfType = 'pial';
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 2;
cfg.elecColors = [coupleMap; coupleMap];
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.elecCoord = [[anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];[anatomyInfo.SgroupAvgCoords,anatomyInfo.SgroupIsLeft]];
cfg.elecNames = [anatomyInfo.PgroupLabels,anatomyInfo.SgroupLabels];
% for ii = 1:length(anatomyInfo.PgroupLabels)
%     cfg.elecNames{ii} = num2str(ii);
%     cfg.elecNames{ii+length(anatomyInfo.PgroupLabels)} = num2str(ii);
% end
% cfg.pairs = pairs;
cfg.lineWidth = 0.4;
cfg.edgeBlack = 'n';
cfg.axis = ax5;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='ri';
cfg.axis = ax6;
cfgOut=plotPialSurf('fsaverage',cfg);


outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');

a = gcf;
set(a,'renderer','zbuffer');
res = 400;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
res = 400;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
