global globalFsDir;
globalFsDir='E:\Data_p\FreeSurferWinMirror';

patients = {'485','486',  '487',  '488', '489', '490', '496', '497', '498', '499', '505',   '510',  '515','520','538','541','544','545'};
%probeChan
probeChans = {'RMH1', 'LEC1', 'RAH1','LPHG2','RAH1','RAH1','RMH1','RMH2','LAH1','LMH1','RAH1','RAH1','RAH1','RMH1','RA1','LMH2','RMH1','REC1'};
%stim channels
stimChannel = {'RpSMA4','LOF4','ROF5','LOF4','RPHG5','RANH5','RTOP7','ROF7','LOF6','LOF6','LOF-AC7','ROF6','ROF2','ROF5','RHSG4','LOF4','ROPR-AI1','ROF4'};

PgroupAvgCoords=[];
PgroupLabels=[];
PgroupIsLeft=[];
SgroupAvgCoords=[];
SgroupLabels=[];
SgroupIsLeft=[];
for ptIndex=1:length(patients)
            fprintf('Working on Participant %s\n',patients{ptIndex});
        
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
        if isempty(probe_ind); warning('stim ind missing'); end
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

save(fullfile(globalFsDir,'anatomyPlotsData','average_PROBE_STIM'),'anatomyInfo')

%% Plotting stim and probe on one brain
anatomyInfo = load(fullfile(globalFsDir,'anatomyPlotsData','average_PROBE_STIM'),'anatomyInfo');
anatomyInfo = anatomyInfo.anatomyInfo;

figName = sprintf('STIM_PROBE_locations');
f0 = figure('Name', figName,'NumberTitle','off');
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','arial');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 19 24.7]); % this size is the maximal to fit on an A4 paper when printing to PDF
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
colormap('jet');

% ax1 = axes('position',[0.1,0.1,.2,.2],'units','centimeters');
% ax2 = axes('position',[0.4,0.1,.2,.2],'units','centimeters');
% ax3 = axes('position',[0.1,0.3,.2,.2],'units','centimeters');
% ax4 = axes('position',[0.4,0.3,.2,.2],'units','centimeters');
% ax5 = axes('position',[0.1,0.5,.2,.2],'units','centimeters');
% ax6 = axes('position',[0.4,0.5,.2,.2],'units','centimeters');
ax7 = axes('position',[0.1,0.7,.2,.2],'units','centimeters');
ax8 = axes('position',[0.4,0.7,.2,.2],'units','centimeters');

% elecColors = repmat([0,0,1],length(anatomyInfo.PgroupLabels),1);
% % brain
% cfg=[]; 
% cfg.view='l';
% cfg.ignoreDepthElec='n';
% cfg.opaqueness=0.3;
% cfg.elecSize = 2;
% cfg.elecColors = elecColors;
% cfg.elecColorScale = [0 64];
% cfg.showLabels='n';
% cfg.title= '';
% cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];  
% cfg.elecNames = anatomyInfo.PgroupLabels;
% cfg.axis = ax1;
% cfgOut=plotPialSurf('fsaverage',cfg);
% cfg.view='r';
% cfg.axis = ax2;
% cfgOut=plotPialSurf('fsaverage',cfg);
% 
% 
% elecColors = repmat([1,0,0],length(anatomyInfo.PgroupLabels),1);
% % brain
% cfg=[]; 
% cfg.view='l';
% cfg.ignoreDepthElec='n';
% cfg.opaqueness=0.3;
% cfg.elecSize = 2;
% cfg.elecColors = elecColors;
% cfg.elecColorScale = [0 64];
% cfg.showLabels='n';
% cfg.title= '';
% cfg.elecCoord = [anatomyInfo.SgroupAvgCoords,anatomyInfo.SgroupIsLeft];  
% cfg.elecNames = anatomyInfo.SgroupLabels;
% cfg.axis = ax3;
% cfgOut=plotPialSurf('fsaverage',cfg);
% cfg.view='r';
% cfg.axis = ax4;
% cfgOut=plotPialSurf('fsaverage',cfg);
% 
% 
% elecColors = [repmat([0.1098,0.0980,0.8431],length(anatomyInfo.PgroupLabels),1);repmat(map(1,:),length(anatomyInfo.SgroupLabels),1)];
% cfg=[]; 
% cfg.view='li';
% cfg.ignoreDepthElec='n';
% cfg.opaqueness=0.3;
% cfg.elecSize = 2;
% cfg.elecColors = elecColors;
% cfg.elecColorScale = [0 64];
% cfg.showLabels='n';
% cfg.title= '';
% cfg.elecCoord = [[anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];[anatomyInfo.SgroupAvgCoords,anatomyInfo.SgroupIsLeft]];  
% cfg.elecNames = [anatomyInfo.PgroupLabels,anatomyInfo.SgroupLabels];
% cfg.axis = ax5;
% cfgOut=plotPialSurf('fsaverage',cfg);
% cfg.view='ri';
% cfg.axis = ax6;
% cfgOut=plotPialSurf('fsaverage',cfg);

elecColors = [repmat([0.1098,0.0980,0.8431],length(anatomyInfo.PgroupLabels),1);repmat([0,0,0],length(anatomyInfo.SgroupLabels),1)];
[subGroup, map] = defineSubGroups_stimAnatomy();

% Make mixed phase a darker shade of gray
map(3,:) = [0.3, 0.3,0.3];
for p_ct = 1:length(anatomyInfo.patients)
    pairs{p_ct,1} = anatomyInfo.PgroupLabels{p_ct};
    pairs{p_ct,2} = anatomyInfo.SgroupLabels{p_ct};
    
    for ii_a = 1:3
        if (ismember(str2num(anatomyInfo.PgroupLabels{p_ct}(1:3)),subGroup{ii_a}))
            stimV = ii_a;
            pairs{p_ct,3} = map(stimV,:);
        end
    end
    pairs{p_ct,4} = anatomyInfo.PgroupLabels{p_ct}(5);
end

% We cannot plot p505 in the same way as other subjects because stim and
% probe are on opposite hemispheres - 
pairs{11,1} = pairs{14,1};
pairs{11,2} = pairs{14,2};
pairs{11,3} = pairs{14,3};
pairs{11,4} = pairs{14,4};

cfg=[]; 
cfg.surfType = 'pial';
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.2;
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
cfg.lineWidth = 0.8;
cfg.edgeBlack = 'n';
cfg.axis = ax7;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='ri';
cfg.axis = ax8;
cfgOut=plotPialSurf('fsaverage',cfg);

outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');

a = gcf;
set(f0,'renderer','zbuffer');
res = 400;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
res = 400;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!

%% Different colors per pt

outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures','perPt');

elecColors = [];
for ii = 1:length(anatomyInfo.patients)
    
    if ii == 11; continue; end

    figName = sprintf('STIM_PROBE_locations_%s',anatomyInfo.patients{ii});
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
    
   
    elecColors(1,:) = [0.1098,0.0980,0.8431];
    elecColors(2,:) = [0,0,0];

    cfg=[];
    cfg.surfType = 'pial';
    cfg.ignoreDepthElec='n';
    cfg.opaqueness=0.2;
    cfg.elecSize = 2;
    cfg.elecColors = elecColors;
    cfg.elecColorScale = [0 64];
    cfg.showLabels='n';
    cfg.title= '';
    cfg.elecCoord = [[anatomyInfo.PgroupAvgCoords(ii,:),anatomyInfo.PgroupIsLeft(ii,:)];...
        [anatomyInfo.SgroupAvgCoords(ii,:),anatomyInfo.SgroupIsLeft(ii,:)]];
    cfg.elecNames{1} = [anatomyInfo.PgroupLabels{ii}];
    cfg.elecNames{2} = [anatomyInfo.SgroupLabels{ii}];
    
    for p_ct = 1:10;
        pair{p_ct,1} = anatomyInfo.PgroupLabels{ii};
        pair{p_ct,2} = anatomyInfo.SgroupLabels{ii};
        for ii_a = 1:3
            if (ismember(str2num(anatomyInfo.PgroupLabels{ii}(1:3)),subGroup{ii_a}))
                stimV = ii_a;
                pair{p_ct,3} = map(stimV,:);
            end
        end
        pair{p_ct,4} = anatomyInfo.PgroupLabels{p_ct}(5);
    end

    cfg.pairs = pair;
    cfg.lineWidth = 0.8;
    cfg.edgeBlack = 'n';
    
    if anatomyInfo.PgroupIsLeft(ii,:) % Electrodes on left hem
        cfg.view='li';
        cfg.axis = ax1;
        cfgOut=plotPialSurf('fsaverage',cfg);

        axes(ax2)
        title([anatomyInfo.patients{ii},cfg.elecNames])
        axis off
    else
        cfg.view='ri';
        cfg.axis = ax2;
        cfgOut=plotPialSurf('fsaverage',cfg);

        axes(ax1)
        title([anatomyInfo.patients{ii},cfg.elecNames])
        axis off
    end

    a = gcf;
    set(f0,'renderer','zbuffer');
    res = 400;
    eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
    close all
end
