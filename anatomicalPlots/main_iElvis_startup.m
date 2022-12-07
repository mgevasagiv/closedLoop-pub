global globalFsDir;
globalFsDir='E:\Data_p\FreeSurferWinMirror';
FirstRun = false; 

% Check iElvis installation works - 
% Depends on download of PT001 and fsaverage from
% their website.
if (FirstRun)
    check_iELVisInstall
end

% PTs with freesurfer-coreg-BioImage ran locally - 
pt_ii = 1; % ran by Maya - mask done, checked 
ptNum{pt_ii} = '485'; 

pt_ii = pt_ii+1; % ran by Maya - mask done, checked 
ptNum{pt_ii} = '486'; 

pt_ii = pt_ii+1; % ran by Maya - mask done, checked 
ptNum{pt_ii} = '487'; 

pt_ii = pt_ii+1; % ran by Maya - mask done, checked  
ptNum{pt_ii} = '488'; 

pt_ii = pt_ii+1;  
ptNum{pt_ii} = '489'; % Data from Hanna - RSTG has only 3 contacts - why?

pt_ii = pt_ii+1;  
ptNum{pt_ii} = '490'; % 

pt_ii = pt_ii+1; % mask done, checked  
ptNum{pt_ii} = '496'; 

pt_ii = pt_ii+1;
ptNum{pt_ii} = '497'; % ran by Maya - mask done, checked 

pt_ii = pt_ii+1;
ptNum{pt_ii} = '498'; % ran by Hanna -  mask done, checked 

pt_ii = pt_ii+1;
ptNum{pt_ii} = '499'; % mask done, checked

pt_ii = pt_ii+1;
ptNum{pt_ii} = '505'; % Data from Hanna - RSTG has only 3 contacts - why?

pt_ii = pt_ii+1; % ran by Maya (VM) - mask done, checked
ptNum{pt_ii} = '510'; 

pt_ii = pt_ii+1; % ran by hanna - mask done, checked
ptNum{pt_ii} = '515'; 

pt_ii = pt_ii+1; % ran by hanna - mask done, checked
ptNum{pt_ii} = '520'; 
   

pt_ii = pt_ii+1; % ran by Maya (VM)
ptNum{pt_ii} = '538'; 
   

pt_ii = pt_ii+1; %  ran by Maya (VM)
ptNum{pt_ii} = '545'; 
   
pt_ii = pt_ii+1; %  ran by Maya (VM)
ptNum{pt_ii} = '544'; 


pt_ii = pt_ii+1; %  ran by Maya (VM)
ptNum{pt_ii} = '541'; 

%%
whatTorun.generateLocTxtFile = false;

% mgrid --> LocTxtFile 
if whatTorun.generateLocTxtFile
    makeIniLocTxtFile(ptNum{pt_ii});
end

% Corrects intracranial electrode locations for brain shift
% Isn't doing anything for depth electrodes but is needed to create the
% files for next parts.
yangWangElecPjct(ptNum{pt_ii});
% Another option - 
% dykstraElecPjct(ptNum{pt_ii});


% Creating images to check if contacts are in the right location -
% ...ptNum\elec_recon\PICS\xxx_LD_area_1Slices
cfg=[]; cfg.printFigs=1;
cfg.markerSize=40;
plotAllDepthsOnSlices(ptNum{pt_ii},'mgrid',cfg)

% Using mLoc CSV file instead
% plotAllDepthsOnSlices_LocFile(ptNum{pt_ii},'mgrid',cfg)

% Reading the txt file with locations into matlab
parcOut=elec2Parc(ptNum{pt_ii},'D');

% Visualization option 1
% Showing 3-D overlay of all electrodes on different directions of the
% brain
cfg=[]; 
cfg.view='omni';
cfg.surfType='pial';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.showLabels='n';
cfg.title=ptNum{pt_ii};
cfgOut=plotPialSurf(ptNum{pt_ii},cfg);
aa = gcf;
res = 200;
outputFigureFolder = fullfile(globalFsDir,ptNum{pt_ii},'elec_recon','PICS');
figName = sprintf('%s_3D_all_electrodes',ptNum{pt_ii});
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(aa.Number),sprintf(' -djpeg  -r%d',res)]); % adding r600 slows down this process significantly!


% Visualization option 2
cfg=[];
cfg.view='omni';
cfg.overlayParcellation='DK';
cfg.showLabels='n';
cfg.surfType='inflated';
cfg.title=sprintf('%s: DK Atlas',ptNum{pt_ii}); 
cfg.elecSize=5;
cfg.ignoreDepthElec='n';
for ii = 1:8
    cfg.elecNames{ii} = parcOut{ii,1};
end
cfg.opaqueness = 0.2;
cfgOut=plotPialSurf(ptNum{pt_ii},cfg);

cfg=[];
cfg.view='l';
cfg.surfType='PIAL';
cfg.ignoreDepthElec='n';
cfg.showLabels='y';
cfgOut=plotPialSurf(ptNum{pt_ii},cfg);


cfg=[];
cfg.view='li';
cfg.surfType='inflated';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.5;
cfg.showLabels='n';
% cfg.elecNames={'LOF1','LOF2','LOF3','LOF4','LOF5','LOF6','LOF7','LOF8','Da1','Da2','Da3','Da4','Da5','Da6','Da7','Da8'};
cfg.title=ptNum{pt_ii};
cfgOut=plotPialSurf(ptNum{pt_ii},cfg);
 
cfg=[];
cfg.view='omni';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.5;
cfg.surfType='inflated';
% cfg.elecNames={'LOF1','LOF2','LOF3','LOF4','LOF5','LOF6','LOF7','LOF8','Da1','Da2','Da3','Da4','Da5','Da6','Da7','Da8'};
for i=1:size(parcOut,1)
elecNames{i}=parcOut{i,1};
end
cfg.title=ptNum{pt_ii};
cfg.showLabels='n';
cfgOut=plotPialSurf(ptNum{pt_ii},cfg)
% % % %  
% % % % 
elecNames=cell(15,1);
for i=1:15 %size(parcOut,1)
elecNames{i}=parcOut{i,1};
end
cfg=[];
cfg.view='omni';
% % cfg.surfType='inflated';

cfg.ignoreDepthElec='n';
cfg.figId=1;
cfg.opaqueness=0.4;
cfg.surfType='inflated';
cfg.elecShape='disk';
cfg.elecColors=linspace(0,1,15)';
cfg.elecColorScale='minmax';
cfg.showLabels='n';
% cfg.elecUnits='r';
cfg.pullOut=1;
cfg.elecNames=elecNames;
cfg.elecSize=6;
cfg.elecCoord='POSTIMPLANT';
cfg.title='PT001: Stimulus Correlations';
cfgOut=plotPialSurf('505',cfg);
 
% Overlaying all electrodes on a avg brain
groupAvgCoords=[];
groupLabels=[];
groupIsLeft=[];
cfg=[];
subs=ptNum;
cfg.plotEm=0;
cfg.isSubdural=0;
for a=1:length(subs),
    fprintf('Working on Participant %s\n',subs{a});
    [avgCoords, elecNames, isLeft]=sub2AvgBrain(subs{a},cfg);
    patient_numoflelectrodes(a)=length(elecNames);
    groupAvgCoords=[groupAvgCoords; avgCoords];
    groupLabels=[groupLabels; elecNames];
    groupIsLeft=[groupIsLeft; isLeft];
end

cfg=[];
cfg.view='l';
% cfg.view='omni';
cfg.elecCoord=[groupAvgCoords groupIsLeft];
cfg.elecNames=groupLabels;
cfg.showLabels='n';
cfg.surfType='inflated';
cfg.opaqueness=0.4;
cfg.elecSize=2;
cfg.title=sprintf('%d pts on Avg. Brain',length(subs));
cfgOut=plotPialSurf('fsaverage',cfg);


% Pending anatomical figures - 
% Generate an image of a stimulating electrode (red), a PROBE electrode (blue),
% recording electrodes (black) - 

% All probe-stimulating couples - 

% Brain wide gain of spindle rate

% Brain wide gain of SW rate

% Generate an image of all MTL-cortical couples

