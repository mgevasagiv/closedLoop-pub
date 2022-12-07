function [elecMatrix, elecLabels, elecRgb] = csvElec2matlab(fsSub)
fsDir=getFsurfSubDir();
mlocFname = sprintf('%s_fsl_based',fsSub);
filename=fullfile(fsDir,fsSub,'elec_recon',[mlocFname '.csv']);
if isempty(dir(filename))
    disp(filename)
    error('csv mLoc file missing')
end

parcOut = csv2Cell(filename);

for ii = 1:length(parcOut)
    if isempty(parcOut{ii,1})
        disp(sprintf('%s - %d electrodes found in CSV',fsSub, ii))
        break
    end
    elecLabels{ii} = sprintf('%s%s_%s_%s',parcOut{ii,6},parcOut{ii,7},parcOut{ii,1},parcOut{ii,2});
    elecMatrix(ii,1:3) = [str2num(parcOut{ii,3}),str2num(parcOut{ii,4}),str2num(parcOut{ii,5})];
    elecRgb(ii,:) = [1,0,0];
end
