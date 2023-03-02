% n = 18 Pts recorded on 2017-2021 in UCLA
% Defining stim vs mixed stimulation cohorts for plots
function [subGroup, cmap] = defineSubGroups_stimAnatomy(type)

if nargin < 1
    type = 1;
end

G_targeted = 1;
G_targeted_diff_anatomy = 2;
G_mixed = 3;
sham_night_index = 4;

cmap(G_targeted,:) = [0.8431,0.0980,0.1098];
cmap(G_targeted_diff_anatomy,:) = [0.6510    0.3373    0.1569];
cmap(G_mixed,:) = [0.6,0.6,0.6];
cmap(sham_night_index,:) = [0.3020    0.6863    0.2902]; 

if type == 1
    subGroup{1} = [486, 487, 488, 498, 515, 520, 545, 541, 544];% frontal,well_targeted_stim_pts

    subGroup{2} = [490,489,496, 538]; % non-frontal,well_targeted_stim_pts
    
    subGroup{3} = [485, 497, 499, 505,510]; % frontal_stim_mixed_phase_pts
    % 485 - two probes were used instead of one
    % 497 - totally random (BR signals unresponsive)
    % 499 - phase targetting is all 'OFF' 
    % 505, 510 - were run with alternating 'ON' vs 'OFF' targetting
end

end