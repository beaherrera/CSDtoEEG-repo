%% calculate dipolar moment in the column coord system
% Author: Beatriz Herrera

clear
clc

%% paths
% folders name: L3_simData_10ms and L5_simData_10ms
pathData = 'L3_simData_10ms'; % name of the folder where simulated data
% is stored

pathCurrent = []; % path to neurons' transmembrane currents (too large
% to include on GitHub repo) generated during biophysical sim

pathSourceCoord = []; % location of the neurons in the cortical column

savePath = fullfile(pathCurrent,'dipoles');
if ~exist(savePath, 'dir') % checks if the folder already exists
    mkdir(savePath);  % creates a folder named 'file'
end

%% parameters

num_neurons = 1000; % number of simulated neurons
% 2200 L3 PCs and 1000 L5 PCs were simulated

load(fullfile(pathData,'simLFP.mat'), "ze");
h = mean(diff(ze)); % [m] inter-electrodes distance
Ne = length(ze); % number of electrodes
a = ze(1); % position of the first electrode

%% calculate the dipoles

d = 0;

for ii=0:(num_neurons - 1)

    sprintf('neuron %d', ii)

    file = fullfile(pathCurrent, ['NeuronsData_r1_n#' num2str(ii) '.mat']);
    load(file, 'It') % nA

    load(fullfile(pathSourceCoord, ['xyzCoord_n' num2str(ii) '.mat']))

    di = ([xc; yc; zc] + [0; 0; (((h*Ne)/2)+a)])*It;
    d = d + di;

    save(fullfile(savePath, ['dipole_n#' num2str(ii) '.mat']), 'di')
end

%% save dipolar moment

save(fullfile(savePath, 'column_eq_dip.mat'), 'd')
