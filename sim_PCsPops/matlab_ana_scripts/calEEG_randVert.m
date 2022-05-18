%% Generate data for Fig 2G

clear
clc

%% Simulation time
t0 = 0; %[ms] start time
tf = 800; % [ms] stop time
dt = 2^(-7);  %[ms] time increment
Fs =1/(dt*1e-3); % [Hz] sampling frequency
ts = t0:dt:tf;  %[ms] time span

%% Folder Paths

pathData = 'simData_10ms'; % name of the folder where simulated data is stored

pathLeadField = []; % path to lead fields

path2CellsData = []; % path to neurons' transmembrane currents (too large
% to include on GitHub repo) generated during biophysical sim

%% load
load(fullfile(path2CellsData, 'dipoles','column_eq_dip.mat'),'d');
% -- downsampling the LFPs to 1kHz
if ~brainstorm('status')
    brainstorm nogui
end
Fs_new = 1e3; % 1kHz new sampling frequency
[ds_In, ts_out] = process_resample('Compute', d, ts*1e-3, Fs_new);

load(fullfile(pathData,'EEG_DipCSD.mat'),'d_CSD')

%% create saving directories
fileSave = fullfile(pathLeadField, 'EEG_In_10ms');
if ~exist(fileSave, 'dir') % checks if the folder already exists
    mkdir(fileSave);  % creates a folder named 'file'
end

fileSaveDip = fullfile(pathLeadField, 'EEG_Dip_10ms');
if ~exist(fileSaveDip, 'dir') % checks if the folder already exists
    mkdir(fileSaveDip);  % creates a folder named 'file'
end

fileSaveDipCSD = fullfile(pathLeadField, 'EEG_DipCSD_10ms');
if ~exist(fileSaveDipCSD, 'dir') % checks if the folder already exists
    mkdir(fileSaveDipCSD);  % creates a folder named 'file'
end

%% calculate EEG for a column center at vert's coordinates

vert = [102, 431, 437, 862, 1835, 3594, 4070, ...
    4117, 4180, 5949, 7441, 14288, 13461, 11787, 11701]; 

% p = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(p)
%     parpool(4);
% end

for jj = 1:length(vert)

    % detailed approach
    cal_EEGIn(vert(jj), pathLeadFieldIN, path2CellsData, fileSave)
    
    % dipolar approaches
    cal_EEGdip(vert(jj), pathLeadField, ds_In, d_CSD, fileSaveDip, ...
        fileSaveDipCSD)

end


%% sub-function

function cal_EEGIn(vert, pathLeadFieldIN,pathIt, fileSave)

% sprintf('vertex %d',vert)

EEGIn = 0;

for ii = 0:999

    % load lead field
    load(fullfile(pathLeadFieldIN,'leadField', ...
        ['leadFieldMonBEMVert' num2str(vert) ...
        'C' num2str(ii) '.mat']),'Keoo');

    % convert to average reference
    Ne = length(Keoo(:,1));
    Hn = eye(Ne) - (ones(Ne,1)*ones(Ne,1)')/Ne;

    Ke = Hn*Keoo;

    % load transmembrane currents
    file = fullfile(pathIt, ['NeuronsData_r1_n#' num2str(ii) '.mat']);
    load(file, 'It') % nA

    % calculate EEG
    EEG_In = Ke*It; % nV
    EEGIn = EEGIn + EEG_In; % nV

    % save EEG
    save(fullfile(fileSave, ...
        ['EEGIn_Vert' num2str(vert) 'C' num2str(ii) '.mat']), ...
        'EEG_In');

end

end

function cal_EEGdip(vertInd, path2leadField, ds, d_CSD, fileSaveDip, ...
    fileSaveDipCSD)

% load lead field
filename = ['leadFieldDipBEMVert' num2str(vertInd) '.mat'];
load(fullfile(path2leadField, 'leadField', filename),'Keoo');

% calculate the average reference
Ne = length(Keoo(:,1));
Hn = eye(Ne) - (ones(Ne,1)*ones(Ne,1)')/Ne;

% change lead fields reference to the average reference
KeDip = Hn*Keoo;

% load normal
load(fullfile(path2leadField, 'sources_coordinates', ...
    ['coordinatesVert' num2str(vertInd) ...
    '_n#' num2str(0) '.mat']), ...
    'xyzCortNorm');

% rotate column
[Rx,Ry] = calRotMatriz(xyzCortNorm);
d_In = Rx'*Ry'*ds;

% calculate the EEG
EEG_CSD = KeDip*(xyzCortNorm'*d_CSD); % nV
EEG_Dip = KeDip*(d_In(:,590:690)); % nV

% save EEGs
save(fullfile(fileSaveDip, ...
    ['EEGDip_Vert' num2str(vertInd) '.mat']), ...
    'EEG_Dip');
save(fullfile(fileSaveDipCSD, ...
    ['EEGDip_Vert' num2str(vertInd) '.mat']), ...
    'EEG_CSD');

end
