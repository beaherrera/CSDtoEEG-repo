%% Fig 2G

clear
clc

%% paths

EEG_folder = []; % path to simulated EEGs for each approach

%% Simulation time
t0 = 0; %[ms] start time
tf = 800; % [ms] stop time
dt = 2^(-7);  %[ms] time increment
Fs =1/(dt*1e-3); % [Hz] sampling frequency
ts = t0:dt:tf;  %[ms] time span
Fs_new = 1e3; % 1kHz new sampling frequency

%% calculate errors

vert =[102, 431, 437, 862, 1835, 3594, 4070, ...
    4117, 4180, 5949, 7441, 14288, 13461, 11787, 11701];

MAG_values = zeros(length(vert),2);
RDM_values = zeros(length(vert),2);

for ii=1:length(vert)

    [MAG_valuesi, RDM_valuesi] = error_fun(vert(ii), EEG_folder, ts, Fs_new);

    MAG_values(ii, :) = MAG_valuesi;

    RDM_values(ii, :) = RDM_valuesi;

end

%% save errors

save('simData_10ms\EEG_errors.mat', 'MAG_values','RDM_values')  

%% generate figure

appName = [repmat({'1'},1,length(RDM_values(:,1))),...
    repmat({'2'},1,length(RDM_values(:,1)))]';
appNameC = categorical(appName);
allRDM = [RDM_values(:,1)', RDM_values(:,2)'];
allMAG = [MAG_values(:,1)', MAG_values(:,2)'];

figure;
subplot(2,1,1); 
b = boxchart(allRDM, 'GroupByColor',appNameC);
b(1, 1).BoxFaceColor = 'k';
b(2, 1).BoxFaceColor = 'r';
ylim([-0.01 0.5])
ylabel('RDM')
xticklabels({})
legend({'In','CSD'},'Location','best','Box','off')
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold')
hold on; 
ax = axes('Position',[.35 .74 .15 .15]);
box on
b = boxchart(allRDM, 'GroupByColor',appNameC);
b(1, 1).BoxFaceColor = 'k';
b(2, 1).BoxFaceColor = 'r';
ylim([0 0.02])
xticklabels({})
set(ax,'linewidth',1.5,'fontweight','bold')

subplot(2,1,2); 
b = boxchart(allMAG, 'GroupByColor',appNameC);
hold on; yline(1,'--k','LineWidth',1)
b(1, 1).BoxFaceColor = 'k';
b(2, 1).BoxFaceColor = 'r';
ylabel('MAG')
xticklabels({})
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold')

%% subfunctions


function [MAG_values, RDM_values] = error_fun(vert_ind, EEG_folder, ts, Fs_new)
% calculate errors
% 
% Inputs:
%   vert_ind: vertex number
%   EEG_folder: path to sim EEG data
%   ts: simulation time in sec
%   Fs_new: new sampling rate to downsample detailed approach data  
%
% Outputs:
%   MAG_values: magnitude error values
%   RDM_values: relative difference error values

sprintf('vertex %d',vert_ind)

EEGIn = 0;

for n = 0:999

    load(fullfile(EEG_folder, 'EEG_In_10ms', ...
        ['EEGIn_Vert' num2str(vert_ind) 'C' num2str(n) '.mat']), ...
        'EEG_In');
    EEGIn = EEGIn + EEG_In;

end
% -- downsampling the LFPs to 1kHz
if ~brainstorm('status')
    brainstorm nogui
end
[EEGIn_s, ~] = process_resample('Compute', EEGIn, ts*1e-3, Fs_new);
EEG_In_s = EEGIn_s(:,590:(1e3/Fs_new):690);

save(fullfile(EEG_folder, 'EEG_In_10ms', ...
        ['EEGIn_Vert' num2str(vert_ind) '.mat']),'EEGIn', 'EEG_In_s')

load(fullfile(EEG_folder, 'EEG_Dip_10ms', ...
    ['EEGDip_Vert' num2str(vert_ind) '.mat']), ...
    'EEG_Dip');
load(fullfile(EEG_folder, 'EEG_DipCSD_10ms', ...
    ['EEGDip_Vert' num2str(vert_ind) '.mat']), ...
    'EEG_CSD');

EEG_In_f = eegfilt(EEG_In_s, Fs_new, 0, 100);
EEG_Dip_f = eegfilt(EEG_Dip, Fs_new, 0, 100);
EEG_CSD_f = eegfilt(EEG_CSD, Fs_new, 0, 100);

MAG_values = [mean(MAG(EEG_In_f, EEG_Dip_f),'omitnan'), ...
    mean(MAG(EEG_In_f, EEG_CSD_f),'omitnan')];

RDM_values = [mean(RDM(EEG_In_f, EEG_Dip_f),'omitnan'),...
    mean(RDM(EEG_In_f, EEG_CSD_f),'omitnan')];

end