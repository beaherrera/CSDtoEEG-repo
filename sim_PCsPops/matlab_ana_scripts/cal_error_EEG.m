%% Fig 2J-K: calculation of RDM and MAG measures.
% Author: Beatriz Herrera

clear
clc

%% paths to simulated EEGs for each neuronal population
EEG_L5_folder = [''];
EEG_L3_folder = [''];

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

MAG_values = zeros(length(vert),3);
RDM_values = zeros(length(vert),3);

for ii=1:length(vert)

    [MAG_valuesi, RDM_valuesi] = error_fun(vert(ii), EEG_L5_folder, ...
        EEG_L3_folder, Fs_new);

    MAG_values(ii, :) = MAG_valuesi;

    RDM_values(ii, :) = RDM_valuesi;

end

%% save errors

save(fullfile(pathData, 'EEG_errors.mat'), 'MAG_values','RDM_values')

%% error Vs columns' distance to the scalp
%% get distance from the columns to the scalp surface

% load scalp surface 
load(fullfile('NMTv2_atlas_surfaces', ...
    'tess_head_bem_1922V_02_fig_02_fig.mat'), 'Faces', 'Vertices')

% calculate the minimal distance between the center of the column and the
% scalp surface
TR.faces = Faces;
TR.vertices = Vertices;
points = coord_vert;
[distances, surface_points] = point2trimesh(TR, 'QueryPoints', points);

% plot the columns location and the shortest vector to the scalp surface
figure; patch(TR,'FaceAlpha',.1, 'EdgeColor', 'none'); 
xlabel('x'); ylabel('y'); zlabel('z'); axis equal; hold on
plot3M = @(XYZ,varargin) plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),varargin{:});
plot3M(points,'*r')
plot3M(surface_points,'*k')
plot3M(reshape([shiftdim(points,-1);shiftdim(surface_points,-1); ...
    shiftdim(points,-1)*NaN],[],3),'k')

%% plot error vs distance from scalp / depth

[fitobject_MAG_In1, gof_MAG_In1] = fit(distances.*100, ...
    MAG_values(:,1), 'poly1');
[fitobject_MAG_CSD1, gof_MAG_CSD1] = fit(distances.*100, ...
    MAG_values(:,2), 'poly1');

[fitobject_RDM_In1, gof_RDM_In1] = fit(distances.*100, ...
    RDM_values(:,1), 'poly1');
[fitobject_RDM_CSD1, gof_RDM_CSD1] = fit(distances.*100, ...
    RDM_values(:,2), 'poly1');

figure;
subplot(1,2,1)
hold on
plot(distances.*100, MAG_values(:,1), '.k', 'MarkerSize', 8)
plot(distances.*100, ...
    feval(fitobject_MAG_In1,distances.*100), '-k')
plot(distances.*100, MAG_values(:,2), '.r', 'MarkerSize', 8)
plot(distances.*100, ...
    feval(fitobject_MAG_CSD1,distances.*100), '-r')
legend({'In', ['slope: ' num2str(round(fitobject_MAG_In1.p1,4)) '' ...
    ' R^2: ' num2str(round(gof_MAG_In1.rsquare,4))], ...
    'CSD', ['slope: ' num2str(round(fitobject_MAG_CSD1.p1,4)) '' ...
    ' R^2: ' num2str(round(gof_MAG_CSD1.rsquare,4))]}, ...
    'Location','east')
xlabel({'Cortcal column depth (cm)'})
ylabel('MAG')
box off
set(gca,'linewidth',1.5,'fontsize',9,'fontweight','bold')

subplot(1,2,2)
hold on
plot(distances.*100, RDM_values(:,1)', '.k', 'MarkerSize', 8)
plot(distances.*100, ...
    feval(fitobject_RDM_In1,distances.*100), '-k')
plot(distances.*100, RDM_values(:,2)', '.r', 'MarkerSize', 8)
plot(distances.*100, ...
    feval(fitobject_RDM_CSD1,distances.*100), '-r')
legend({'In', ['slope: ' num2str(round(fitobject_RDM_In1.p1,4)) '' ...
    ' R^2: ' num2str(round(gof_RDM_In1.rsquare,4))], ...
    'CSD', ['slope: ' num2str(round(fitobject_RDM_CSD1.p1,4)) '' ...
    ' R^2: ' num2str(round(gof_RDM_CSD1.rsquare,4))]}, ...
    'Location','east')
xlabel({'Cortcal column depth (cm)'})
ylabel('RDM')
box off
set(gca,'linewidth',1.5,'fontsize',9,'fontweight','bold')

%% subfunctions

function [MAG_values, RDM_values] = error_fun(vert_ind, EEG_L5_folder, EEG_L3_folder, Fs_new)

sprintf('vertex %d',vert_ind)

% load EEG ground-truth
load(fullfile(EEG_L3_folder, 'EEG_In_10ms', ...
    ['EEGIn_Vert' num2str(vert_ind) '.mat']), 'EEG_In_s')
EEG_L3_In_s = EEG_In_s;
load(fullfile(EEG_L5_folder, 'EEG_In_10ms', ...
    ['EEGIn_Vert' num2str(vert_ind) '.mat']), 'EEG_In_s')
EEG_L5_In_s = EEG_In_s;

% load EEG STC dipole
load(fullfile(EEG_L3_folder, 'EEG_Dip_10ms', ...
    ['EEGDip_Vert' num2str(vert_ind) '.mat']), ...
    'EEG_Dip');
EEG_L3_Dip = EEG_Dip(:,1:101);

load(fullfile(EEG_L5_folder, 'EEG_Dip_10ms', ...
    ['EEGDip_Vert' num2str(vert_ind) '.mat']), ...
    'EEG_Dip');
EEG_L5_Dip = EEG_Dip;

% load EEG CSD dipole
load(fullfile(EEG_L3_folder, 'EEG_DipCSD_10ms', ...
    ['EEGDip_Vert' num2str(vert_ind) '.mat']), ...
    'EEG_CSD');
EEG_L3_CSD = EEG_CSD(:,1:101);

load(fullfile(EEG_L5_folder, 'EEG_DipCSD_10ms', ...
    ['EEGDip_Vert' num2str(vert_ind) '.mat']), ...
    'EEG_CSD');
EEG_L5_CSD = EEG_CSD;

% filter EEG at 100Hz
EEG_In_f = eegfilt(EEG_L3_In_s + EEG_L5_In_s, Fs_new, 0, 100);
EEG_Dip_f = eegfilt(EEG_L3_Dip + EEG_L5_Dip, Fs_new, 0, 100);
EEG_CSD_f = eegfilt(EEG_L3_CSD + EEG_L5_CSD, Fs_new, 0, 100);

% calculate error measures
MAG_values = [mean(MAG(EEG_In_f, EEG_Dip_f),'omitnan'), ...
    mean(MAG(EEG_In_f, EEG_CSD_f),'omitnan'),...
    mean(MAG(EEG_Dip_f, EEG_CSD_f),'omitnan')];

RDM_values = [mean(RDM(EEG_In_f, EEG_Dip_f),'omitnan'),...
    mean(RDM(EEG_In_f, EEG_CSD_f),'omitnan'),...
    mean(RDM(EEG_Dip_f, EEG_CSD_f),'omitnan')];

end