
clear
clc

%% load EEGs
pathData = 'L3andL5_EEG_1'; % name of the folder where simulated data is stored

% -- L3 PCs --
load(fullfile('L3_simData_10ms','filtered_simEEG.mat'), 'EEG_In_f', ...
    "EEG_Dip_f","EEG_CSD_f")
EEG_In_f_L3 = EEG_In_f; clearvars EEG_In_f
EEG_Dip_f_L3 = EEG_Dip_f; clearvars EEG_Dip_f
EEG_CSD_f_L3 = EEG_CSD_f; clearvars EEG_CSD_f

load(fullfile('L3_simData_10ms','simLFP.mat'), 'VeD_1');
VeD_1_L3 = VeD_1; clearvars VeD_1

% -- L5 PCs --
load(fullfile('L5_simData_10ms','filtered_simEEG.mat'), 'EEG_In_f', ...
    "EEG_Dip_f","EEG_CSD_f")
EEG_In_f_L5 = EEG_In_f; clearvars EEG_In_f
EEG_Dip_f_L5 = EEG_Dip_f; clearvars EEG_Dip_f
EEG_CSD_f_L5 = EEG_CSD_f; clearvars EEG_CSD_f

load(fullfile('L5_simData_10ms', pathData,'simLFP.mat'), 'VeD_1', 'ze', ...
    "ts_out", "Fs_new");
VeD_1_L5 = VeD_1; clearvars VeD_1

%% total LFP

VeD_1 = VeD_1_L3 + VeD_1_L5;

%% plot the LFP
font = 10;
Ne = length(ze); % number of electrodes in the shank

figure;
LFPmax = 8.5*max(abs(VeD_1.*1e3),[],'all'); % scale for the LFP plot
plot(ts_out.*1e3, (ze.*1e3.*ones(length(ts_out), length(ze)))' + 0.5,'-', ...
    'Color', [0.5 0.5 0.5], 'linewidth', 1)
hold on
for ii = 1:Ne
    plot(ts_out.*1e3,VeD_1(Ne-ii+1,:).*1e3./LFPmax + ze(ii).*1e3 + 0.5, ...
        'color', 'r', 'clipping','on', 'linewidth', 1.5)
end
hold on
plot((599-10+80)-[20 40], [-0.03; -0.03] + 0.5, '-k', ...
    'LineWidth', 1.5);  % time scale
hold on
plot((599-10+80)-[20 20], [-3e-2 ; 2e-2]./LFPmax+ ze(1).*1e3 + 0.5, '-k', ...
    'LineWidth', 1.5); % voltage scale
hold off
text((599-10+80)-mean([20 40]), -0.08 + 0.5, '20 ms', ...
    'HorizontalAlignment','center', 'fontsize', font, ...
    'FontWeight','bold') % time scale label
text((599-10+80)-18, mean([-2.5e-2 ; 0e-2]./LFPmax+ ze(1).*1e3) + 0.5, ...
    '0.05 mV', 'HorizontalAlignment','left', 'fontsize', font ...
    ,'FontWeight','bold') % voltage scale label
ylabel('Cortical Depth (mm)','fontsize',font,'fontweight','bold')
yticks((0.05:0.2:1.65) + 0.5)
values_y = sort(0.05:0.2:1.65,'descend');
tickslabLFP = {};
for ii=values_y
    tickslabLFP = cat(2, tickslabLFP, num2str(ii,'%4.2f'));
end
yticklabels(tickslabLFP)
ax = gca;
ax.FontSize = font;
ylim([-0.05 1.7+0.05] + 0.5)
xlim([(599-10) (599-10+80)])
box 'off'
set(ax,'XTickLabel',[],'XTick',[],'XColor',[1 1 1],'FontWeight','bold',...
    'FontSize',12,'LineWidth',2)

%% calculate the CSD
%% ----- Parameters
el_pos = ze; % depth coordinates in m
a = el_pos(1); % position of the first electrode
h = 0.1e-3; % [m] inter-electrodes distance
% convert mm -> m
diam = 3e-3; % [m] cylinder diameter,
% convert mm -> m
rc = diam/2; % [m] radius of the cortical column
cond = 0.33; %[S/m] gray matter conductance
cond_top = cond;
gauss_sigma = 0.1e-3;   %[m] Gaussian
% filter std, convert mm -> m
filter_range = 5*gauss_sigma; % numeric filter must be finite in
% extent

%% ---- solve Pettersen model
Fcs = F_cubic_spline(el_pos,diam,cond,cond_top); % cubic spline
% function

[zs_1,CSD_cs_1] = make_cubic_splines(el_pos,VeD_1,Fcs); % compute CSD
if ~isempty(gauss_sigma) && gauss_sigma~=0 % filter iCSD
    [zs_1,CSD_cs_1] = gaussian_filtering(zs_1,CSD_cs_1,gauss_sigma,...
        filter_range);
end
iCSD_1 = CSD_cs_1; % [nA/mm3] current source density
zs_1 = zs_1*1e3; % depth in m, convert to mm

%%
save(fullfile(pathData,'simCSD.mat'), "zs_1","iCSD_1")

%% figure: supra-CF stim

tspan = (590:(1e3/Fs_new):699) - 590; % generate time span of the data
font = 14; % font size
figure;
imagesc(tspan, zs_1, iCSD_1(:,590:(1e3/Fs_new):699));
yticks(0.05:0.2:1.65)
c = colorbar;
colormap(flip(jet));
c.Label.String = 'nA/mm^{3}';
max_CSD = max(abs(iCSD_1), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
c.Label.FontSize = font;
c.FontSize = font;
c.FontWeight = 'bold';
ax = gca;
ax.FontSize = font;
xlabel('Time (ms)', 'FontSize',font);
ylabel('Cortical Depth (mm)','FontSize',font);
xlim([0 80])
set(ax,'fontweight','bold','FontSize',12,'LineWidth',2)

%% CSD dipole

d_CSD_1 = (-(zs_1.*1e-3 - (((h*Ne)/2)+a))'.*mean(diff(zs_1.*1e-3)))'*...
    iCSD_1(:,590:(1e3/Fs_new):699)*(pi*(rc^2)).*1e9; % A*m -> nA*m
% baseline correction
d_CSD_1 = d_CSD_1 - mean(d_CSD_1(1:9));

%% comparison of the dipoles in the column's coord system

load(fullfile([],'column_eq_dip.mat'),'d');
dL5 = d; clearvars d
load(fullfile([],'column_eq_dip.mat'),'d');
dL3 = d; clearvars d

d_f = eegfilt(dL3 + dL5, Fs, 0, 100);
[ds_f, ~] = process_resample('Compute', d_f, ts*1e-3, Fs_new);

d_CSD_xyz = [0 0 1]'*d_CSD_1;

%% --- errors
MAG(ds_f(:,590:(1e3/Fs_new):699) - ...
    mean(ds_f(:,590:(1e3/Fs_new):600), 2), d_CSD_xyz)
RDM(ds_f(:,590:(1e3/Fs_new):699) - ...
    mean(ds_f(:,590:(1e3/Fs_new):600), 2),d_CSD_xyz)

%%
figure;
subplot(3,1,1)
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
    ds_f(1,590:(1e3/Fs_new):699) - ...
    mean(ds_f(1,590:(1e3/Fs_new):600)), '-k','LineWidth',1)
hold on;
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
    d_CSD_xyz(1,:), '-r','LineWidth',1)
ylabel({'x'})
ylim([min(d_CSD_xyz(3,:))-0.01 max(d_CSD_xyz(3,:))+0.01])
title('Current Dipole Moment (nA*m)')
box off
xlim([0 80])
set(gca,'fontweight','bold','FontSize',9,'LineWidth',1.5, ...
    'XTickLabel',[],'XTick',[],'XColor',[1 1 1])

subplot(3,1,2)
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
    ds_f(2,590:(1e3/Fs_new):699) - ...
    mean(ds_f(2,590:(1e3/Fs_new):600)), '-k','LineWidth',1)
hold on;
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
    d_CSD_xyz(2,:), '-r','LineWidth',1)
ylim([min(d_CSD_xyz(3,:))-0.01 max(d_CSD_xyz(3,:))+0.01])
ylabel({'y'})
% xlabel('Time (ms)')
legend({'STC','CSD'},'Box','off')
box off
xlim([0 80])
set(gca,'fontweight','bold','FontSize',9,'LineWidth',1.5, ...
    'XTickLabel',[],'XTick',[],'XColor',[1 1 1])

subplot(3,1,3)
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
    ds_f(3,590:(1e3/Fs_new):699) - ...
    mean(ds_f(3,590:(1e3/Fs_new):600)), '-k','LineWidth',1)
hold on;
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
    d_CSD_xyz(3,:), '-r','LineWidth',1)
ylim([min(d_CSD_xyz(3,:))-0.01 max(d_CSD_xyz(3,:))+0.01])
ylabel({'z'})
xlabel('Time (ms)')
box off
xlim([0 80])
set(gca,'fontweight','bold','FontSize',9,'LineWidth',1.5)

%% V4 EEG

EEG_In_f = EEG_In_f_L3 + EEG_In_f_L5;
EEG_Dip_f = EEG_Dip_f_L3 + EEG_Dip_f_L5;
EEG_CSD_f = EEG_CSD_f_L3 + EEG_CSD_f_L5;

%% save
save(fullfile(pathData,'filtered_simEEG.mat'), 'EEG_In_f', ...
    "EEG_Dip_f","EEG_CSD_f")

%% EEG plot

figure;
subplot(1,3,1)
plot(ts_out(590:(1e3/Fs_new):699).*1e3-590, (EEG_In_f - ...
    mean(EEG_In_f(:,1:10),2))')
hold on;
xline(20,'--')
hold on;
xline(55,'--')
ylabel('EEG (\muV)')
box off
xlim([0 80])
ylim([-0.3 0.3])
set(gca,'fontweight','bold','FontSize',10,'LineWidth',2)

subplot(1,3,2)
plot(ts_out(590:(1e3/Fs_new):699).*1e3-590, EEG_Dip_f')
hold on;
xline(20,'--')
hold on;
xline(55,'--')
xlabel('Time (ms)')
box off
xlim([0 80])
ylim([-0.3 0.3])
set(gca,'fontweight','bold','FontSize',10,'LineWidth',2)

subplot(1,3,3)
plot(ts_out(590:(1e3/Fs_new):699).*1e3-590, EEG_CSD_f')
hold on;
xline(20,'--')
hold on;
xline(55,'--')
box off
ylim([-0.3 0.3])
xlim([0 80])
set(gca,'fontweight','bold','FontSize',10,'LineWidth',2)

%% calculate the errors for V4 columns

mean(MAG(EEG_In_f - mean(EEG_In_f(:,1:10),2),EEG_Dip_f))
mean(RDM(EEG_In_f - mean(EEG_In_f(:,1:10),2),EEG_Dip_f))

mean(MAG(EEG_In_f - mean(EEG_In_f(:,1:10),2),EEG_CSD_f))
mean(RDM(EEG_In_f - mean(EEG_In_f(:,1:10),2),EEG_CSD_f))

mean(MAG(EEG_Dip_f,EEG_CSD_f))
mean(RDM(EEG_Dip_f,EEG_CSD_f))

