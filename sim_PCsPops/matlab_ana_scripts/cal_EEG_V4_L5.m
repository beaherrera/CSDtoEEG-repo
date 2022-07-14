%% Fig 2 D-F

clear
clc

%% Folder Paths

pathData = 'L5_simData_10ms'; % name of the folder where simulated data is
% stored (see folder for all simulation results found on the paper)

pathLeadFieldIN = []; % path to point-sources lead field matrices (data
% too large to include on GitHub repo)

path2CellsData_1 = []; % path to supra-thershold
% stim sim data (too large to include on GitHub repo)
path2CellsData_2 = []; % path to threshold
% stim sim data (too large to include on GitHub repo)

path2leadField_V4 = 'leadfields_NMTv2_atlas'; % path to
% dipolar lead field matrix

%% load LFP data
% supra-CF stimulus
load(fullfile(pathData,'SimData_r1_PS1000.mat'), 'LFP','ze')
LFP_1 = LFP;
% threshold CF stimulus
load(fullfile(pathData,'SimData_r1_PS1000_c2.mat'), 'LFP')
LFP_2 = LFP;

%% Simulation time
t0 = 0; %[ms] start time
tf = 800; % [ms] stop time
dt = 2^(-7);  %[ms] time increment
Fs =1/(dt*1e-3); % [Hz] sampling frequency
ts = t0:dt:tf;  %[ms] time span

%% processing the extracellular potentials

% - low pass filter at 100 Hz
Ve_1 = eegfilt(LFP_1, Fs, 0, 100); % [mv] local field potentials
Ve_2 = eegfilt(LFP_2, Fs, 0, 100);

% -- downsampling the LFPs to 1kHz
if ~brainstorm('status')
    brainstorm nogui
end

Fs_new = 1e3; % 1kHz new sampling frequency
[VeD_1, ts_out] = process_resample('Compute', Ve_1, ts*1e-3, Fs_new);
VeD_1 = VeD_1.*1e-3; % convert from mV to V

[VeD_2, ~] = process_resample('Compute', Ve_2, ts*1e-3, Fs_new);
VeD_2 = VeD_2.*1e-3; % convert from mV to V
ze = double(-ze).*1e-6; % change z axis direction and convert from um to m
Ne = length(ze); % number of electrodes in the shank

%% save processed data
save(fullfile(pathData,'simLFP.mat'), 'VeD_1', 'VeD_2', ...
    "ze", "ts_out", "Fs_new");

%% visualization
%% figure: supra CF LFP
font = 10;
figure;
LFPmax = 7*max(abs(VeD_1.*1e3),[],'all'); % scale for the LFP plot
plot(ts_out.*1e3, (ze.*1e3.*ones(length(ts_out), length(ze)))' + 0.5,'-', ...
    'Color', [0.5 0.5 0.5], 'linewidth', 1)
hold on
for ii = 1:Ne
    plot(ts_out.*1e3,VeD_1(Ne-ii+1,:).*1e3./LFPmax + ze(ii).*1e3 + 0.5, ...
        'color', 'r', 'clipping','on', 'linewidth', 1.5)
end
hold on
plot((599-10+80)-[20 40], [-0.03; -0.03] + 0.5, '-k', 'LineWidth', 1.5);  % time scale
hold on
plot((599-10+80)-[20 20], [-3e-2 ; 2e-2]./LFPmax+ ze(1).*1e3 + 0.5, '-k', ...
    'LineWidth', 1.5); % voltage scale
hold off
text((599-10+80)-mean([20 40]), -0.08 + 0.5, '20 ms', ...
    'HorizontalAlignment','center', 'fontsize', font,'FontWeight','bold') % time scale label
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
ylim([-0.05 1.7] + 0.5)
xlim([(599-10) (599-10+80)])
box 'off'
set(ax,'XTickLabel',[],'XTick',[],'XColor',[1 1 1],'FontWeight','bold',...
    'FontSize',12,'LineWidth',2)

%% figure: threshold stim LFP

figure;
LFPmax = 7*max(abs(VeD_2.*1e3),[],'all'); % scale for the LFP plot
plot(ts_out.*1e3, (ze.*1e3.*ones(length(ts_out), length(ze)))' + 0.5,'-', ...
    'Color', [0.5 0.5 0.5], 'linewidth', 1)
hold on
for ii = 1:Ne
    plot(ts_out.*1e3,VeD_2(Ne-ii+1,:).*1e3./LFPmax + ze(ii).*1e3 + 0.5, ...
        'color', 'r', 'clipping','on', 'linewidth', 1.5)
end
hold on
plot((599-10+80)-[20 40], [-0.03; -0.03] + 0.5, '-k', 'LineWidth', 1.5);  % time scale
hold on
plot((599-10+80)-[20 20], [-3e-2 ; 2e-2]./LFPmax+ ze(1).*1e3 + 0.5, '-k', ...
    'LineWidth', 1.5); % voltage scale
hold off
text((599-10+80)-mean([20 40]), -0.08 + 0.5, '20 ms', ...
    'HorizontalAlignment','center', 'fontsize', font,'FontWeight','bold') % 
% time scale label
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
ylim([-0.05 1.7] + 0.5)
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

[zs_2,CSD_cs_2] = make_cubic_splines(el_pos,VeD_2,Fcs); % compute CSD
if ~isempty(gauss_sigma) && gauss_sigma~=0 % filter iCSD
    [zs_2,CSD_cs_2] = gaussian_filtering(zs_2,CSD_cs_2,gauss_sigma,...
        filter_range);
end
iCSD_2 = CSD_cs_2; % [nA/mm3] current source density
zs_2 = zs_2*1e3; % depth in m, convert to mm

%%
save(fullfile(pathData,'simCSD.mat'), "zs_1","iCSD_1", "zs_2","iCSD_2")

%% ---- plot CSD map
%% figure: supra-CF stim

tspan = (590:(1e3/Fs_new):699) - 590; % generate time span of the data
font = 14; % font size
figure;
imagesc(tspan, zs_1, iCSD_1(:,590:(1e3/Fs_new):699).*1e-3);
yticks(0.05:0.2:1.65)
c = colorbar;
colormap(jet);
c.Label.String = '\muA/mm^{3}';
max_CSD = max(abs(iCSD_1.*1e-3), [], 'all');
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

%% figure: threshold stim

figure;
imagesc(tspan, zs_2, iCSD_2(:,590:(1e3/Fs_new):699).*1e-3);
yticks(0.05:0.2:1.65)
c = colorbar;
colormap(jet);
c.Label.String = '\muA/mm^{3}';
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

%% calculate the dipolar moment from the CSD
d_CSD_1 = (-(zs_1.*1e-3 - (((h*Ne)/2)+a))'.*mean(diff(zs_1.*1e-3)))'*...
    iCSD_1(:,590:(1e3/Fs_new):699)*(pi*(rc^2)).*1e9; % A*m -> nA*m
% baseline correction
d_CSD_1 = d_CSD_1 - mean(d_CSD_1(1:9));

d_CSD_2 = (-(zs_2.*1e-3 - (((h*Ne)/2)+a))'.*mean(diff(zs_2.*1e-3)))'*...
    iCSD_2(:,590:(1e3/Fs_new):699)*(pi*(rc^2)).*1e9; % A*m -> nA*m
% baseline correction
d_CSD_2 = d_CSD_2 - mean(d_CSD_2(1:9));

%% calculate the dipolar moment from the transmembrane currents

vertInd_right = 3831;
vertInd_left = 3786;

% load dipolar moments (column coordinates system)
path2dipoles = []; % 
load(fullfile(path2dipoles, 'column_eq_dip.mat'),'d');
d_1 = d;
load(fullfile(path2CellsData_2, 'dipoles','column_eq_dip.mat'),'d');
d_2 = d;

% load dipoles' normal
load(fullfile(pathLeadFieldIN, 'sources_coordinates', ...
    ['coordinatesVert' num2str(vertInd_right) ...
    '_n#' num2str(0) '.mat']), ...
    'xyzCortNorm');
xyzCortNormR = xyzCortNorm; clearvars xyzCortNorm
load(fullfile(pathLeadFieldIN, 'sources_coordinates', ...
    ['coordinatesVert' num2str(vertInd_left) ...
    '_n#' num2str(0) '.mat']), ...
    'xyzCortNorm');
xyzCortNormL = xyzCortNorm; clearvars xyzCortNorm

% calculate rotation matrices for each dipole
[Rx_right,Ry_right] = calRotMatriz(xyzCortNormR);
[Rx_left,Ry_left] = calRotMatriz(xyzCortNormL);

% rotate dipolar moments to the head coordinate system
dR_In = Rx_right'*Ry_right'*d_2;
dL_In = Rx_left'*Ry_left'*d_1;

%% downsample d_In for comparison with d_CSD
% -- downsampling the LFPs to 1kHz
if ~brainstorm('status')
    brainstorm nogui
end
Fs_new = 1e3; % 1kHz new sampling frequency
[dRs_In, ~] = process_resample('Compute', dR_In, ts*1e-3, Fs_new);
[dLs_In, ts_out] = process_resample('Compute', dL_In, ts*1e-3, Fs_new);

%% visual comparison of the dipoles

figure;
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
    xyzCortNormL*dLs_In(:,590:(1e3/Fs_new):699) - mean(xyzCortNormL* ...
    dLs_In(:,590:(1e3/Fs_new):600)), '-k','LineWidth',1)
hold on;
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
    d_CSD_1, '-r','LineWidth',1)
ylabel({'Current Dipole','Moment (nA*m)'})
xlabel('Time (ms)')
legend({'In','CSD'},'Box','off')
box off
xlim([0 110])
set(gca,'fontweight','bold','FontSize',12,'LineWidth',2)

MAG(xyzCortNormL*dLs_In(:,590:(1e3/Fs_new):699) - mean(xyzCortNormL* ...
    dLs_In(:,590:(1e3/Fs_new):600)),d_CSD_1)
RDM(xyzCortNormL*dLs_In(:,590:(1e3/Fs_new):699) - mean(xyzCortNormL* ...
    dLs_In(:,590:(1e3/Fs_new):600)),d_CSD_1)

figure;
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
xyzCortNormR*dRs_In(:,590:(1e3/Fs_new):699) - mean(xyzCortNormR* ...
dRs_In(:,590:(1e3/Fs_new):600)), '-k','LineWidth',1)
hold on;
plot((ts_out(:,590:(1e3/Fs_new):699)-0.589).*1e3, ...
d_CSD_2, '-r','LineWidth',1)
ylabel({'Current Dipole','Moment (nA*m)'})
xlabel('Time (ms)')
legend({'In','CSD'},'Box','off')
box off
xlim([0 110])
set(gca,'fontweight','bold','FontSize',12,'LineWidth',2)

%% calculate the EEG produced by each dipole

% load the lied field field matrix
filename = ['leadFieldDipBEMVert' num2str(vertInd_left) '.mat'];
load(fullfile(path2leadField_V4, filename));

% calculate average reference
Ne = length(Keoo(:,1));
Hn = eye(Ne) - (ones(Ne,1)*ones(Ne,1)')/Ne;

% change lead fields reference to the average reference
KeDip = Hn*Keoo;

% calculate EEGs
EEG_CSD = KeDip*[xyzCortNormL'*d_CSD_1; xyzCortNormR'*d_CSD_2].*1e-3; % uV
EEG_Dip = KeDip*[dLs_In(:,590:(1e3/Fs_new):699) - ...
    mean(dLs_In(:,590:(1e3/Fs_new):600),2);...
    dRs_In(:,590:(1e3/Fs_new):699) - ...
    mean(dRs_In(:,590:(1e3/Fs_new):600),2)].*1e-3; % uV

%% save 
save(fullfile(pathData,'EEG_DipIn.mat'),'EEG_Dip','dR_In', 'dL_In')
save(fullfile(pathData,'EEG_DipCSD.mat'),'EEG_CSD','d_CSD_1','d_CSD_2')

%% calculate EEG In (point sources per neuron)

EEG_In_R = 0;
EEG_In_L = 0;

for ii = 0:999

    sprintf('neuron %d',ii)

    % load lead field for point sources in neuron 'ii'
    load(fullfile(pathLeadFieldIN,'leadField', ...
        ['leadFieldMonBEMVert' num2str(vertInd_right) ...
        'C' num2str(ii) '.mat']));
    KeooR =Keoo; clearvars Keoo

    load(fullfile(pathLeadFieldIN,'leadField', ...
        ['leadFieldMonBEMVert' num2str(vertInd_left) ...
        'C' num2str(ii) '.mat']));
    KeooL =Keoo; clearvars Keoo

    % change lead field reference to avg reference
    Ne = length(KeooR(:,1));
    Hn = eye(Ne) - (ones(Ne,1)*ones(Ne,1)')/Ne;

    KeR = Hn*KeooR;
    KeL = Hn*KeooL;

    % load neuron 'ii' transmembrane currents
    file = fullfile(path2CellsData_1, ['NeuronsData_r1_n#' num2str(ii) '.mat']);
    load(file, 'It') % nA
    It_1 = It; 
    file = fullfile(path2CellsData_2, ['NeuronsData_r1_n#' num2str(ii) '.mat']);
    load(file, 'It') % nA
    It_2 = It;

    % calculate EEG for each hemisphere and add them
    EEG_InR = KeR*It_2;
    EEG_InL = KeL*It_1;
    EEG_In_R = EEG_In_R + EEG_InR; % nV
    EEG_In_L = EEG_In_L + EEG_InL; % nV

    % save EEGs
    save(fullfile(pathLeadFieldIN, 'EEG_In_V4_Fig2', ...
        ['EEGIn_Vert' num2str(vertInd_right) 'C' num2str(ii) '.mat']), ...
        'EEG_InR');
    save(fullfile(pathLeadFieldIN, 'EEG_In_V4_Fig2', ...
        ['EEGIn_Vert' num2str(vertInd_left) 'C' num2str(ii) '.mat']), ...
        'EEG_InL');

end

%% downsample EEG_In
if ~brainstorm('status')
    brainstorm nogui
end
[EEG_InL_ds, ~] = process_resample('Compute', EEG_In_L, ts*1e-3, Fs_new);
[EEG_InR_ds, ~] = process_resample('Compute', EEG_In_R, ts*1e-3, Fs_new);

%%  select window of interest 
EEG_InR_ds = EEG_InR_ds(:,590:(1e3/Fs_new):699).*1e-3; % nV -> uV
EEG_InL_ds = EEG_InL_ds(:,590:(1e3/Fs_new):699).*1e-3; % nV -> uV
EEG_In_ds = EEG_InR_ds + EEG_InL_ds; % uV

%% save EEG_In
save(fullfile(pathData,'EEG_In.mat'),'EEG_In_ds','EEG_InL_ds','EEG_InR_ds',...
    'EEG_In_L','EEG_In_R')

%% filter EEG at 100 Hz

EEG_In_f = eegfilt(EEG_In_ds, Fs_new, 0, 100);
EEG_Dip_f = eegfilt(EEG_Dip, Fs_new, 0, 100);
EEG_CSD_f = eegfilt(EEG_CSD, Fs_new, 0, 100);

%% save filtered data
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

mean(MAG(EEG_In_f - mean(EEG_In_f(:,1:10),2),EEG_Dip))
mean(RDM(EEG_In_f - mean(EEG_In_f(:,1:10),2),EEG_Dip))

mean(MAG(EEG_In_f - mean(EEG_In_f(:,1:10),2),EEG_CSD))
mean(RDM(EEG_In_f - mean(EEG_In_f(:,1:10),2),EEG_CSD))

mean(MAG(EEG_Dip,EEG_CSD))
mean(RDM(EEG_Dip,EEG_CSD))

