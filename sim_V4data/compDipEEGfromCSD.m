%% Compute the dipolar moment and the resulting EEG from the CSD
% author: Beatriz Herrera
% Distributed under GPL-3.0 License
%   Created: November 9th, 2020

clear
clc

%% Load current source densities of contra and ipsi lateral conditions
[namefile, pathfile] =uigetfile('*.mat',['File -> CSD matrix (nA/mm3)' ...
    ' (Channels x TimePoints)']);
load([pathfile namefile], 'contra', 'ipsi'); clearvars namefile pathfile
CSD_contra = contra.*1e-3; % nA/mm3 -> uA/mm3
CSD_ipsi = ipsi.*1e-3; % nA/mm3 -> uA/mm3
Ne = length(CSD_contra(:,1)); % total number of electrodes

%% -- enter cortical parameters
opts.Interpreter = 'tex';
definput = {'0.1','3','0.1'}; % default input values
Parameters = inputdlg({['Enter the inter-electrodes distance in mm ' ...
    '(must correspond with that of the CSD matrix)'], ...
    'Enter the diameter of the cortical column in mm',...
    ['Enter the relative distance of the first electrode with respect' ...
    ' to the pial matter in mm']}...
    ,'',[1,90],definput,opts);
h = str2double(Parameters{1}); % inter-electrodes distance, convert from
% string to double
rc = (str2double(Parameters{2}))/2; % cortical column radius, convert from
% string to double
a = str2double(Parameters{3}); % inter-electrodes distance, convert from
% string to double
zs = (a:h:((Ne-1)*h + a))'; % position of the electrodes along z with
% respect to pia matter, convert from string to double 

%% Compute Dipolar Moment
disp('Estimating the current dipole moment')
d_CSDcontra = 1e-6.*(-(zs-median(zs)).*h)'*CSD_contra*(pi*(rc^2)); % mA*m
d_CSDipsi = 1e-6.*(-(zs-median(zs)).*h)'*CSD_ipsi*(pi*(rc^2)); % mA*m

% visualize the current dipole moment
figure;
plot(d_CSDcontra,'-b'); hold on;
plot(d_CSDipsi, '-r')
legend('Contra', 'Ipsi')
xlabel('Time-Points')
ylabel('Current Dipole Moment [mA*m]')

%% Compute the EEG
opts.Interpreter = 'tex';
definput = {'left'};
hemisphere = inputdlg({['Enter the hemisphere where the contralateral dipolar '...
    'moment will be located. (Options: left | right)']}...
    ,'',[1,90],definput,opts);

% concatenate the current dipole moments based on the location of the
% target.
switch hemisphere{1}
    case 'left'
        de = [d_CSDcontra; d_CSDipsi];
    case 'right'
        de = [d_CSDipsi; d_CSDcontra];
    otherwise
        
end

opts.Interpreter = 'tex';
definput = {'1'};
Param = inputdlg({['Enter the number of dipole locations (lead field matrices)'... 
    ' that will be used to estimate the EEG.']}...
    ,'',[1,90],definput,opts);
numKe = str2double(Param{1});

% -- lead field matrix
for ii=1:numKe
    [namefile, pathfile] =uigetfile('*.mat',sprintf('File -> Lead Field Matrix # %d', ii));
    Ke = load([pathfile namefile]); clearvars namefile pathfile
    varName = fieldnames(Ke);
    Ke = Ke.(varName{1});
    
    EEG = Ke*de; % mV
    
    outMatUnits = {'mA*m', 'mV'};
    [namefile, pathfile] = uiputfile('*.mat',['Select folder path to save '...
        'the EEG and dipolar moment matrices and enter file name.']);
    save([pathfile namefile],'de', 'EEG', 'outMatUnits')
    
    clearvars Ke EEG
end
