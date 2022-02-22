%% program compute the CSD using the spline-iCSD method from Pettersen et al. 2006
% Copyright 2020-, Beatriz Herrera 
% Distributed under GPL-3.0 License
%   Created: November 9th, 2020

clear
clc

%% Input files
% -- local field potential
[namefile, pathfile] =uigetfile('*.mat',['File -> Local Field ' ...
    'Potential Matrix (Channels x TimePoints)']);
lfp = load([pathfile namefile]); clearvars namefile pathfile
varName = fieldnames(lfp);
lfp = lfp.(varName{1});
% -- asking for the units
unitsLFP = inputdlg(['Enter the units of the LFP (Options: mV,' ...
    ' uV and V)'],'',[1,50]);
unitsLFP = unitsLFP{1};
% convert LFP to appropriate units if necessary
switch unitsLFP
    case 'mV'
        lfp = lfp.*1e-3; % mV -> V
    case 'uV'
        lfp = lfp.*1e-6; % uV -> V
    case 'V'
        
    otherwise
        warning('Please enter one of the following units: uV, mV or V.')
end

% -- values of z (depth)
[namefile, pathfile] =uigetfile('*.mat',['File -> z-coordinates' ...
    ' to compute the CSD (Channels x 1)']);
ze = load([pathfile namefile]); clearvars namefile pathfile
varName = fieldnames(ze);
ze = ze.(varName{1});
% -- asking for the units
unitsz = inputdlg('Enter the units of z (Options: um, mm or m)','',[1,50]);
unitsz = unitsz{1};
% convert z to appropriate units if necessary
switch unitsz
    case 'mm'
        ze = ze.*1e-3; % mm -> m
    case 'um'
        ze = ze.*1e-6; % um -> m
    case 'm'
        
    otherwise
        warning('Please enter one of the following units: um, mm or m.')
end
% -- iCSD parameters
opts.Interpreter = 'tex';
definput = {'0.1','3','0.33','0.1'}; % default parameters
Parameters = inputdlg({'Enter the inter-electrodes distance in mm', ...
    'Enter the diameter of the cortical column in mm', ...
    'Enter the grey matter conductance in S/m',...
    'Enter the standard deviation of the Gaussian Filter (\sigma) in mm'}...
    ,'Spline iCSD parameters',[1,90],definput,opts);
Ne = length(ze); % number of electrodes in the shank
el_pos = ze; % depth coordinates in m
h = str2double(Parameters{1}).*1e-3; % inter-electrodes distance, convert mm -> m
% ze = a:h:((Ne-1)*h + a); % [m] electrode positions with respect to the pia surface
diam = (str2double(Parameters{2})*1e-3); % [m] cylinder diameter, convert mm -> m
cond = str2double(Parameters{3}); %[S/m] gray matter conductance
cond_top = cond; 
gauss_sigma = str2double(Parameters{3})*1e-3;   %[m] Gaussian filter std, convert mm -> m
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent

%% ------------  compute CSD
%% solve Pettersen model
Fcs = F_cubic_spline(el_pos,diam,cond,cond_top); % cubic spline function

[zs,CSD_cs] = make_cubic_splines(el_pos,lfp,Fcs); % compute CSD
if ~isempty(gauss_sigma) && gauss_sigma~=0 % filter iCSD
    [zs,CSD_cs] = gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
end
iCSD = CSD_cs; % [nA/mm3] current source density -> multiply by .*1e-3 
                % for converting to uA/mm3
zs = zs*1e3; % depth in m, convert to mm
outMatUnits = {'mm','nA/mm3'}; % output units for the CSD and depth

%% save CSD matrix

[namefile, pathfile] = uiputfile('*.mat',['Select folder path to save' ...
    ' CSD matrix and enter file name.']);
save([pathfile namefile],'zs', 'iCSD', 'outMatUnits')

%% Visualization
visualization = inputdlg(['Do you want to plot the CSD map? Enter: yes -> y or '...
    'no -> n.'],'',[1,50]);
visualization = visualization{1};
switch visualization
    case 'y'
        timeSpan = inputdlg({'Enter sampling rate in Hz', ['Pre-stimulus' ...
            ' time window in ms'],...
            'Post-stimulus time window in ms'},'',[1,50]);
        Fs = str2double(timeSpan{1});
        preWin = str2double(timeSpan{2});
        postWin = str2double(timeSpan{3});
    otherwise
        disp('Results saved in output folder.')
end

tspan = preWin:(1e-3/Fs):postWin; % generate time span of the data
tspan = tspan(1:end-1);
font = 14; % font size
figure;
imagesc(tspan, zs, iCSD);
c = colorbar;
colormap(jet);
c.Label.String = 'nA/mm^{3}';
c.Label.FontSize = font;
max_CSD = max(abs(iCSD), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
ax = gca;
ax.FontSize = font;
xlabel('Time (ms)', 'FontSize',font);
ylabel('z (mm)','FontSize',font);
