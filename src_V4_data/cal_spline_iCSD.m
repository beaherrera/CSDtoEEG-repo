%% calculate iCSD maps from laminar recordings
% Author: Beatriz Herrera

clear
clc

%% Parameters
Ne = 16; % number of electrodes in the shank
a = 0.1; % [mm] position of the first electrode
elec_spacing = 0.15; % [mm] electrode spacing
ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % electrode positions with
% respect to the pia surface
% WARNING: position of the first electrode must be different from zero
el_pos = ze*1e-3;  % mm to m
cond = 0.33; %[S/m] gray matter conductance | NOTE: if length(cond)==1, the
% function iCSDspline_Pettersen() considers con_top = cond (conductance at
% the top of the cylinder or above the pia matter)
gauss_sigma = 0.1e-3;   %[m] Gaussian filter std
diam = 3e-3; % [m] cylinder diameter

%%
% Ve == LOCAL FIELD POTENTIALS IN VOLTS
% Outputs: zs ~ mm & iCSD ~ nA/mm3
[zs, iCSD] = iCSDspline_Pettersen(Ve,el_pos,diam,cond,gauss_sigma);
