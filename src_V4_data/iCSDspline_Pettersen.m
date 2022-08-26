function [zs, iCSD] = iCSDspline_Pettersen(Ve,el_pos,diam,cond,gauss_sigma)
% FUNCTION ICSDSPLINE_PETTERSEN calculate the CSD using spline iCSD method
% 
% method by Pettersen et al. 2006. Method implementation obtained from
% the iCSDplotter toolbox.
% written by Beatriz Herrera, 2019.
%
% Inputs:
%   - Ve: local field potentials [Nchannels x Ntimepts] in volts
%   - el_pos: electrodes position relative to the pia matter [meters]
%   - diam: diameter of the cylinder considered to calculate the CSD [meters]
%   - cond: conductance in S/m. If length(cond)==1 cond_top=cond= grey
%   matter conductance. otherwise [cond, cond_top] = [grey matter 
%   conductance, conductance at the top of the cylinder]
%   - gauss_sigma: Gaussian filter standard deviation in meters for the
%   spatial smoothing. 
%
% Outputs:
%   - zs: [mm] depth relative to the pia matter at which the CSD was
%   calculated
%   - iCSD: [nA/mm3] CSD matrix [Nelectrodes x Ntimepts] 
%
%% Parameters
if length(cond) > 1
    cond_1 = cond(1); % [S/m] grey matter conductance
    cond_top = cond(2); %[S/m] conductance at the top (cylinder)
else
    cond_1 = cond;
    cond_top = cond;
end
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent

%% solve Pettersen model
Fcs = F_cubic_spline(el_pos,diam,cond_1,cond_top);

[zs,CSD_cs] = make_cubic_splines(el_pos,Ve,Fcs);
[zs,CSD_cs] = gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range); %
% current source density and electrodes' position
iCSD = CSD_cs; % [nA/mm3] current source density
zs = zs*1e3; % m -> mm

end

