function dz = cal_dip_CSD(zs, iCSD, rc)
%CAL_DIP_CSD calculate the current dipole moment from the CSD
%   calculates the current dipole moment from the laminar current source
%   density (CSD) assuming the dipole is located at the center of a
%   cylindrical cortical column of 'rc' mm radius, which also corresponds
%   to the center of the solumn's coordinate system. 
%
% Inputs:
%   zs: [mm] depth relative to the pia matter at which the CSD was
%       calculated [Nelectrodes x 1]
%   iCSD: [nA/mm3] current source density [Nelectrodes x Ntimepts]
%   z: cortical depth -> m (zero coordinate at pia matter)
%   rc: [mm] cortical column radius
% 
% Output:
%   dz: [nA*m] current dipole moment amplitude
%
% Author: Beatriz Herrera
%
%%

dz = 1e-3.*(-(zs-median(zs)).*mean(diff(zs)))'*iCSD*(pi*(rc^2)); % nA*m

end