function mag_value = MAG(y_exact, y_num)
%% Computes the Magnitude ratio (MAG) between ana and num estimations
% The MAG reveals the variations in potential amplitude or, in other words,
% alterations in source strength.
% Pursiainen et al. (2016), Phys. Med. Biol. 61  8502–8520
%
% Input:
%   y_ana: analitic, ground-truth or exact solution. Nchannels x
%   Ntime_points
%   y_num: numerical or approximate solution. Nchannels x Ntime_points.
% Output:
%  mag_value: magnitude measurement. 1 x Nchannels.
%
% Author: Beatriz Herrera, October 15, 2020
%
%% calculations

mag_value = (vecnorm(y_num,2,2)./vecnorm(y_exact,2,2));

end