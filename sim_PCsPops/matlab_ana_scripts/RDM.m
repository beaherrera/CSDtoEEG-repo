function rdm_value = RDM(y_exact, y_num)
%% Computes the Relative Difference (RDM) between ana and num estimations
% The RDM reflects the topographical forward modeling error 
% in terms of location and orientation.
% Pursiainen et al. (2016), Phys. Med. Biol. 61  8502–8520
%
% Input:
%   y_ana: analitic, ground-truth or exact solution. Nchannels x
%   Ntime_points
%   y_num: numerical or approximate solution. Nchannels x Ntime_points.
% Output:
%  rdm_value: relative difference between measurements. 1 x Nchannels.
%  Values between 0 and 2.
%
% Author: Beatriz Herrera, October 15, 2020
%
%% calculations

rdm_value = vecnorm((normalize(y_num,2,'norm') - normalize(y_exact,2,'norm')),2,2);

end