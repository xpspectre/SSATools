function init_ssatools
% Initialize SSATools
% Just adds the SSATools folder to the path
ssatools_path = fileparts(mfilename('fullpath'));
addpath(ssatools_path);