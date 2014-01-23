function ssa_sim(filename)

% Variables:
%   N : number of species
%   M : number of reactions
%   t : tsteps x 1 array of times
%   s : tstepx x N array of species, s(i,:) is an 1 X N array of species at
%       time i
%   c : struct of constants
%   

%% Parse input file
[settings,c,species,reactions] = parse_input(filename);

%% Initialize

% Required settings
tstart = settings.tstart;
tend   = settings.tend;
tsteps = settings.tsteps;
V      = settings.volume;

% Initial state of system/molecule counts
s0 = struct2array(species);

% Generate propensity functions
%   Rates are assumed by mass action kinetics unless the arbitrary rate
%   flag is set on a reaction

%% Simulate
% [t,s] = direct_solver(s0,tstart,tend,tsteps,c)


disp(1)