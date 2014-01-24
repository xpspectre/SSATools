function ssa_sim_test

%% Parse input file
[settings,constants,species,reactions] = parse_input('test2test.txt');

%% Initialize
% Generate reaction propensity solver
get_propensities_gen(constants,reactions);

% Generate reaction updating block
get_updates_gen(reactions);

% Initial state of system/molecule counts
s0 = struct2array(species);

%% Simulate
[t,s] = simulate_SSA(settings,constants,species,reactions);

function [t,s] = simulate_SSA(settings,constants,species,reactions)
% Get constants and species names so we don't have to use cell arrays in
% compiled functions
constant_names = char(fieldnames(constants));
species_names = char(fieldnames(species));
reaction_names = char(fieldnames(reactions));

% Call compiled solver
[t_out,s_out] = solve_direct(settings,species,species_names,reactions,reaction_names);

% Convert structs to arrays for output
t = 1;
s = 1;