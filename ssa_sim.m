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
[settings,constants,species,reactions] = parse_input(filename);

%% Initialize
% Generate reaction propensity solver
get_propensities_gen(constants,reactions);

% Generate reaction updating block
get_updates_gen(reactions);

% Generate species unpacking (struct->array) function
species_unpacker_gen(species);

% Dependency graph
dependency = get_dependency(reactions);

%% Simulate

[t,s] = simulate_SSA(settings,constants,species,reactions);


disp(1)

function [t,s] = simulate_SSA(settings,constants,species,reactions)
% Get constants and species names so we don't have to use cell arrays in
% compiled functions
constant_names = char(fieldnames(constants));
species_names = char(fieldnames(species));
reaction_names = char(fieldnames(reactions));

% Call compiled solver
[t_out,s_out] = solve_direct(settings,species,species_names,reactions,reaction_names);

% Diagnostic Plots
plot(t_out,s_out())

% Convert structs to arrays for output
t = 1;
s = 1;