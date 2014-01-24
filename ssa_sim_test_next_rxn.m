function ssa_sim_test_next_rxn

%% Parse input file
[settings,constants,species,reactions] = parse_input('test.txt');

%% Initialize
% Generate reaction propensity solver
get_propensities_gen(constants,reactions);

% Generate reaction updating block
get_updates_gen(reactions);

% Generate species unpacking (struct->array) function
species_unpacker_gen(species);

% Make dependency graph function
get_dependency_gen(reactions);

% Decrease simulation tsteps for test function
settings.tsteps = 3;

%% Simulate
% Call compiled solver
[t_out,s_out] = solve_next_reaction(settings,species,reactions);
