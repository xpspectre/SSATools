function ssa_sim_test

%% Parse input file
[settings,constants,species,reactions] = parse_input('test.txt');

% Get reduced stoichiometry reactions
reduced_reactions = reduce_stoichimetry(reactions);

%% Initialize
% Generate reaction propensity solver
get_propensities_gen(constants,reactions);

% Generate reaction updating block
get_updates_gen(reactions);

% Generate species unpacking (struct->array) function
species_unpacker_gen(species);

% Decrease simulation tsteps for test function
settings.tsteps = 3;

%% Simulate
% Call compiled solver
[t_out,s_out] = solve_direct(settings,species,reactions);
