function ssa_sim_test_next_rxn

[settings,constants,species,reactions] = parse_input('test2.txt');

% Get reduced stoichiometry reactions
reduced_reactions = reduce_stoichimetry(reactions);

%% Initialize
% Generate reaction propensity solver
get_propensities_gen(constants,reactions);

% Generate reaction updating block
get_updates_gen(reduced_reactions);

% Generate species unpacking (struct->array) function
species_unpacker_gen(species);

% Make dependency graph function
get_dependency_gen(reduced_reactions);

% Decrease simulation tsteps for test function
settings.tsteps = 3;

% Reload directories
rehash

%% Simulate
% Call compiled solver
[~,~] = solve_next_reaction(settings,species,reactions);
