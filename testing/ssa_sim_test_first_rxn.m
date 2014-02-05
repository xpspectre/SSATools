function ssa_sim_test_first_rxn

%% Parse input file
[settings,constants,species,reactions] = parse_input('test.txt');

% Get reduced stoichiometry reactions
reduced_reactions = reduce_stoichimetry(reactions);

%% Initialize
% Generate reaction propensity solver
get_propensities_gen(constants,reactions);

% Generate reaction updating block
get_updates_gen(reduced_reactions);

% Generate species unpacking (struct->array) function
species_unpacker_gen(species);

% Decrease simulation tsteps for test function
settings.tsteps = 3;

% Reload directories
rehash

%% Simulate
% Call compiled solver
[~,~] = solve_first_reaction(settings,species,reactions);
