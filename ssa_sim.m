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

%TEST
reduced = reduce_stoichimetry(reactions);

%% Initialize
% Generate reaction propensity solver
get_propensities_gen(constants,reactions);

% Generate reaction updating block
get_updates_gen(reactions);

% Generate species unpacking (struct->array) function
species_unpacker_gen(species);

% Make dependency graph function
get_dependency_gen(reactions);

%% Simulate
% Call compiled solver
% [t_out,s_out] = solve_direct(settings,species,reactions);
% [t_out,s_out] = solve_first_reaction(settings,species,reactions);
[t_out,s_out] = solve_next_reaction(settings,species,reactions);

plot(t_out,s_out())
legend('A','B','C')

disp(0)
