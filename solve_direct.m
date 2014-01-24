function [t_out,s_out] = solve_direct(settings,species,species_names,reactions,reaction_names)
    % 
    % species : struct with species names as fields and initial counts as values
    % c       : struct
    % reactions : struct
    % s_out : tsteps x N array of species concentrations over time
    
    % Required settings
    tstart = settings.tstart;
    tend   = settings.tend;
    tsteps = settings.tsteps;
    V = settings.volume;

    % Calculate number of species
    N = length(structfun(@isempty,species));

    % Calculate number of reactions
    M = length(structfun(@isempty,reactions));

    % Preallocate and initialize t and s
    t_store = [tstart; zeros(tsteps-1,1)];
    s_store = [species_unpacker(species); zeros(tsteps-1,N)];
    t = tstart;
    s_current = species;
    step = 1;

    % Main solver loop
    while (t<tend && step<tsteps)

        % Get reaction propensities
        a = get_propensities(s_current,V);

        % Propensity normalization
        a0 = sum(a);

        % Calculate waiting time tau with 1 random number
        tau = (1/a0)*log(1/rand);

        % Find index u of next reaction such that a1 + a2+ ... + au-1 < ao*r2 <
        % a1 + a2 + ..... + au
        % a is first normalized to create a vector of probabilities for each rxn
        u = pick_rxn(a/a0);

        % Update species
        s_current = get_updates(s_current,u);

        % Increment step and time
        step = step + 1;
        t = t + tau;

        % Store species counts and time
        s_store(step,:) = species_unpacker(s_current);
        t_store(step) = t;

    end

    % Output
    t_out = t_store(1:step);
    s_out = s_store(1:step,:);


function state = pick_rxn(probs)
    % Vectorized function to pick a state to jump to based on a list of
    % probabilities. All we're doing is creating a list of bins
    % that a randomly chosen number from 0-1 can fall into, and then
    % searching for that bin using vector operations. Same thing can be
    % accomplished via a for loop.
    %
    % This vectorized version is significantly faster than the double for
    % loop implementation for large reaction lists.
    
    selection = rand;
    
    P = cumsum([0, probs]);
    P(end) = 1;
    
    C = P >= selection;
    D = P < selection;
    
    state = norm(find(~(C(2:end) - D(1:end-1))));
    