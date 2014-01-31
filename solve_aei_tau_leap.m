function [t_out,s_out] = solve_aei_tau_leap(settings,species,reactions)
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
    
    % Set seed for random number generator
    rng(settings.seed);
    
    % Settings specific to adaptive explicit-implicit tau-leaping algorithm
    na = settings.na;
    nb = settings.nb;
    nc = settings.nc;
    nd = settings.nd;
    delta = settings.delta;
    epsilon = settings.epsilon;

    % Calculate number of species
    N = length(structfun(@isempty,species));

    % Calculate number of reactions
    M = length(reactions);

    % Preallocate and initialize t and s
    t_store = [tstart; zeros(tsteps-1,1)];
    s_store = [species_unpacker(species); zeros(tsteps-1,N)];
    t = tstart;
    s_current = species;
    step = 1;
    
    % Counts down from nb for single-reaction ssa
    ssa_countdown = 0;
    
    % Store h(i) and n(i) for each species
    species_names = fieldnames(species);
    for i = 1:N
        [his.(species_names{i}), nis.(species_names{i})] = highest_order(species,reactions);
    end

    % Main solver loop
    while (t<tend && step<tsteps)

        % Get reaction propensities
        a = get_propensities(s_current,0,V);

        % Propensity normalization
        a0 = sum(a);
        
        % Decide between modes
        
        % Get critical reactions (rnxs are ordered/numbered in here)
        crits = get_critical_rxns(a,s_current,nc);
        
        % Compute candidate step sizes
        % Explicit tau-leaping gets step size from noncritical rxns
        noncrits = array_invert_pos(crits);
        
        % Implicit tau-leaping gets step size from rxns that are
        % noncritical and not in partial equilibrium
        
        
        % SSA Mode
        if ssa_countdown > 0
            % Calculate waiting time tau with 1 random number
            tau = (1/a0)*log(1/rand);

            % Find index u of next reaction such that a1 + a2+ ... + au-1 < ao*r2 <
            % a1 + a2 + ..... + au
            % a is first normalized to create a vector of probabilities for each rxn
            u = pick_rxn(a/a0);

            % Update species
            s_current = get_updates(s_current,u);
            
            % Increment time
            t = t + tau;
            
            % Decrement ssa counter
            ssa_countdown = ssa_countdown - 1;
        end
        
        % Explicit tau-leaping mode
        
        
        % Implicit tau-leaping mode        
        
        
        

        % Increment step
        step = step + 1;

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
    
function crits = get_critical_rxns(a,s,reactions,nc)
    % Inputs:
    %   a : vector of propensities
    %   s : vector of current counts of species
    %   reactions : struct array of reactions
    %   nc : threshold
    % Outputs:
    %   crits : array of critical rxn indices
    M = length(reactions);
    for i = 1:M
        % Get reactant names and stoichiometries
        reactants = fieldnames(reactions(i).reactants);
        num_reactants = length(reactants);
        reactant_stoichs = structfun(@sum,reactions(i).reactants);
        
        dists = zeros(num_reactants); % preallocate
        for j = 1:num_reactants
            dists(j) = floor(s.(reactants{j})/reactant_stoichs(j));
        end
        min_dist = min(dists);
        
        if a(i) > 0 && min_dist < nc
            crits = [crits i];
        end
    end

function inv_arr = array_invert_pos(arr,N)
    % Returns vector of positions inv_arr not included in arr, total size N
    total_arr = 1:N;
    idxs = ismember(total_arr,arr);
    inv_arr = total_arr(~idxs);
    
function g = bound_change(s,species,his,nis)
    %   s : state
    %   species : current species
    sigma = arrayfun(@(j) j/(s.(species)-j), 1:nis.(species)-1);
    g = his.(species) + his.(species)/nis.(species)*sum(sigma);
    
function [h,n] = highest_order(species,reactions)
    % 
    % h : highest order of reactions where species is a reactant
    % n : max number of species molecules required by any of the
    %     highest-order reactions
    
    rxn_stoichs = []; % stoichiometry pairs
    rxn_orders = []; % order of containing rxn
    
    M = length(reactions);
    for i = 1:M
        reactants = fieldnames(reactions(i).reactants);
        if ismember(species,reactants)
            % Overall stoichiometry of reaction
            stoichs = structfun(@sum,reactions(i).reactants);
            order = sum(stoichs);
            rxn_orders = [rxn_orders; order];
            
            % Stoichiometry of species
            rxn_stoichs = [rxn_stoichs; reactions(i).reactants.(species)];
        end
    end
    
    % Get highest order rxn species is a reactant in
    h = max(rxn_order);
    
    % Get max number of molecules of species required in a reactant of order h
    max_order_idxs = (rxn_orders==h);
    n = max(rxn_stoichs(max_order_idxs));
    
    