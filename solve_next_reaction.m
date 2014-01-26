function [t_out,s_out] = solve_next_reaction(settings,species,reactions)
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

    % Get all reaction propensities
    a = get_propensities(s_current,0,V); % note updated syntax
    
    % Get putative times
    taus = (1./a).*log(1./rand(1,M)) + t; % include starting t for absolute times
    
    % Store putative times in indexed priority queue
    ipq = Indexed_Priority_Queue([1:M;taus]);

    % Main solver loop
    while (t<tend && step<tsteps)
        
        % Get smallest tau_u and corresponding reaction u
        [u,tau_u] = ipq.get_minimum();
        
        % Determine new state
        t = tau_u; % update time
        s_current = get_updates(s_current,u); % update species
        step = step + 1;
        
        % Store state and time
        s_store(step,:) = species_unpacker(s_current);
        t_store(step) = t;

        % Get affected reactions
        affected_nodes = get_dependency(u);
        
        % Update affected nodes
        for i = 1:length(affected_nodes) % select reaction
            
            n = affected_nodes(i);
            
            % Get new propensities
            an_old = a(n);
            a(n) = get_propensities(s_current,n,V);
            
            % Change times
            % Big thanks to the creators of Dizzy (http://magnet.systemsbiology.net/software/Dizzy/)
            % I had to look at their source code to figure out what happens
            % when a rxn goes from 0 propensity to non-0 propensity as a
            % result of another rxn. It's not very clear what Gibson and
            % Bruck's paper or any of the written sources do from their
            % texts. Instead, I got a solution from
            % \src\org\systemsbiology\chem\SimulatorStochasticGibsonBruck.java,
            % updateReactionRateAndTime function. They just generated a new
            % random time...
            if n ~= u
                if a(n) > 0 && an_old > 0 % non-0 propensities
                    tau_a = an_old/a(n) * (ipq.get_rxn(n) - t) + t;
                else % became 0 propensity or was 0 propensity 
                    if a(n) > 0 % switched to non-0 propensity
                        tau_a = 1/a(n)*log(1/rand) + t;
                    else % remains/became 0 propensity
                        tau_a = Inf;
                    end
                end
            else % n == u
                tau_a = 1/a(n)*log(1/rand) + t;
            end
            ipq.update_rxn(n,tau_a);
        end
        
    end
    
    % Output
    t_out = t_store(1:step);
    s_out = s_store(1:step,:);
