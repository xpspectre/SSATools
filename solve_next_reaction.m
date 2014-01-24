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
    a = get_propensities(s_current,V);
    
    % Get putative times
    taus = (1./a).*log(1./rand(1,M));
    
    % Store putative times in indexed priority queue
    ipq = Indexed_Priority_Queue([1:M;taus]);

    % Main solver loop
    while (t<tend && step<tsteps)
        
        % Get smallest tau
        [u,tau] = ipq.get_minimum();
        
        %Set time
        t_current = tau;
        
        % Update species
        s_current = get_updates(s_current,u);

        % Get dependency
        affected_nodes = get_dependency(u);
        
        % Update propensitites for affected nodes
%         for i = 1:length(affected_nodes)
%             n = affected_nodes(i);
%             a(n) = 
%         end
        
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
