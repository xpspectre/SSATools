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
    
%     % Indicates whether a propensity is 0
%     zero_p = (a==0);
%     became_zero = ones(1,length(a))*tstart;
%     
%     % Time since any propensitiy went to 0, Inf if nonzero
%     t_since_0 = zeros(1,length(a));
    
    % Get putative times
    taus = (1./a).*log(1./rand(1,M));
    
    % Store putative times in indexed priority queue
    ipq = Indexed_Priority_Queue([1:M;taus]);

    % Main solver loop
    while (t<tend && step<tsteps)
        
        % Get smallest tau and corresponding reaction u
        [u,tau] = ipq.get_minimum();
        
        %Set time
        t = tau;
        
        % Update species
        s_current = get_updates(s_current,u);

        % Get dependency
        affected_nodes = get_dependency(u);
        
        % Update affected nodes
        for i = 1:length(affected_nodes)
            
            n = affected_nodes(i);
            
            % Get new propensities
            an_old = a(n);
            a(n) = get_propensities(s_current,n,V);
            
            % Change times
            if n ~= u
                if a(n) == 0
%                     zero_p(n) = 1; % note that rxn n propensity is 0
%                     
                    tau_a = Inf;
                else
                    tau_a = an_old/a(n) * (ipq.get_rxn(n) - t) + t;
                end
            else % n == u
                tau_a = 1/a(n)*log(1/rand) + t;
            end
            ipq.update_rxn(n,tau_a);
        end
        
        % Diagnostics
%         sprintf('%d\t%f\n',step,t)
        
        % Increment step and time
        step = step + 1;
%         t = t + tau;

        % Store species counts and time
        s_store(step,:) = species_unpacker(s_current);
        t_store(step) = t;

    end

    % Output
    t_out = t_store(1:step);
    s_out = s_store(1:step,:);
