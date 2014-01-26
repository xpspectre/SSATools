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
    
    % TEST
%     tdiff = zeros(tsteps,1);
    
    % Get all reaction propensities
    a = get_propensities(s_current,0,V); % note updated syntax
    
    % Get putative times
    taus = (1./a).*log(1./rand(1,M)) + t; % include starting t for absolute times
    
    % Mark whether a reaction was 0 previously
    was_zero = (a==0);
    % Mark last time reaction was nonzero (irrevelent for nonzero propensity rxns)
    last_nonzero_t = was_zero*t;
    % Mark last propensity
    last_nonzero_a = was_zero*0;
    
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
        for n = affected_nodes % select reaction
            
            % Get new propensities
            an_old = a(n);
            a(n) = get_propensities(s_current,n,V);
            
            % Change times
            if n ~= u
                if an_old == 0 && a(n) ~= 0 % switched back to nonzero propensity
                    disp('ghaghfawhfoawghw')
%                     was_zero(n) = 0;
                    t_last_nonzero = last_nonzero_t(n);
                    a_last_nonzero = last_nonzero_a(n);
                    tau_a = a_last_nonzero/a(n) * (t - t_last_nonzero) + t;
                elseif an_old ~= 0 && a(n) == 0 % becomes zero propensity
%                     was_zero(n) = 1;
                    last_nonzero_t(n) = t;
                    last_nonzero_a(n) = an_old;
                    tau_a = Inf;
                else
                    tau_a = an_old/a(n) * (ipq.get_rxn(n) - t) + t;
                end
                %TEST
%                 tdiff(step)=tau_a<t;
            else % n == u
                tau_a = 1/a(n)*log(1/rand) + t;
                %TEST
%                 tdiff(step)=tau_a<t;
            end
            ipq.update_rxn(n,tau_a);
        end
        
    end
    
    %TEST
%     find(tdiff<0)

    % Output
    t_out = t_store(1:step);
    s_out = s_store(1:step,:);
