function [t_out,s_out] = solve_first_reaction(settings,species,reactions)
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
        a = get_propensities(s_current,0,V);

        % Calculate M independent random numbers
        taus = (1./a).*log(1./rand(1,M));
        
        % Get first reaction (minimum tau)
        [tau, u] = min(taus);

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
