function get_critical_rxns_gen(reactions)
    % Generates a get_critical_rxns.m function that returns critical rxns
    % (those whose concs. may go negative) for a tau-leaping algorithm
    % Inputs:
    %   a : vector of propensities
    %   s : vector of current counts of species
    %   reactions : struct array of reactions
    %   nc : threshold
    
    fid = fopen('get_critical_rxns.m','w');

    fprintf(fid, 'function crits = get_critical_rxns(a,s,nc)\n');
    
    fprintf(fid, 'crits = [];\n');
    
    M = length(reactions);
    for i = 1:M
        % Get reactant names and stoichiometries
        reactants = fieldnames(reactions(i).reactants);
        reactant_stoichs = structfun(@sum,reactions(i).reactants);
        
        fprintf(fid, 'dists = zeros(%d,1);\n',length(reactants));
        
        for j = 1:length(reactants)
            fprintf(fid, 'dists(%d) = floor(s.%s/%d);\n',j,reactants{j},reactant_stoichs(j));
        end
        
        fprintf(fid, 'min_dist = min(dists);\n');
        
        fprintf(fid, 'if a(%d) > 0 && min_dist < nc\n',i);
        fprintf(fid, '\tcrits = [crits %d];\n',i);
        fprintf(fid, 'end\n');
    end
    
    fclose(fid);