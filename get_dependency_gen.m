function get_dependency_gen(reactions)
% Create a file that calculates the reaction propensities
% OR CHANGE TO GENERATE ALL NEEDED FILES, can change with method

fid = fopen('get_dependency.m','w');

fprintf(fid, 'function affected_nodes = get_dependency(u)\n');

reaction_names = fieldnames(reactions);
M = length(reaction_names); % number of reactions

% Reduce rxns to stoichimetric changes only - so A+B->A+C doesn't change A
% - low priority - also need to do this in get_updates

% Augment reactions struct with DependsOn and Affects
for i = 1:M
    rxn_i = reaction_names{i};
    
    % DependsOn
    if isempty(reactions.(rxn_i).reactants) % no reactants
        reactions.(rxn_i).DependsOn = {};
    else
        reactions.(rxn_i).DependsOn = fieldnames(reactions.(rxn_i).reactants);
    end
    
    % Affects
    if isempty(reactions.(rxn_i).products) % no products
        reactions.(rxn_i).Affects = unique(reactions.(rxn_i).DependsOn);
    else
        reactions.(rxn_i).Affects = unique([reactions.(rxn_i).DependsOn; fieldnames(reactions.(rxn_i).products)]);
    end
    
end

% Preallocate dependency graph
for i = 1:M
    rxn_i = reaction_names{i};
    dependency.(rxn_i).In = {};
    dependency.(rxn_i).Out = {};
end

% Populate dependency graph
for i = 1:M
    rxn_i = reaction_names{i};
    for j = 1:M % Cycle thru all rxns' Affects
        rxn_j =reaction_names{j};
        
        % rxn_i dependent species found in rxn_j affected species?
        % Add index of reaction
        if ~isempty(ismember(reactions.(rxn_i).DependsOn,reactions.(rxn_j).Affects))
            dependency.(rxn_i).In = [dependency.(rxn_i).In; j]; % Depends on species affected by these rxns
            dependency.(rxn_j).Out = [dependency.(rxn_j).Out; i]; % Affects species that these rxns depend on
        end
    end
end
dep = struct2cell(dependency); % convert rxn names to indices in a cell array


% Write cases to file
fprintf(fid, 'switch u\n');
for i = 1:M
    fprintf(fid, '\tcase %d\n',i);
    affected_nodes = cell2mat(dep{i}.Out);
    fprintf(fid, '\t\taffected_nodes = [%d',affected_nodes(end)); % minimum 1 affected node, itself
    for j = 1:length(affected_nodes)-1
        fprintf(fid, ',%d',affected_nodes(j));
    end
    fprintf(fid, '];\n');
end
fprintf(fid, '\totherwise\n');
fprintf(fid, '\t\taffected_nodes = [];\n');
fprintf(fid, 'end\n');

fclose(fid);
