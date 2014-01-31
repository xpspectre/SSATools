function get_dependency_gen(reactions)
% Accessible dependency graph
% Takes in a reduced reaction struct

fid = fopen('get_dependency.m','w');

fprintf(fid, 'function affected_nodes = get_dependency(u)\n');

M = length(reactions); % number of reactions

% Augment reactions struct with DependsOn and Affects
for i = 1:M
    % DependsOn reactants
    if isempty(reactions(i).reactants) % no reactants
        reactions(i).DependsOn = {};
    else
        reactions(i).DependsOn = fieldnames(reactions(i).reactants);
    end
    
    % Affects DependsOn and products
    if isempty(reactions(i).products) % no products
        reactions(i).Affects = reactions(i).DependsOn;
    else
        reactions(i).Affects = unique([reactions(i).DependsOn; fieldnames(reactions(i).products)]);
    end
    
end

% Preallocate dependency graph
for i = 1:M
    dependency(i).In = {};
    dependency(i).Out = {};
end

% Populate dependency graph
for i = 1:M
    for j = 1:M % Cycle thru all rxns' Affects
        % rxn_i dependent species found in rxn_j affected species?
        % Add index of reaction
        if ~isempty(ismember(reactions(i).DependsOn,reactions(j).Affects))
            dependency(i).In = [dependency(i).In; j]; % Depends on species affected by these rxns
            dependency(j).Out = [dependency(j).Out; i]; % Affects species that these rxns depend on
        end
    end
end
% dep = struct2cell(dependency); % convert rxn names to indices in a cell array


% Write cases to file
fprintf(fid, 'switch u\n');
for i = 1:M
    fprintf(fid, '\tcase %d\n',i);
    affected_nodes = cell2mat(dependency(i).Out);
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
