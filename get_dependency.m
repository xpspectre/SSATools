function dependency =  get_dependency(reactions)
% Dependency graph of reactions
% Store dependency graph as a struct:
%   dependency :
%       Reaction 1 :
%           DependsOn : Reaction 2
%           Affects :  Reaction N
%       Reaction 2 :
%           Depends On : Reaction N
%           Affects : 
% Self edges included trivially

rxns = fieldnames(reactions); % list or reactions
M = length(rxns);

% Reduce rxns to stoichimetric changes only - so A+B->A+C doesn't change A
% - low priority

% Augment reactions struct with DependsOn and Affects
for i = 1:M
    rxn_i = rxns{i};
    reactions.(rxn_i).DependsOn = fieldnames(reactions.(rxn_i).reactants);
    reactions.(rxn_i).Affects = unique([reactions.(rxn_i).DependsOn; fieldnames(reactions.(rxn_i).products)]);
end

% Preallocate dependency graph
for i = 1:M
    rxn_i = rxns{i};
    dependency.(rxn_i).In = {};
    dependency.(rxn_i).Out = {};
end

% Populate dependency graph
for i = 1:M
    rxn_i = rxns{i};
    for j = 1:M % Cycle thru all rxns' Affects
        rxn_j =rxns{j};
        
        % rxn_i dependent species found in rxn_j affected species?
        if sum(ismember(reactions.(rxn_i).DependsOn,reactions.(rxn_j).Affects))
            dependency.(rxn_i).In = [dependency.(rxn_i).In; rxn_j]; % Depends on species affected by these rxns
            dependency.(rxn_j).Out = [dependency.(rxn_j).Out; rxn_i]; % Affects species that these rxns depend on
        end
    end
end

% Convert dependency graph to struct of char arrays
for i = 1:M
    rxn_i = rxns{i};
    dependency.(rxn_i).In = char(dependency.(rxn_i).In);
    dependency.(rxn_i).Out = char(dependency.(rxn_i).Out);
end
