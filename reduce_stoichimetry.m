function reactions =  reduce_stoichimetry(reactions)

reaction_names = fieldnames(reactions);
M = length(reaction_names); % number of reactions

for i = 1:M
    
    % Ignore rxns that have either no reactants or no products
    if isempty(reactions.(reaction_names{i}).reactants)
        continue
    elseif isempty(reactions.(reaction_names{i}).products)
        continue
    end
    
    reactants = fieldnames(reactions.(reaction_names{i}).reactants);
    products = fieldnames(reactions.(reaction_names{i}).products);
    
    for j = 1:length(reactants)
        for k = 1:length(products)
            if reactants{j} == products{k}
                while reactions.(reaction_names{i}).reactants.(reactants{j}) > 0 && reactions.(reaction_names{i}).products.(products{k}) > 0
                    reactions.(reaction_names{i}).reactants.(reactants{j}) = reactions.(reaction_names{i}).reactants.(reactants{j}) - 1;
                    reactions.(reaction_names{i}).products.(products{k}) = reactions.(reaction_names{i}).products.(products{k}) - 1;
                end
                
                % Cull 0 stoich reactions
                if reactions.(reaction_names{i}).reactants.(reactants{j}) == 0
                    reactions.(reaction_names{i}).reactants = rmfield(reactions.(reaction_names{i}).reactants, (reactants{j}));
                end
                if reactions.(reaction_names{i}).products.(products{k}) == 0
                    reactions.(reaction_names{i}).products = rmfield(reactions.(reaction_names{i}).products, (reactants{j}));
                end
            end
        end
    end
end
