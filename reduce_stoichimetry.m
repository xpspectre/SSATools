function reactions =  reduce_stoichimetry(reactions)

M = length(reactions); % number of reactions

for i = 1:M
    
    % Ignore rxns that have either no reactants or no products
    if isempty(reactions(i).reactants) || isempty(reactions(i).products)
        continue
    end
    
    reactants = fieldnames(reactions(i).reactants);
    products = fieldnames(reactions(i).products);
    
    for j = 1:length(reactants)
        for k = 1:length(products)
            if reactants{j} == products{k}
                while reactions(i).reactants.(reactants{j}) > 0 && reactions(i).products.(products{k}) > 0
                    reactions(i).reactants.(reactants{j}) = reactions(i).reactants.(reactants{j}) - 1;
                    reactions(i).products.(products{k}) = reactions(i).products.(products{k}) - 1;
                end
                
                % Cull 0 stoich reactions
                if reactions(i).reactants.(reactants{j}) == 0
                    reactions(i).reactants = rmfield(reactions(i).reactants, (reactants{j}));
                end
                if reactions(i).products.(products{k}) == 0
                    reactions(i).products = rmfield(reactions(i).products, (reactants{j}));
                end
            end
        end
    end
end
