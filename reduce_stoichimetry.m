function reduced =  reduce_stoichimetry(reactions)

reaction_names = fieldnames(reactions);
M = length(reaction_names); % number of reactions

for i = 1:M
    reactants = fieldnames(reactions.(reaction_names{i}).reactants);
    reactant_stoichs = structfun(@sum,reactions.(reaction_names{i}).reactants);
    products = fieldnames(reactions.(reaction_names{i}).products);
    product_stoichs = structfun(@sum,reactions.(reaction_names{i}).products);
    
    for j = 1:length(reactants)
        for k = 1:length(products)
            if reactants
        end
    end
    
    disp(1)
    
end

% for species i in reactants
%     for species j in products
%         if species i == species j 
%             while both nonzero
%             decrement stoichs for both by 1
        