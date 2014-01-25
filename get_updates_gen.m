function get_updates_gen(reactions)
% Create a file that calculates the species updates

fid = fopen('get_updates.m','w');

fprintf(fid, 'function s = get_updates(s,u)\n');

fprintf(fid, 'switch u\n');

reaction_names = fieldnames(reactions);
M = length(reaction_names); % number of reactions

for i = 1:M
    fprintf(fid, '\tcase %d\n',i);
    
    reactants = fieldnames(reactions.(reaction_names{i}).reactants);
    reactant_stoichs = structfun(@sum,reactions.(reaction_names{i}).reactants);
    products = fieldnames(reactions.(reaction_names{i}).products);
    product_stoichs = structfun(@sum,reactions.(reaction_names{i}).products);
    
    % Decrement reactants
    for j = 1:length(reactants) % number of reactants
        fprintf(fid, '\t\ts.%s = s.%s - %d;\n',reactants{j},reactants{j},reactant_stoichs(j));
    end
    
    % Increment products
    for j = 1:length(products) % number of products
        fprintf(fid, '\t\ts.%s = s.%s + %d;\n',products{j},products{j},product_stoichs(j));
    end
end

fprintf(fid, 'end\n');

fclose(fid);