function get_updates_gen(reactions)
% Create a file that calculates the species updates

fid = fopen('get_updates.m','w');

fprintf(fid, 'function s = get_updates(s,u)\n');

fprintf(fid, 'switch u\n');

reaction_names = char(fieldnames(reactions));
M = length(structfun(@isempty,reactions)); % number of reactions

for i = 1:M
    fprintf(fid, '\tcase %d\n',i);
    
    reactants = char(fieldnames(reactions.(reaction_names(i,:)).reactants));
    products = char(fieldnames(reactions.(reaction_names(i,:)).products));
    
    % Decrement reactants
    for j = 1:size(reactants,1) % number of reactants
        fprintf(fid, '\t\ts.%s = s.%s - 1;\n',reactants(j,:),reactants(j,:));
    end
    
    % Increment products
    for j = 1:size(products,1) % number of products
        fprintf(fid, '\t\ts.%s = s.%s + 1;\n',products(j,:),products(j,:));
    end
end

fprintf(fid, 'end\n');

fclose(fid);