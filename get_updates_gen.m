function get_updates_gen(reactions)
    % Create a file that calculates the species updates

    fid = fopen('get_updates.m','w');

    fprintf(fid, 'function s = get_updates(s,u)\n');

    fprintf(fid, 'switch u\n');

    M = length(reactions); % number of reactions

    for i = 1:M
        fprintf(fid, '\tcase %d\n',i);

        % Decrement reactants
        if isempty(reactions(i).reactants) % 0th order/no reactants
            % skip decrementing reactions
        else
            reactants = fieldnames(reactions(i).reactants);
            reactant_stoichs = structfun(@sum,reactions(i).reactants);
            for j = 1:length(reactants) % number of reactants
                fprintf(fid, '\t\ts.%s = s.%s - %d;\n',reactants{j},reactants{j},reactant_stoichs(j));
            end
        end

        % Increment products
        if isempty(reactions(i).products) % degradations/no products
            % skip incrementing reactions
        else
            products = fieldnames(reactions(i).products);
            product_stoichs = structfun(@sum,reactions(i).products);
            for j = 1:length(products) % number of products
                fprintf(fid, '\t\ts.%s = s.%s + %d;\n',products{j},products{j},product_stoichs(j));
            end
        end
    end

    fprintf(fid, 'end\n');

    fclose(fid);