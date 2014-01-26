function get_propensities_gen(constants,reactions)
% Create a file that calculates the reaction propensities
% OR CHANGE TO GENERATE ALL NEEDED FILES, can change with method

fid = fopen('get_propensities.m','w');

fprintf(fid, 'function a = get_propensities(s,u,V)\n');

constant_names = char(fieldnames(constants));
reaction_names = char(fieldnames(reactions));
L = length(structfun(@isempty,constants)); % number of constants
M = length(structfun(@isempty,reactions)); % number of reactions

% Assign all the constants
for i = 1:L
    fprintf(fid,'%s = %f;\n',constant_names(i,:),constants.(constant_names(i,:)));
end

% Return a single propensity
fprintf(fid, 'switch u\n');
for i = 1:M
    fprintf(fid, 'case %d\n',i);
    fprintf(fid, 'a = ');
    print_propensity(fid,i,reactions);
end

% Otherwise print entire propensity array (next_reaction initialization or other methods each step)
fprintf(fid, 'otherwise\n');
fprintf(fid,'a = zeros(1,%d);\n',M); % Preallocate reaction propensity array

for i = 1:M
    fprintf(fid,'a(%d) = ',i); % reaction number
    print_propensity(fid,i,reactions);
end

fprintf(fid, 'end\n');

fclose(fid);

function print_propensity(fid,i,reactions)
    reaction_names = fieldnames(reactions);
    rate_str = reactions.(reaction_names{i}).rate_str;
    
    % 0th order rxns will have empty reactant fields
    if isempty(reactions.(reaction_names{i}).reactants)
        fprintf(fid,'%s * V;\n',rate_str);
        return
    end
    
    % Higher order rxns
    reactants = fieldnames(reactions.(reaction_names{i}).reactants);
    stoichs = structfun(@sum,reactions.(reaction_names{i}).reactants);
    order = sum(stoichs);
    switch order
        case 1
            fprintf(fid,'s.%s * %s;\n',reactants{1},rate_str);
        case 2
            if size(unique(reactants),1)==order % 2 different molecules reacting
                fprintf(fid,'s.%s * s.%s * %s / V;\n',reactants{1},reactants{2},rate_str);
            else % 2 same molecules reacting
                fprintf(fid,'s.%s * (s.%s-1) * %s / V;\n',reactants{1},reactants{1},rate_str);
            end
        otherwise % ignore higher order rxns
    end


