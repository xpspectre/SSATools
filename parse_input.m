function [settings,constants,species,reactions] =  parse_input(filename)
% Parse input file, JSON format
    % Read file into a string and parse
    data = parse_json(fileread(filename));
    
    % Parse Settings
    settings = data.settings;
    
    % Parse Constants
    constants = data.constants;
    
    % Parse Species
    species = data.species;
    
    % Parse Reactions
    reactions_strs = data.reactions;
    r_names = fieldnames(reactions_strs);
    for i = 1:length(r_names)
        rxn = reactions_strs.(r_names{i});
        
        % Get [reactants_str] ->{[rate_str]} [products_str]
        [~,tok_i,~] = regexp(rxn,'(.*?)->\{(.*)}(.*)','match','tokens');
        reactants_str = tok_i{1,1}{1,1};
        rate_str = tok_i{1,1}{1,2};
        products_str = tok_i{1,1}{1,3};
        
        clear reactants products
        
        % Get reactant species and stoichimetries
        [~,tok_j,~] = regexp(reactants_str,'([0-9]*)([A-Za-z0-9_])','match','tokens');
        for j = 1:length(tok_j)
            if strcmpi(tok_j{1,j}{1,1},'') % empty coefficient => 1
                tok_j{1,j}{1,1} = '1';
            end
            reactants.(tok_j{1,j}{1,2}) = str2double(tok_j{1,j}{1,1});
        end
        
        % Get product species and stoichimetries
        [~,tok_j,~] = regexp(products_str,'([0-9]*)([A-Za-z0-9_])','match','tokens');
        for j = 1:length(tok_j)
            if strcmpi(tok_j{1,j}{1,1},'') % empty coefficient => 1
                tok_j{1,j}{1,1} = '1';
            end
            products.(tok_j{1,j}{1,2}) = str2double(tok_j{1,j}{1,1});
        end
        
        % Build reactions struct
        reactions.(r_names{i}).reactants = reactants;
        reactions.(r_names{i}).products = products;
        reactions.(r_names{i}).rate_str = rate_str;
        
    end
