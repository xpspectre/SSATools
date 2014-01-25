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
        
        reactants = [];
        products = [];
        rxn = reactions_strs.(r_names{i});
                
        % Get [reactants_str] ->{[rate_str]} [products_str]
        p1 = '^(?<reactants>.*?)';
        p2 = '(?<f_rate>->\{[^}]*\})';
        p3 = '(?<r_rate><-\{[^}]*\})?';
        p4 = '(?<products>.*?)$';
        expr = [p1 p2 p3 p4];
        parts = regexp(rxn,expr,'names');
        
        % Get reactant species and stoichimetries
        tok = regexp(parts.reactants,'([0-9]*)([A-Za-z0-9_])*','tokens');
        for j = 1:length(tok)
            if strcmpi(tok{1,j}{1,1},'') % empty coefficient => stoich 1
                tok{1,j}{1,1} = '1';
            end
            reactants.(tok{1,j}{1,2}) = str2double(tok{1,j}{1,1});
        end
        
        % Get product species and stoichimetries
        tok = regexp(parts.products,'([0-9]*)([A-Za-z0-9_])*','tokens');
        for j = 1:length(tok)
            if strcmpi(tok{1,j}{1,1},'') % empty coefficient => stoich 1
                tok{1,j}{1,1} = '1';
            end
            products.(tok{1,j}{1,2}) = str2double(tok{1,j}{1,1});
        end
        
        % Process forward rxn ( of the form: ->{[rate]} )
        %   Assumes a forward rate is always present
        f_rate_str = parts.f_rate(4:end-1);
        if strcmpi(f_rate_str(1),'''') == 1 && strcmpi(f_rate_str(end),'''') == 1
            f_rate_type = 'raw';
        else % default
            f_rate_type = 'massaction';
        end
        
        % Process reverse rate ( of the form: <-{[rate]} ) if present
        %   Assumes a forward rate is always present
        if strcmpi(parts.r_rate,'') || isempty(parts.r_rate(4:end-1)) % forward rxn only
            r_rate_str = '';
            r_rate_type = 'none';
        else % reversible rxn
            r_rate_str = parts.r_rate(4:end-1);
            if strcmpi(r_rate_str(1),'''') == 1 && strcmpi(r_rate_str(end),'''') == 1
                r_rate_type = 'raw';
            else % default
                r_rate_type = 'massaction';
            end
        end
        
        % Build reactions struct
        reactions.(r_names{i}).reactants = reactants;
        reactions.(r_names{i}).products = products;
        reactions.(r_names{i}).f_rate_str = f_rate_str;
        reactions.(r_names{i}).f_rate_type = f_rate_type;
        reactions.(r_names{i}).r_rate_str = r_rate_str;
        reactions.(r_names{i}).r_rate_type = r_rate_type;
    end
