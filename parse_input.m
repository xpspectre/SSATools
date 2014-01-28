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
        rate_str = parts.f_rate(4:end-1);
        if strcmpi(rate_str(1),'''') == 1 && strcmpi(rate_str(end),'''') == 1
            rate_str = rate_str(2:end-1); % strip quotes
            rate_type = 'raw';
        else % default
            rate_type = 'massaction';
        end
        
        % Build  forward reaction struct
        reactions.(r_names{i}).reactants = reactants;
        reactions.(r_names{i}).products = products;
        reactions.(r_names{i}).rate_str = rate_str;
        reactions.(r_names{i}).rate_type = rate_type;
        
        % Process reverse rxn ( of the form: <-{[rate]} ) if present
        %   Assumes a forward rate is always present
        %   Makes a new rxn titled [forwardrxn]r, otherwise the same
        %   [forwardrxn] reactants and products reversed
        if strcmpi(parts.r_rate,'') || isempty(parts.r_rate(4:end-1)) % forward rxn only
            % do nothing
        else % reversible rxn
            rrxn_name = [r_names{i} 'r'];

            rate_str = parts.r_rate(4:end-1);
            if strcmpi(rate_str(1),'''') == 1 && strcmpi(rate_str(end),'''') == 1
                rate_str = rate_str(2:end-1); % strip quotes
                rate_type = 'raw';
            else % default
                rate_type = 'massaction';
            end
            
            reactions.(rrxn_name).reactants = products; % reverse of the fwd rxn
            reactions.(rrxn_name).products = reactants;
            reactions.(rrxn_name).rate_str = rate_str;
            reactions.(rrxn_name).rate_type = rate_type;
            
        end
    end
