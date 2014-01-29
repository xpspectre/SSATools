function reactions = add_reaction(reactions,rxn_name,rxn_str)
    % Add reaction (and its reverse, if applicable) to reactions struct
    % Args: 
    %   reactions : struct containing reactions
    %   rxn_name : string containing reaction name
    %   rxn_str : string containing reaction
    % Returns:
    %   reactions : struct containing reactions with added rxn

    % Get [reactants_str] ->{[rate_str]} [products_str]
    p1 = '^(?<reactants>.*?)';
    p2 = '(?<f_rate>->\{[^}]*\})';
    p3 = '(?<r_rate><-\{[^}]*\})?';
    p4 = '(?<products>.*?)$';
    expr = [p1 p2 p3 p4];
    parts = regexp(rxn_str,expr,'names');
    
    % Initialize reactant and product structs
    reactants = [];
    products = [];

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
    reactions.(rxn_name).reactants = reactants;
    reactions.(rxn_name).products = products;
    reactions.(rxn_name).rate_str = rate_str;
    reactions.(rxn_name).rate_type = rate_type;

    % Process reverse rxn ( of the form: <-{[rate]} ) if present
    %   Assumes a forward rate is always present
    %   Makes a new rxn titled [forwardrxn]r, otherwise the same
    %   [forwardrxn] reactants and products reversed
    if strcmpi(parts.r_rate,'') || isempty(parts.r_rate(4:end-1)) % forward rxn only
        % do nothing
    else % reversible rxn
        rrxn_name = [rxn_name 'r'];

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