function [settings,constants,species,reactions] =  parse_input(filename)
    % Parse input file, JSON format
    % Read file into a string and parse
    data = parse_json(fileread(filename));
    
    % Parse Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings = data.settings;
    settings_fields = fieldnames(settings);
    
    % Make a random seed if one isn't supplied
    if ~ismember('seed',settings_fields) || isempty(settings.seed)
        settings.seed = int32(RandStream.shuffleSeed);
    end
    
    % Parse Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constants = data.constants;
    
    % Parse Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    species = data.species;
    species_list = fieldnames(species); % for checking reactions
    
    % Parse Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rxn_names = fieldnames(data.reactions);
    num_rxns = length(rxn_names);
    reactions = [];
    for i = 1:num_rxns
        % Add reaction
        % Get next available rxn_idx
        rxn_idx = length(reactions) + 1;
        reactions = add_reaction(reactions,rxn_idx,data.reactions.(rxn_names{i}));
        
        % Make sure all species found in reaction are also initialized in
        % species; if not, add and initialize to 0
        if ~isempty(reactions(i).reactants)
            reactants = fieldnames(reactions(i).reactants);
            for j = 1:length(reactants)
                if ~ismember(reactants{j},species_list)
                    warning('parse_input: species %s not found in species list ... adding',reactants{j})
                    species = add_species(species,reactants{j},0);
                end
            end
        end
        
        if ~isempty(reactions(i).products)
            products = fieldnames(reactions(i).products);
            for j = 1:length(products)
                if ~ismember(products{j},species_list)
                    warning('parse_input: species %s not found in species list ... adding',products{j})
                    species = add_species(species,products{j},0);
                end
            end
        end
    end
    