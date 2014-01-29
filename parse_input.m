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
    species_list = fieldnames(species); % for checking reactions
    
    % Parse Reactions
    rxn_names = fieldnames(data.reactions);
    num_rxns = length(rxn_names);
    reactions = [];
    for i = 1:num_rxns
        % Add reaction
        reactions = add_reaction(reactions,rxn_names{i},data.reactions.(rxn_names{i}));
        
        % Make sure all species found in reaction are also initialized in
        % species; if not, add and initialize to 0
        if ~isempty(reactions.(rxn_names{i}).reactants)
            reactants = fieldnames(reactions.(rxn_names{i}).reactants);
            for j = 1:length(reactants)
                if ~ismember(reactants{j},species_list)
                    warning('parse_input: species %s not found in species list ... adding',reactants{j})
                    species = add_species(species,reactants{j},0);
                end
            end
        end
        
        if ~isempty(reactions.(rxn_names{i}).products)
            products = fieldnames(reactions.(rxn_names{i}).products);
            for j = 1:length(products)
                if ~ismember(products{j},species_list)
                    warning('parse_input: species %s not found in species list ... adding',products{j})
                    species = add_species(species,products{j},0);
                end
            end
        end
    end
    