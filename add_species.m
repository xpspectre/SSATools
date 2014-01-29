function species = add_species(species,species_name,init_count)
    % Add species and initial count
    % Args: 
    %   species : struct containing species
    %   species_name : string containing new species
    %   init_count : string containing initial count
    % Returns:
    %   species : struct containing species with added species
    
    species.(species_name) = init_count;