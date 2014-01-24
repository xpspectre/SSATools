function species_unpacker_gen(species)
% Create a file that unpacks species in struct from to species in array form

fid = fopen('species_unpacker.m','w');

fprintf(fid, 'function s_array = species_unpacker(s_struct)\n');

species_list = char(fieldnames(species));
num_species = size(species_list,1);

fprintf(fid, 's_array = zeros(1,%d);\n',num_species);

for i = 1:num_species
    fprintf(fid, 's_array(%d) = s_struct.%s;\n',i,species_list(i,:));
end

fclose(fid);