function [t,s] = solve_direct(s0,tstart,tend,tsteps,V,c)

% Calculate number of species
N = length(s0);

% Preallocate t and s
t_store = zeros(tsteps,1);
s_store = zeros(tsteps,N);

% Initialize stored time and species
t_store(1) = tstart;
s_store(1,:) = s0;

% Main solver loop
while (t<tend && step<tsteps)
    
end