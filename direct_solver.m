function [t,s] = direct_solver(t0,s0,c,t_end,t_steps)

% Calculate number of species
N = length(s0);

% Preallocate t and s
t_store = zeros(tsteps,1);
s_store = zeros(tsteps,N);

% Initialize stored time and species
t_store(1) = t0;
s_store(1,:) = s0;

% Main solver loop
while (t<t_end)
    
end