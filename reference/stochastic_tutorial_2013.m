function stochastic_tutorial_2013
% Sample implementations of Gilespie algorithm, with comparison to
% deterministic solutions

% Conversion reaction: A <--> B ---> C at rates 'k1', 'k2' and 'k3'
% respectively

clc; clear all; close all

tic % start timer

%% Set parameters and run deterministic/stochastic conversion rxn

% Parameters
k1 = 0.01;
k2 = 1.6e-3;
k3 = 5e-4;
params_conv = [k1,k2,k3];

% Initial conditions
A = 100;
B = 100;
C = 0;
y0_conv = [A,B,C];

% Run time
tspan = 1e3;

% Run deterministic
[tdc,ydc] = ode15s(@deterministic_conversion,[0 tspan],y0_conv,[],params_conv);

% Run stochastic
[tsc,ysc] = stochastic_conversion(tspan,params_conv,y0_conv);

% Plot both
figure()
plot(tdc,ydc,tsc,ysc)
xlabel('Time (s)')
ylabel('Species')
legend('A Deterministic','B Deterministic','C Determinisitc','A Stochastic','B Stochastic', 'C Stochastic','Location','best')
ylim([0 200])

toc % end timer

%% ODE function for conversion reaction

function dydt = deterministic_conversion(~,init_conv,params)

% A -> B with rate k1
% B -> A with rate k2

% unpack and collect params
pCell = num2cell(params);
[k1,k2,k3] = pCell{:};

% Unpack and collect species 
yCell = num2cell(init_conv);
[A,B,C] = yCell{:};

% ODEs
dA = -k1*A + k2*B;
dB = k1*A - k2*B -k3*B;
dC = k3*B;

dydt = [dA;dB;dC];


%% Stochastic conversion 
function [t,y] = stochastic_conversion(tspan,params,init)

% Rate constants (numbered 1,2,...M)
M = length(params);

% unpack and collect params
pCell = num2cell(params);
[k1,k2,k3] = pCell{:};

% Number of reactant molecules (numbered 1,2,...N)
N = length(init);

% Unpack and collect species 
yCell = num2cell(init);
[A,B,C] = yCell{:};

% Initialize time
t = 0;

% Set up storage variables for output 
count = 2; % count is a counter for storage variables, count = 1 corresponds to the initial values

% Pre-allocate values in matrices to reduce run time 
%(adding rows and columns to matrices on-the-fly will increase run time)
tstore = zeros(1e6,1); 
Xstore = zeros(1e6,N);

% Assign initial values to time and species
tstore(1) = t;
Xstore(1,:) = init;

% Run iteration
while (t<=tspan) 
    % Calculate au values
   
    a(1) = k1*A;  %a1
    a(2) = k2*B;  %a2
    a(3) = k3*B;  %a3
    
    % Calculate ao: ao = sum(a1,a2,a3...aM)
    ao = sum(a);
    
    % Generate random numbers r1 and r2, uniformly distributed on [0 1]
    % rand function pulls from uniform distribution on [0 1]
    r1 = rand(1,1);
    r2 = rand(1,1);
    
    % Calculate waiting time, tau
    tau = (1/ao)*log(1/r1);
    
    % Find index u of next reaction such that a1 + a2+ ... + au-1 < ao*r2 <
    % a1 + a2 + ..... + au
    
   
    yr2 = ao*r2;
    
    for k = 1:M
        sum1 = sum(a(1:k-1));
        sum2 = sum(a(1:k));
        
        if yr2 > sum1 && yr2 <= sum2
            u = k;
        end
    end
    
    % Update species
    if u==1
        A = A-1;
        B = B+1;
    elseif u==2
        B = B-1;
        A = A+1;
    elseif u ==3
        B = B-1;
        C = C+1;
       
    end
    
    % Store outputs
    tstore(count) = t;
    Xstore(count,:) = [A;B;C];
    
    % Update time and counter
    t = t+tau;
    count = count + 1;

end

% Assign final outputs
t = tstore(1:count-1);
y = Xstore(1:count-1,:);