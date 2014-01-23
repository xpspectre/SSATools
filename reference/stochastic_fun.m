function [t,s] = stochastic_fun(tend,p,s0)
% function should be compiled as a mex file

    % Number of reactant molecules (numbered 1,2,...N)
    N = length(pack_species(s0));

    % Initialize time
    t = 0;
    
    % Set up storage variables for output 
    time_idx = 2; % time_idx is a counter for storage variables, time_idx = 1 corresponds to the initial values

    % Pre-allocate values in matrices to reduce run time 
    % (adding rows and columns to matrices on-the-fly will increase run time)
    %
    % Now, at the end of the Gillespie simulation, the algorithm prints the
    % length of a species concentration vector to the command window. If
    % this size is greater than 1e7, you need to increase the value in the
    % zeros() commands in order to speed up the simulation.
    
    tstore = zeros(8e6,1); 
    Xstore = zeros(8e6,N);
    
    % Assign initial values to time and species
    tstore(1) = t;
    Xstore(1,:) = pack_species(s0);
    s = s0;

    % Run iteration
    while (t<=tend) 
    
        num_rxns = 75; % For preallocation of reaction rate vector
        
        % Calculate au values
        % assume constant volume = 1
        a = zeros(1,num_rxns);
        
        fX = p.gamma/(p.ce+s.X); % ssrA degradation
        
        a(1) = p.ka  * s.p_00_a * s.a2; % AraC dimers (a2) binding to activator plasmid
        a(2) = p.k_a * s.p_10_a;
        a(3) = p.ka  * s.p_01_a * s.a2;
        a(4) = p.k_a * s.p_11_a;
        a(5) = p.ka  * s.p_02_a * s.a2;
        a(6) = p.k_a * s.p_12_a;
        
        a(7)  = p.ka  * s.p_00_r * s.a2; % AraC dimers (a2) binding to repressor plasmid
        a(8)  = p.k_a * s.p_10_r;
        a(9)  = p.ka  * s.p_01_r * s.a2;
        a(10) = p.k_a * s.p_11_r;
        a(11) = p.ka  * s.p_02_r * s.a2;
        a(12) = p.k_a * s.p_12_r;
        
        a(13) = 2*p.kr * s.p_00_a * s.r4; % LacI tetramers (r4) binding to activator plasmid
        a(14) = p.k_r  * s.p_01_a;
        a(15) = 2*p.kr * s.p_10_a * s.r4;
        a(16) = p.k_r  * s.p_11_a;
        
        a(17) = 2*p.kr * s.p_00_r * s.r4; % LacI tetramers (r4) binding to repressor plasmid
        a(18) = p.k_r  * s.p_01_r;
        a(19) = 2*p.kr * s.p_10_r * s.r4;
        a(20) = p.k_r  * s.p_11_r;
        
        a(21) = p.kr    * s.p_01_a * s.r4; % LacI tetramers (r4) binding to activator plasmid
        a(22) = 2*p.k_r * s.p_02_a;
        a(23) = p.kr    * s.p_11_a * s.r4;
        a(24) = 2*p.k_r * s.p_12_a;
        
        a(25) = p.kr    * s.p_01_r * s.r4; % LacI tetramers (r4) binding to repressor plasmid
        a(26) = 2*p.k_r * s.p_02_r;
        a(27) = p.kr    * s.p_11_r * s.r4;
        a(28) = 2*p.k_r * s.p_12_r;
        
        a(29) = p.kl  * s.p_12_a; % looped promoter binding on activator plasmid
        a(30) = p.kl  * s.p_02_a;
        a(31) = p.kul * s.p_L0_a;
        
        a(32) = p.kl  * s.p_12_r; % looped promoter binding on repressor plasmid
        a(33) = p.kl  * s.p_02_r;
        a(34) = p.kul * s.p_L0_r;
        
        a(35) = p.ba * s.p_00_a; % transcription of activator plasmid
        a(36) = p.alpha*p.ba * s.p_10_a;
        
        a(37) = p.br * s.p_00_r; % transcription of repressor plasmid
        a(38) = p.alpha*p.br * s.p_10_r;
        
        a(39) = p.ta * s.ma; % translation
        a(40) = p.tr * s.mr;
        
        a(41) = p.kfa * s.auf; % folding
        a(42) = p.kfr * s.ruf;
        
        a(43) = p.kda  * s.a * (s.a-1); % multimerization
        a(44) = p.k_da * s.a2;
        a(45) = p.kdr  * s.r * (s.r-1);
        a(46) = p.k_dr * s.r2;
        a(47) = p.kt  * s.r2 * (s.r2-1);
        a(48) = p.k_t * s.r4;
        
        a(49) = p.da * s.ma; % free species decay
        a(50) = p.dr * s.mr;
        a(51) = p.lambda*fX * s.auf;
        a(52) = fX * s.ruf;
        a(53) = p.lambda*fX * s.a;
        a(54) = fX * s.r;
        a(55) = p.lambda*fX * s.a2;
        a(56) = fX * s.r2;
        a(57) = fX * s.r4;
        
        a(58) = fX * s.p_10_a; % promoter bound species decay on activator plasmid
        a(59) = fX * s.p_11_a;
        a(60) = fX * s.p_12_a;
        a(61) = fX * s.p_01_a;
        a(62) = fX * s.p_11_a;
        a(63) = 2*fX * s.p_02_a;
        a(64) = 2*fX * s.p_12_a;
        a(65) = 2*p.epsilon*fX * s.p_L2_a;
        a(66) = p.epsilon*fX * s.p_L1_a;
        
        a(67) = fX * s.p_10_r; % promoter bound species decay on repressor plasmid
        a(68) = fX * s.p_11_r;
        a(69) = fX * s.p_12_r;
        a(70) = fX * s.p_01_r;
        a(71) = fX * s.p_11_r;
        a(72) = 2*fX * s.p_02_r;
        a(73) = 2*fX * s.p_12_r;
        a(74) = 2*p.epsilon*fX * s.p_L2_r;
        a(75) = p.epsilon*fX * s.p_L1_r;

        
        % Sanity check simulation: no rate equations should be negative
%         if any(a < 0)
%             s
%             error('One or more rate equations in the simulation is negative.');
%         end
        
        % Calculate ao: ao = sum(a1,a2,a3...aM)
        ao = sum(a);

        % Generate random number r1, uniformly distributed on [0 1] to
        % calculate waiting time.
        % rand function pulls from uniform distribution on [0 1]
        r1 = rand;

        % Calculate waiting time, tau
        tau = (1/ao)*log(1/r1);

        % Find index u of next reaction such that a1 + a2+ ... + au-1 < ao*r2 <
        % a1 + a2 + ..... + au
        % a is first normalized to create a vector of probabilities for each rxn
        u = norm(pick_rxn(a / ao));
        
        % Update species
        %
        % I've found that including the rate equation in a comment above
        % the species helps to keep track of what is going on.
        
        switch u 
            case 1 % a(1) = p.ka  * s.p_00_a * s.a2;
                s.p_00_a = s.p_00_a - 1;
                s.a2     = s.a2     - 1;
                s.p_10_a = s.p_10_a + 1;
            case 2 % a(2) = p.k_a * s.p_10_a;
                s.p_00_a = s.p_00_a + 1;
                s.a2     = s.a2     + 1;
                s.p_10_a = s.p_10_a - 1;
            case 3 % a(3) = p.ka  * s.p_01_a * s.a2;
                s.p_01_a = s.p_01_a - 1;
                s.a2     = s.a2     - 1;
                s.p_11_a = s.p_11_a + 1;
            case 4 % a(4) = p.k_a * s.p_11_a;
                s.p_01_a = s.p_01_a + 1;
                s.a2     = s.a2     + 1;
                s.p_11_a = s.p_11_a - 1;
            case 5 % a(5) = p.ka  * s.p_02_a * s.a2;
                s.p_02_a = s.p_02_a - 1;
                s.a2     = s.a2     - 1;
                s.p_12_a = s.p_12_a + 1;
            case 6 % a(6) = p.k_a * s.p_12_a;
                s.p_02_a = s.p_02_a + 1;
                s.a2     = s.a2     + 1;
                s.p_12_a = s.p_12_a - 1;
                
            case 7 % a(7) = p.ka  * s.p_00_r * s.a2;
                s.p_00_r = s.p_00_r - 1;
                s.a2     = s.a2     - 1;
                s.p_10_r = s.p_10_r + 1;
            case 8 % a(8) = p.k_r * s.p_10_r;
                s.p_00_r = s.p_00_r + 1;
                s.a2     = s.a2     + 1;
                s.p_10_r = s.p_10_r - 1;
            case 9 % a(9) = p.ka  * s.p_01_r * s.a2;
                s.p_01_r = s.p_01_r - 1;
                s.a2     = s.a2     - 1;
                s.p_11_r = s.p_11_r + 1;
            case 10 % a(10) = p.k_r * s.p_11_r;
                s.p_01_r = s.p_01_r + 1;
                s.a2     = s.a2     + 1;
                s.p_11_r = s.p_11_r - 1;
            case 11 % a(11) = p.ka  * s.p_02_r * s.a2;
                s.p_02_r = s.p_02_r - 1;
                s.a2     = s.a2     - 1;
                s.p_12_r = s.p_12_r + 1;
            case 12 % a(12) = p.k_r * s.p_12_r;
                s.p_02_r = s.p_02_r + 1;
                s.a2     = s.a2     + 1;
                s.p_12_r = s.p_12_r - 1;
                
            case 13 % a(13) = 2*p.kr * s.p_00_a * s.r4;
                s.p_00_a = s.p_00_a - 1;
                s.r4     = s.r4     - 1;
                s.p_01_a = s.p_01_a + 1;
            case 14 % a(14) = p.k_r  * s.p_01_a;
                s.p_00_a = s.p_00_a + 1;
                s.r4     = s.r4     + 1;
                s.p_01_a = s.p_01_a - 1;
            case 15 % a(15) = 2*p.kr * s.p_10_a * s.r4;
                s.p_10_a = s.p_10_a - 1;
                s.r4     = s.r4     - 1;
                s.p_11_a = s.p_11_a + 1;
            case 16 % a(16) = p.k_r  * s.p_11_a;
                s.p_10_a = s.p_10_a + 1;
                s.r4     = s.r4     + 1;
                s.p_11_a = s.p_11_a - 1;
                
            case 17 % a(17) = 2*p.kr * s.p_00_r * s.r4;
                s.p_00_r = s.p_00_r - 1;
                s.r4     = s.r4     - 1;
                s.p_01_r = s.p_01_r + 1;
            case 18 % a(18) = p.k_r  * s.p_01_r;
                s.p_00_r = s.p_00_r + 1;
                s.r4     = s.r4     + 1;
                s.p_01_r = s.p_01_r - 1;
            case 19 % a(19) = 2*p.kr * s.p_10_r * s.r4;
                s.p_10_r = s.p_10_r - 1;
                s.r4     = s.r4     - 1;
                s.p_11_r = s.p_11_r + 1;
            case 20 % a(20) = p.k_r  * s.p_11_r;
                s.p_10_r = s.p_10_r + 1;
                s.r4     = s.r4     + 1;
                s.p_11_r = s.p_11_r - 1;
                
            case 21 % a(21) = p.kr    * s.p_01_a * s.r4;
                s.p_01_a = s.p_01_a - 1;
                s.r4     = s.r4     - 1;
                s.p_02_a = s.p_02_a + 1;
            case 22 % a(22) = 2*p.k_r * s.p_02_a;
                s.p_01_a = s.p_01_a + 1;
                s.r4     = s.r4     + 1;
                s.p_02_a = s.p_02_a - 1;
            case 23 % a(23) = p.kr    * s.p_11_a * s.r4;
                s.p_11_a = s.p_11_a - 1;
                s.r4     = s.r4     - 1;
                s.p_12_a = s.p_12_a + 1;
            case 24 % a(24) = 2*p.k_r * s.p_12_a;
                s.p_11_a = s.p_11_a + 1;
                s.r4     = s.r4     + 1;
                s.p_12_a = s.p_12_a - 1;
                
            case 25 % a(25) = p.kr    * s.p_01_r * s.r4;
                s.p_01_r = s.p_01_r - 1;
                s.r4     = s.r4     - 1;
                s.p_02_r = s.p_02_r + 1;
            case 26 % a(26) = 2*p.k_r * s.p_02_r;
                s.p_01_r = s.p_01_r + 1;
                s.r4     = s.r4     + 1;
                s.p_02_r = s.p_02_r - 1;
            case 27 % a(27) = p.kr    * s.p_11_r * s.r4;
                s.p_11_r = s.p_11_r - 1;
                s.r4     = s.r4     - 1;
                s.p_12_r = s.p_12_r + 1;
            case 28 % a(28) = 2*p.k_r * s.p_12_r;
                s.p_11_r = s.p_11_r + 1;
                s.r4     = s.r4     + 1;
                s.p_12_r = s.p_12_r - 1;
                
            case 29 % a(29) = p.kl  * s.p_12_a;
                s.p_12_a = s.p_12_a - 1;
                s.p_L2_a = s.p_L2_a + 1;
                s.a2     = s.a2     + 1;
            case 30 % a(30) = p.kl  * s.p_02_a;
                s.p_02_a = s.p_02_a - 1;
                s.p_L2_a = s.p_L2_a + 1;
            case 31 % a(31) = p.kul * s.p_L0_a;
                s.p_L0_a = s.p_L0_a - 1;
                s.p_00_a = s.p_00_a + 1;
                
            case 32 % a(32) = p.kl  * s.p_12_r;
                s.p_12_r = s.p_12_r - 1;
                s.p_L2_r = s.p_L2_r + 1;
                s.a2     = s.a2     + 1;
            case 33 % a(33) = p.kl  * s.p_02_r;
                s.p_02_r = s.p_02_r - 1;
                s.p_L2_r = s.p_L2_r + 1;
            case 34 % a(34) = p.kul * s.p_L0_r;
                s.p_L0_r = s.p_L0_r - 1;
                s.p_00_r = s.p_00_r + 1;
                
            case 35 % a(35) = p.ba * s.p_00_a;
                s.ma = s.ma + 1;
            case 36 % a(36) = p.alpha*p.ba * p_10_a;
                s.ma = s.ma + 1;
                
            case 37 % a(37) = p.ba * s.p_00_r;
                s.mr = s.mr + 1;
            case 38 % a(38) = p.alpha*p.ba * p_10_r;
                s.mr = s.mr + 1;
                
            case 39 % a(39) = p.ta * s.ma;
                s.auf = s.auf + 1;
            case 40 % a(40) = p.tr * s.mr;
                s.ruf = s.ruf + 1;
                
            case 41 % a(41) = p.kfa * s.auf;
                s.auf = s.auf - 1;
                s.a   = s.a   + 1;
            case 42 % a(42) = p.kfr * s.ruf;
                s.ruf = s.ruf - 1;
                s.r   = s.r   + 1;
                
            case 43 % a(43) = p.kda  * s.a * (s.a-1); 
                s.a  = s.a  - 2;
                s.a2 = s.a2 + 1;
            case 44 % a(44) = p.k_da * s.a2;
                s.a  = s.a  + 2;
                s.a2 = s.a2 - 1;
            case 45 % a(45) = p.kdr  * s.r * (s.r-1);
                s.r  = s.r  - 2;
                s.r2 = s.r2 + 1;
            case 46 % a(46) = p.k_dr * s.r2;
                s.r  = s.r  + 2;
                s.r2 = s.r2 - 1;
            case 47 % a(47) = p.kt  * s.r2 * (s.r2-1);
                s.r2 = s.r2 - 2;
                s.r4 = s.r4 + 1;
            case 48 % a(48) = p.k_t * s.r4;
                s.r2 = s.r2 + 2;
                s.r4 = s.r4 - 1;
                
            case 49 % a(49) = p.da * s.ma;
                s.ma = s.ma - 1;
            case 50 % a(50) = p.dr * s.mr;
                s.mr = s.mr - 1;
            case 51 % a(51) = p.lambda*f(X) * s.auf;
                s.auf = s.auf - 1;
            case 52 % a(52) = f(X) * s.ruf;
                s.ruf = s.ruf - 1;
            case 53 % a(53) = p.lambda*f(X) * s.a;
                s.a = s.a - 1;
            case 54 % a(54) = f(X) * s.r;
                s.r = s.r - 1;
            case 55 % a(55) = p.lambda*f(X) * s.a2;
                s.a2 = s.a2 - 1;
            case 56 % a(56) = f(X) * s.r2;
                s.r2 = s.r2 - 1;
            case 57 % a(57) = f(X) * s.r4;
                s.r4 = s.r4 - 1;
                
            case 58 % a(58) = f(X) * s.p_10_a;
                s.p_10_a = s.p_10_a - 1;
                s.p_00_a = s.p_00_a + 1;
            case 59 % a(59) = f(X) * s.p_11_a;
                s.p_11_a = s.p_11_a - 1;
                s.p_01_a = s.p_01_a + 1;
            case 60 % a(60) = f(X) * s.p_12_a;
                s.p_12_a = s.p_12_a - 1;
                s.p_02_a = s.p_02_a + 1;
            case 61 % a(61) = f(X) * s.p_01_a;
                s.p_01_a = s.p_01_a - 1;
                s.p_00_a = s.p_00_a + 1;
            case 62 % a(62) = f(X) * s.p_11_a;
                s.p_11_a = s.p_11_a - 1;
                s.p_10_a = s.p_10_a + 1;
            case 63 % a(63) = 2*f(X) * s.p_02_a;
                s.p_02_a = s.p_02_a - 1;
                s.p_01_a = s.p_01_a + 1;
            case 64 % a(64) = 2*f(X) * s.p_12_a;
                s.p_12_a = s.p_12_a - 1;
                s.p_11_a = s.p_11_a + 1;
            case 65 % a(65) = 2*p.epsilon*f(X) * s.p_L2_a;
                s.p_L2_a = s.p_L2_a - 1;
                s.p_L1_a = s.p_L1_a + 1;
            case 66 % a(66) = p.epsilon*f(X) * s.p_L1_a;
                s.p_L1_a = s.p_L1_a - 1;
                s.p_L0_a = s.p_L0_a + 1;
                
            case 67 % a(67) = f(X) * s.p_10_r;
                s.p_10_r = s.p_10_r - 1;
                s.p_00_r = s.p_00_r + 1;
            case 68 % a(68) = f(X) * s.p_11_r;
                s.p_11_r = s.p_11_r - 1;
                s.p_01_r = s.p_01_r + 1;
            case 69 % a(69) = f(X) * s.p_12_r;
                s.p_12_r = s.p_12_r - 1;
                s.p_02_r = s.p_02_r + 1;
            case 70 % a(70) = f(X) * s.p_01_r;
                s.p_01_r = s.p_01_r - 1;
                s.p_00_r = s.p_00_r + 1;
            case 71 % a(71) = f(X) * s.p_11_r;
                s.p_11_r = s.p_11_r - 1;
                s.p_10_r = s.p_10_r + 1;
            case 72 % a(72) = 2*f(X) * s.p_02_r;
                s.p_02_r = s.p_02_r - 1;
                s.p_01_r = s.p_01_r + 1;
            case 73 % a(73) = 2*f(X) * s.p_12_r;
                s.p_12_r = s.p_12_r - 1;
                s.p_11_r = s.p_11_r + 1;
            case 74 % a(74) = 2*p.epsilon*f(X) * s.p_L2_r;
                s.p_L2_r = s.p_L2_r - 1;
                s.p_L1_r = s.p_L1_r + 1;
            case 75 % a(75) = p.epsilon*f(X) * s.p_L1_r;
                s.p_L1_r = s.p_L1_r - 1;
                s.p_L0_r = s.p_L0_r + 1;
                
        end
        
        % Calculate total ssrA tags X
        s.X = getTags(s);

        % Store outputs
        tstore(time_idx) = t;
        Xstore(time_idx,:) = pack_species(s);

        % Update time and counter
        t = t+tau;
        time_idx = time_idx + 1;
        
%         fprintf('%d\t%f\t%f\t%f\n',int32(time_idx),t,s.a2,s.r4)

%         if ~mod(time_idx,2e4)
%             disp(['Progress: ' num2str(t/tend*100) '%']);
%         end
        
    end
    
    % Assign final outputs
    t = tstore(1:time_idx-1);
    s = unpack_species(Xstore(1:time_idx-1,:));
    
    length(s.p_00_a)
end

function state = pick_rxn(probs)
    % Vectorized function to pick a state to jump to based on a list of
    % probabilities. All we're doing is creating a list of bins
    % that a randomly chosen number from 0-1 can fall into, and then
    % searching for that bin using vector operations. Same thing can be
    % accomplished via a for loop.
    %
    % This vectorized version is significantly faster than the double for
    % loop implementation for large reaction lists.
    
    selection = rand;
    
    M = cumsum([0, probs]);
    M(end) = 1;
    
    C = M >= selection;
    D = M < selection;
    
    state = find(~(C(2:end) - D(1:end-1)));
end

function X = getTags(s)
    % assume bound promoters have accessible tags
    % assume unfolded proteins have active ssrA tags
    X = sum([0*s.p_00_a, 4*s.p_01_a, 4*s.p_02_a, 2*s.p_10_a, 6*s.p_11_a, 10*s.p_12_a, ...
        0*s.p_00_r, 4*s.p_01_r, 4*s.p_02_r, 2*s.p_10_r, 6*s.p_11_r, 10*s.p_12_r, ...
        0*s.p_L0_a, 4*s.p_L1_a, 8*s.p_L2_a, ...
        0*s.p_L0_r, 4*s.p_L1_r, 8*s.p_L2_r, ...
        0*s.ma, 0*s.mr, s.auf, s.ruf, s.a, s.r, 2*s.a2, 2*s.r2, 4*s.r4]);
end

%%% These species packing and unpacking functions automate the process of
%%% collecting vectors of species and differential equations. Putting this
%%% code in a separate function cuts down on errors from mixing up
%%% parameter / species order and copying and pasting variable collection
%%% code throughout your program.

function y = pack_species(s)
    y = [s.p_00_a, s.p_01_a, s.p_02_a, s.p_10_a, s.p_11_a, s.p_12_a, ...
        s.p_00_r, s.p_01_r, s.p_02_r, s.p_10_r, s.p_11_r, s.p_12_r, ...
        s.p_L0_a, s.p_L1_a, s.p_L2_a, ...
        s.p_L0_r, s.p_L1_r, s.p_L2_r, ...
        s.ma, s.mr, s.auf, s.ruf, s.a, s.r, s.a2, s.r2, s.r4, ...
        s.X];
end

function s = unpack_species(y)
    s.p_00_a = y(:,1);
    s.p_01_a = y(:,2);
    s.p_02_a = y(:,3);
    s.p_10_a = y(:,4);
    s.p_11_a = y(:,5);
    s.p_12_a = y(:,6);
    s.p_00_r = y(:,7);
    s.p_01_r = y(:,8);
    s.p_02_r = y(:,9);
    s.p_10_r = y(:,10);
    s.p_11_r = y(:,11);
    s.p_12_r = y(:,12);
    s.p_L0_a = y(:,13);
    s.p_L1_a = y(:,14);
    s.p_L2_a = y(:,15);
    s.p_L0_r = y(:,16);
    s.p_L1_r = y(:,17);
    s.p_L2_r = y(:,18);
    s.ma = y(:,19);
    s.mr = y(:,20);
    s.auf = y(:,21);
    s.ruf = y(:,22);
    s.a = y(:,23);
    s.r = y(:,24);
    s.a2 = y(:,25);
    s.r2 = y(:,26);
    s.r4 = y(:,27);
    s.X = y(:,28);
end

