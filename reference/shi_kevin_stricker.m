function shi_kevin_stricker
%%% Calls the external file stochastic_fun.m OR stochastic_fun_mex.mexw64
%%% Calls the external file stochastic_mod.m OR stochastic_mod_mex.mexw64


%%% Data generation and figure plotting are separated into different
%%% functions. When a particular function is finished running, it saves its
%%% results to a .mat file so you can quickly access the data again without
%%% rerunning the model.

    clc; close all

    %%% Generate data %%%
    
    run_fig4b();
    run_fig4d();
    run_figs17();
    run_figs19();
    run_figm0();
    run_figm1();
    run_figm2();
    run_figm3();

    %%% Create figures %%%
     
    plot_fig4b();
    plot_fig4d();
    plot_figs17();
    plot_figs19();

    %%% Modification
    
    fig4b_analysis();
    plot_figm0();
    plot_figm1();
    plot_figm2();
    plot_figm3();
    
end

%%
function run_fig4b(seed)
    disp('Generating Figure 4b data...');
    
    % These routines take in and store an optional random seed for the Gillespie
    % simulations. If this function is called again with the same seed, it
    % will produce identical data each time it is run, which can be useful
    % for debugging and testing your implementation.
    
    if ~exist('seed','var')
        % Generate a random seed based on the current time. This seed will
        % be saved to fig1_data.mat
        seed = now;
        rng(seed);
    else
        rng(seed);
    end
    
    IPTG = 2;
    ara = 0.7;
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_fun_mex(200,p,s0);
    toc
    
    save('fig4b_data.mat');
end

function plot_fig4b()
    % Load the previously generated data from the .mat file
    load('fig4b_data.mat');
    
    % Plot the figure
    figure
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
        data.t_stoch,data.s_stoch.r4./1000,'r', ...
        data.t_stoch,data.s_stoch.mr./1000,'k')
    xlabel('Time (min)')
    ylabel('Molecules (x1000)')
    legend('AraC dimers','LacI tetramers', 'lacI mRNA','Location','Best')
    title('Fig 4B, 0.7% arabinose, 2 mM IPTG')
end

%%
function run_fig4d(seed)
    disp('Generating Figure 4d data...');
    
    if ~exist('seed','var')
        seed = now;
        rng(seed);
    else
        rng(seed);
    end
    
    ara = 0.7; % keep this constant
    
    periods = {};
%     IPTGconcs = [2, 10, 25]; % mM
    IPTGconcs = [0, 0.5, 1, 2, 5, 8, 10, 15, 20, 25]; % mM
    for i = 1:length(IPTGconcs)
        IPTG = IPTGconcs(i);
        [p,s0] = get_defaults(ara,IPTG);
        
        iperiods = [];
        
        for j = 1:10 % multiple runs/IPTGconc
            tic
            data.IPTG = IPTG;
            [data.t,data.s] = stochastic_fun_mex(200,p,s0);
            toc
            
            % Just do the data processing here - the mex code runs fast enough
            
            % Coarsen data
            t = data.t(1:1000:end);
            r4 = data.s.r4(1:1000:end);
            
            r4 = smooth(r4);
            [pks,locs] = findpeaks(r4,'minpeakheight',150,'minpeakdistance',35); % good compromise thresholds
            
            % Manually check peak detection
%             figure
%             plot(t,r4,'r')
%             hold on
%             plot(t(locs),pks+10,'k^','markerfacecolor',[1 0 0]);
            
            % Get periods
            tp = t(locs);
            per = tp(2:end) - tp(1:end-1);
            iperiods = vertcat(iperiods,per);
            
            
        end
        periods{i} = iperiods;
    end
    clear data
    save('fig4d_data.mat');
end

function plot_fig4d()
    % Load the previously generated data from the .mat file
    load('fig4d_data.mat');
    
    % get 1 sd norm dist for each IPTG conc
    figure
    hold on
    for i = 1:length(IPTGconcs)
        per = cell2mat(periods(i));
        mu = mean(per);
        sd = std(per);
        
        plot(IPTGconcs(i),mu,'o')
        errorbar(IPTGconcs(i),mu,sd)
        xlim([-2 30])
        ylim([10 50])
        
    end
    xlabel('IPTG (mM)')
    ylabel('Period (min)')
    title('Fig 4D, IPTG Period Tuning, 0.7% arabinose')
    hold off

end

%%
function run_figs17(seed)
    disp('Generating Figure S17 data...');
    
    if ~exist('seed','var')
        seed = now;
        rng(seed);
    else
        rng(seed);
    end
    
    IPTG = 2;
    ara = 3.5; % high
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_fun_mex(1000,p,s0);
    toc
    
    save('figs17_data.mat');
end

function plot_figs17()
    load('figs17_data.mat');
    
    figure
    plot(data.t_stoch,data.s_stoch.a2,'g', ...
        data.t_stoch,data.s_stoch.r4,'r', ...
        data.t_stoch,data.s_stoch.mr,'k')
    xlabel('Time (min)')
    ylabel('Molecules')
    legend('AraC dimers','LacI tetramers', 'lacI mRNA','Location','Best')
    title('Fig S17, Bistable Oscillations, 3.5% arabinose, 2 mM IPTG')
end

%% negative feedback only
function run_figs19(seed)
    disp('Generating Figure S19 data...');
    
    if ~exist('seed','var')
        seed = now;
        rng(seed);
    else
        rng(seed);
    end
    
    IPTG = 2;
    ara = 0; % no positive feedback
    [p,s0] = get_defaults(ara,IPTG);
    
    % kl = 0
    p.kl = 0; % 1/min
    p.br = 10; % 1/min
    s0.p_00_a = 0; % no positive feedback
    s0.r4 = 500; % start with some LacI tetramers as in figure, prevents initial spike
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_fun_mex(600,p,s0);
    toc
    
    save('figs19A_data.mat');
    
    % kl = 0.36
    p.kl = 0.36; % 1/min
    p.br = 10; % 1/min
    s0.p_00_a = 0; % no positive feedback
    s0.r4 = 500; % start with some LacI tetramers as in figure, prevents initial spike
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_fun_mex(600,p,s0);
    toc
    
    save('figs19B_data.mat');
    
    % kl = 360
    p.kl = 360; % 1/min
    p.br = 10; % 1/min
    s0.p_00_a = 0; % no positive feedback
    s0.r4 = 500; % start with some LacI tetramers as in figure, prevents initial spike
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_fun_mex(600,p,s0);
    toc
    
    save('figs19C_data.mat');
end

function plot_figs19()

    figure

    load('figs19A_data.mat');
    subplot(3,1,1)
    plot(data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    legend('A','Location','NorthWest')
    title('Fig S19, Negative Feedback Only, Free LacI tetramer')
    
    load('figs19B_data.mat');
    subplot(3,1,2)
    plot(data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    legend('B','Location','NorthWest')
    
    load('figs19C_data.mat');
    subplot(3,1,3)
    plot(data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    legend('C','Location','NorthWest')
    xlabel('Time (min)')
end

%% Modification

function fig4b_analysis()
    load('fig4b_data.mat');
    
    % also plot f(X) here?
    fX = p.gamma./(p.ce + data.s_stoch.X);
    
    figure
    subplot(2,1,1)
    plot(data.t_stoch,data.s_stoch.X./1000,...
        data.t_stoch,data.s_stoch.a2./1000)
    ylabel('Molecules (x1000)')
    legend('Total ssrA tags','AraC dimers','Location','Best')
    title('Degradation Rates, 0.7% arabinose, 2 mM IPTG')
    
    subplot(2,1,2)
    semilogy(data.t_stoch,fX./1000)
    xlabel('Time (min)')
    ylabel('f(X) (min^{-1})')
end

function run_figm0()
    disp('Generating Figure M0 data...');
    
    if ~exist('seed','var')
        seed = now;
        rng(seed);
    else
        rng(seed);
    end
    
    IPTG = 2;
    ara = 0.7;
    [p,s0] = get_defaults(ara,IPTG);
    s0.a2 = 100;
    s0.r4 = 100;
    
    % Make ssrA degradation f(X) = const.
    p.gamma = 1e7;
    p.ce = 1e9;
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_fun_mex(200,p,s0);
    toc
    
    save('figm0_data1.mat');
    
    p.gamma = 1e7;
    p.ce = 1e8;
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_fun_mex(200,p,s0);
    toc
    
    save('figm0_data2.mat');
    
    p.gamma = 1e7;
    p.ce = 1e7;
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_fun_mex(200,p,s0);
    toc
    
    save('figm0_data3.mat');
    
    
end

function plot_figm0()
    figure
    
    load('figm0_data1.mat');
    subplot(3,1,1)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
    data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('Fig M1, Constant Degradation, k_d = 0.01')
    
    load('figm0_data2.mat');
    subplot(3,1,2)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
        data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('k_d = 0.1')
    
    load('figm0_data3.mat');
    subplot(3,1,3)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
        data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('k_d = 1')
    xlabel('Time (min)')
end

function run_figm1()
    disp('Generating Figure M1 data...');
    
    if ~exist('seed','var')
        seed = now;
        rng(seed);
    else
        rng(seed);
    end
    
    IPTG = 2;
    ara = 0.7;
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_mod_mex(200,p,s0);
    toc
    
    save('figm1_data.mat');
end

function plot_figm1()
    load('figm1_data.mat');
    
    figure
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
        data.t_stoch,data.s_stoch.r4./1000,'r', ...
        data.t_stoch,data.s_stoch.mr./1000,'k')
    xlabel('Time (min)')
    ylabel('Molecules (x1000)')
    legend('AraC dimers','LacI tetramers', 'lacI mRNA','Location','Best')
    title('Fig M1, Ninfa oscillator,  0.7% arabinose, 2 mM IPTG')
end

function run_figm2()
    disp('Generating Figure M2 data...');
    
    if ~exist('seed','var')
        seed = now;
        rng(seed);
    else
        rng(seed);
    end
    
    IPTG = 0;
    ara = 0.7;
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_mod_mex(200,p,s0);
    toc
    
    save('figm2_data1.mat');
    
    IPTG = 10;
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_mod_mex(200,p,s0);
    toc
    
    save('figm2_data2.mat');
    
    IPTG = 25;
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_mod_mex(200,p,s0);
    toc
    
    save('figm2_data3.mat');
end

function plot_figm2()
    figure
    
    load('figm2_data1.mat');
    subplot(3,1,1)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
    data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('Fig M2, Ninfa oscilaltor, IPTG = 0 mM')
    
    load('figm2_data2.mat');
    subplot(3,1,2)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
        data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('IPTG = 10 mM')
    
    load('figm2_data3.mat');
    subplot(3,1,3)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
        data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('IPTG = 25 mM')
    xlabel('Time (min)')
end

function run_figm3()
    disp('Generating Figure M2 data...');
    
    if ~exist('seed','var')
        seed = now;
        rng(seed);
    else
        rng(seed);
    end
    
    IPTG = 2;
    ara = 0.0;
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_mod_mex(200,p,s0);
    toc
    
    save('figm3_data1.mat');
    
    ara = 2.5;
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_mod_mex(200,p,s0);
    toc
    
    save('figm3_data2.mat');
    
    ara = 3.5;
    [p,s0] = get_defaults(ara,IPTG);
    
    tic
    [data.t_stoch,data.s_stoch] = stochastic_mod_mex(200,p,s0);
    toc
    
    save('figm3_data3.mat');
end

function plot_figm3()
    figure
    
    load('figm3_data1.mat');
    subplot(3,1,1)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
    data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('Fig M2, Ninfa oscilaltor, arabinose = 0%')
    
    load('figm3_data2.mat');
    subplot(3,1,2)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
        data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('arabinose = 2.5%')
    
    load('figm3_data3.mat');
    subplot(3,1,3)
    plot(data.t_stoch,data.s_stoch.a2./1000,'g', ...
        data.t_stoch,data.s_stoch.r4,'r')
    ylabel('Molecules')
    title('arabinose = 3.5%')
    xlabel('Time (min)')
end

%% Utility functions

function [p,s] = get_defaults(ara,IPTG)
    
    % IPTG dependence (IPTG in mM)
    k_r = 1.8; % 1/min
    crmin = 0.01; % 1/molecules
    crmax = 0.2; % 1/molecules
    kr1 = 0.035; % mM
    b1 = 2;
    kr = k_r*((crmax-crmin)/(1+(IPTG/kr1)^b1)+crmin);
    
    % arabinose dependence (in %w/v)
    k_a = 1.8; % 1/min
    camin = 0; % 1/molecules
    camax = 1; % 1/molecules
    ka1 = 2.5; % %
    kr2 = 1.8; % mM
    c1 = 2;
    b2 = 2;
    ka = k_a*((camax-camin)*ara^c1/(ka1^c1+ara^c1)/(1+(IPTG/kr2)^b2)+camin);
    
    %%% Kinetic parameters and model constants
    p.ka = ka;
    p.k_a = k_a;
    p.kr = kr;
    p.k_r = k_r;
    p.kl = 0.36; % 1/min
    p.kul = 0.18; % 1/min
    p.ba = 0.36; % 1/min
    p.br = 0.36; % 1/min
    p.alpha = 20;
    p.ta = 90; % 1/min
    p.tr = 90; % 1/min
    p.kfa = 0.9; % 1/min
    p.kfr = 0.9; % 1/min
    p.kda = 0.018; % 1/min
    p.k_da = 0.00018; % 1/min
    p.kdr = 0.018; % 1/min
    p.k_dr = 0.00018; % 1/min
    p.kt = 0.018; % 1/min
    p.k_t = 0.00018; % 1/min
    p.da = 0.54; % 1/min
    p.dr = 0.54; % 1/min
    p.lambda = 2.5;
    p.epsilon = 0.2;
    
    p.gamma = 1080; % molecules/min
    p.ce = 0.1; % molecules
    
    %%% Initial species concentrations
    
    Na = 50; % molecules, activator plasmids
    Nr = 25; % molecules, repressor plasmids
    
    s.p_00_a = Na; % unlooped activator promoter
    s.p_01_a = 0;
    s.p_02_a = 0;
    s.p_10_a = 0;
    s.p_11_a = 0;
    s.p_12_a = 0;
    
    s.p_00_r = Nr; % unlooped repressor promoter
    s.p_01_r = 0;
    s.p_02_r = 0;
    s.p_10_r = 0;
    s.p_11_r = 0;
    s.p_12_r = 0;
    
    s.p_L0_a = 0; % looped activator promoter
    s.p_L1_a = 0;
    s.p_L2_a = 0;
    
    s.p_L0_r = 0; % looped repressor promoter
    s.p_L1_r = 0;
    s.p_L2_r = 0;
    
    s.ma = 0; % mRNAs
    s.mr = 0;
    
    s.auf = 0; % unfolded proteins
    s.ruf = 0;
    s.a = 0; % proteins
    s.r = 0;
    s.a2 = 0; % multimers
    s.r2 = 0;
    s.r4 = 0;
    
    s.X = getTags(s);
    
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

