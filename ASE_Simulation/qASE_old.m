% qASE.m
%
% Quantitative Asymmetric Spin Echo Sequence Simulation. Generates data to
% be used in Bayesian inference on parameters (Asymmetric_Bayes.m)
%
% Requires qASE_model.m, which should be in the same folder.
%
% 
%       Copyright (C) University of Oxford, 2016-2018
%
% 
% Created by MT Cherukara, 17 May 2016
%
% CHANGELOG:
%
% 2018-08-23. Various changes.

clear; 

plot_fig = 1;       % set to 1 in order to plot a figure
save_data = 1;      % set to 1 in order to save out ASE data


%% Model Parameters

% constants 
params.B0   = 3.0;          % T         - static magnetic field
params.dChi = 2.64e-7;      % parts     - susceptibility difference
params.gam  = 2.67513e8;    % rad/s/T   - gyromagnetic ratio

% scan parameters 
params.TE   = 0.082;        % s         - echo time
params.TR   = 3.000;        % s         - repetition time
params.TI   = 0;            % s         - FLAIR inversion time

% model fitting parameters
params.S0   = 100;          % a. units  - signal
params.R2t  = 1/0.087;      % 1/s       - rate constant, tissue
params.R2e  = 4;            % 1/s       - rate constant, extracellular
params.dF   = 5;            % Hz        - frequency shift
params.lam0 = 0.0;          % no units  - ISF/CSF signal contribution
params.zeta = 0.03;         % no units  - deoxygenated blood volume
params.OEF  = 0.40;         % no units  - oxygen extraction fraction
params.Hct  = 0.400;        % no units  - fractional hematocrit
params.T1t  = 1.200;        % s         - tissue T1
params.T1b  = 1.580;        % s         - blood T1
params.T1e  = 3.870;        % s         - CSF T1

% analysis parameters
params.tc_man = 0;          % BOOL      - should Tc be defined manually?
params.tc_val = 0.0;        % s         - manual Tc (if tc_man = 1)

% noise
params.SNR = 100;


%% Compute Model

% define tau values that we want to simulate
tau = (-28:4:64)/1000; % for testing
% tau = linspace(-0.028,0.064,1000); % for visualising

np = length(tau);

% call qASE_model
[S_total,params] = qASE_model(tau,params.TE,params);



%% Add Noise
S_sample = S_total + (S_total.*randn(1,np)./params.SNR);
S_sample(S_sample < 0) = 0;
S_sample = S_sample./max(S_total);

% calculate maximum data standard deviaton
params.sig = min(S_sample)/params.SNR;


%% Plot Figure
if plot_fig
    
    % create a figure
    figure(1); hold on; box on;
    
    % plot the signal
    S_log = log((S_total)./max(S_total));
    l.s = plot(1000*tau,S_log,'-');
    xlim([(1000*min(tau))-4, (1000*max(tau))+4]);
    
    % labels on axes
    xlabel('Spin Echo Displacement \tau (ms)');
    ylabel('Log (Signal)');

    
end % if plot_fig


%% Save Data
if save_data
    dat_title = strcat('ASE_signal_',date);
    
    % pull out values of TE and tau
    T_sample = tau;
    
    if length(params.TE) ~= length(tau)
        TE_sample(1:length(tau)) = params.TE;
    else
        TE_sample = params.TE;
    end
    
    % Save the data out
    save(dat_title,'T_sample','S_sample','TE_sample','params');
end % if save_data
