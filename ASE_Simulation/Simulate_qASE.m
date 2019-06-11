% Simulate_qASE.m
%
% Quantitative Asymmetric Spin Echo Sequence Simulation. Generates data to
% be used in Bayesian inference on parameters (gridSearchBayesian.m and
% others). Based on MTC_qBOLD.m. Requires genParams.m and qASE_model.m, which
% must be in the PATH.
%
% 
%       Copyright (C) University of Oxford, 2016-2019
%
% 
% Created by MT Cherukara, 17 May 2016
%
% CHANGELOG:
%
% 2019-06-11 (MTC). Updated for re-submitted model-fitting paper.
%
% 2018-10-10 (MTC). Added the option to specify which model to use, and whether
%       to include T1 effects (including FLAIR), as PARAMS. Added genParams.m as
%       a function to generate the PARAMS structure, rather than hard-coding it
%       in this script
%
% 2018-04-05 (MTC). Added the option to set the critical time (TC) manually as a
%       pair of parameters. Also changed MTC_ASE_tissue.m accordingly.
%
% 2018-01-12 (MTC). Change the way the noise standard deviation params.sig
%       is calculated and applied, so that it actually results in the SNR
%       that we want. Also some cleanup.
%
% 2017-08-07 (MTC). Added a call to MTC_qASE_model, rather than repeating
%       its contents here, so that the whole thing is more modular. Removed
%       the automatic saving of the plot, as it's kind of unnecessary
%
% 2017-04-04 (MTC). Modified the plotting parameters, removed a bunch of
%       unnecessary variables from the params struct.
%
% 2016-11-08 (MTC). Changed tau range to go up to 64 ms, changed some
%       other constants, added saving out of data sampled at key points,
%       with noise added in too.
% 
% 2016-05-27 (MTC). Changed the way the compartments are added together,
%       which was apparently causing problems (I don't know why).

clear; 
% close all;

plot_fig = 1;       
save_data = 1;      % set to 1 in order to save out ASE data


%% Model Parameters

% Create a parameter structure
params = genParams;

% Physiology
params.zeta = 0.03;         % no units  - deoxygenated blood volume
params.OEF  = 0.40;         % no units  - oxygen extraction fraction

% Simulation
params.model  = 'Full';     % STRING    - model type: 'Full','Asymp',
params.contr  = 'OEF';      % STRING    - contrast source: 'OEF','R2p'
params.incT1  = 1;          % BOOL      - should T1 differences be considered?
params.incT2  = 1;          % BOOL      - should T2 differences be considered?
params.incIV  = 1;          % BOOL      - should blood compartment be included?

% noise
params.SNR = 100;


%% Compute Model

% define tau values that we want to simulate
% tau = (-16:8:64)/1000; % for testing
tau = linspace(-0.028,0.064,1000); % for visualising

np = length(tau);

% call MTC_qASE_model
[S_total,params] = qASE_model(tau,params.TE,params);

% Normalize to the spin-echo
SEind = find(tau > -1e-9,1);
S_total = S_total./S_total(SEind);


%% Optionally add noise 
S_sample = S_total + (randn(1,np)./params.SNR);
S_sample(S_sample < 0) = 0;

% calculate maximum data standard deviaton
params.sig = min(S_sample)/params.SNR;


%% Plot Figure
if plot_fig
    
    % create a figure
    figure(1); hold on; box on; 
    
    % plot the signal
    S_log = log(S_total);
    l.s = plot(1000*tau,S_log,'-');
    xlim([(1000*min(tau)), (1000*max(tau))]);
    
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
    save(dat_title,'T_sample','S_sample','S_total','TE_sample','params');
end % if save_data
