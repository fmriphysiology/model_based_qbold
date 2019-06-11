function params = genParams(varargin)
% Generate a parameter structure PARAMS with the standard set of values to be 
% used in MTC_qASE.m, and derived qBOLD model scripts. Allows the user to
% specify the following parameters as name-value pair arguments:
%
%       OEF, DBV, vCSF, TE, SNR, Model, incIV, incT2
%
% Other parameters can be changed manually after the PARAMS structure has been
% created
%
% MT Cherukara
% 2018-10-11
%
% CHANGELOG:
%
% 2019-06-11 (MTC). Updated for re-submitted model-fitting paper.


% Read user-input parameters
q = inputParser;
addParameter(q, 'TE'    , 0.082 , @isnumeric);      % Echo time
addParameter(q, 'OEF'   , 0.40  , @isnumeric);      % Oxygen extraction fraction
addParameter(q, 'DBV'   , 0.03  , @isnumeric);      % Deoxygenated blood volume
addParameter(q, 'vCSF'  , 0.00  , @isnumeric);      % CSF volume
addParameter(q, 'SNR'   , inf   , @isnumeric);      % Signal to noise ratio
addParameter(q, 'Model' , 'Full' );                 % Simulated qBOLD model
addParameter(q, 'incIV' , true  , @islogical);      % Include Blood compartment
addParameter(q, 'incT2' , true  , @islogical);      % Include T2 weightings

parse(q,varargin{:});
r = q.Results;


% Physical constants
params.dChi = 0.264e-6;     % parts     - susceptibility difference
params.gam  = 2.67513e8;    % rad/s/T   - gyromagnetic ratio
params.B0   = 3.0;          % T         - static magnetic field

% Sequence parameters
params.TE   = r.TE;         % s         - echo time
params.TR   = 3.000;        % s         - repetition time
params.TI   = 1.210;        % s         - FLAIR inversion time

% Relaxometry parameters
params.T1t  = 1.200;        % s         - tissue T1
params.T1b  = 1.580;        % s         - blood T1
params.T1e  = 3.870;        % s         - CSF T1

params.R2t  = 11.5;         % 1/s       - rate constant, tissue
params.R2e  = 4;            % 1/s       - rate constant, extracellular

% Physiological parameters
params.OEF  = r.OEF;        % no units  - oxygen extraction fraction
params.zeta = r.DBV;        % no units  - deoxygenated blood volume
params.Hct  = 0.400;        % no units  - fractional hematocrit
params.S0   = 100;          % a. units  - signal 


% Simulation Parameters
params.model  = r.Model;    % STRING    - model type: 'Full','Asymp','Phenom'
params.SNR    = r.SNR;      % no units  - simulated signal to noise ratio
params.contr  = 'OEF';      % STRING    - contrast source: 'OEF','R2p','dHb',...
params.tc_man = 0;          % BOOL      - should Tc be defined manually?
params.tc_val = 0.0;        % s         - manual Tc (if tc_man = 1)
params.incT1  = 0;          % BOOL      - should T1 differences be considered?
params.incIV  = r.incIV;    % BOOL      - should the blood compartment be added?
params.incT2  = r.incT2;    % BOOL      - should blood compartment be included?

% Derived parameters
params.dw = (4/3)*pi*params.gam*params.B0*params.dChi*params.Hct*params.OEF;