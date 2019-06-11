function Linear_sqBOLD(data_dir,sub_num)
% Linear_sqBOLD.m usage:
%
%       Linear_sqBOLD(data_dir,sub_num)
%
% Solves the log-linear version of the qBOLD model (Yablonskiy & Haacke, 1994)
% for ASE data located in DATA_DIR (string), with subject number SUB_NUM (which
% should be supplied as an integer). Produces and saves out NIFTY files
% containing maps of R2', DBV, and OEF.
%
% This script is designed to be used with the following dataset:
%       Cherukara MT, Stone AJ, Chappell MA, Blockley NP. Data acquired to
%       demonstrate model-based Bayesian inference of brain oxygenation using
%       quantitative BOLD, Oxford University Research Archive 2018. doi: <Please
%       see ORA entry for DOI> 
%
% It is adapted from ase_qbold_3d.m (Stone & Blockley, 2017), and requires
% read_avw.m and save_avw.m
%
% 
%       Copyright (C) University of Oxford, 2016-2018
%
% 
% Created by AJ Stone, 30 June 2016
%
% CHANGELOG:
%
% 2018-08-27 (MT Cherukara). Adapted for use with the model-based qBOLD
%       investigation, for comparison with a Bayesian analysis.


% Hardcoded scan parameters (in s)
start_tau = -0.028;
delta_tau =  0.004;
end_tau   =  0.064;

% Model constants
Hct = 0.40;             % hct ratio in small vessels
dChi0 = 0.264*10^-6;    % ppm, sus difference between fully oxy & deoxy rbc's
gamma = 2.675*10^4;     % rads/(secs.Gauss) gyromagnetic ratio
B0 = 3*10^4;            % Gauss, Field strength


% Load ASE data
[V, ~, scales] = read_avw([data_dir,'/sub-0',num2str(sub_num),'/func/sub-0',num2str(sub_num),'_acq-flair_ase.nii.gz']);
[x, y, z, v] = size(V);


% Fit R2'

% X
tau = start_tau:delta_tau:end_tau;

% Check that tau matches the number of volumes 
if (length(tau) ~= v)
    disp('List of Tau values doesn''t match the number of volumes') 
    sprintf('Number of volumes = %1.0f', v)
    
else

    % Y
    ln_Sase = log(V); 
    ln_Sase(isnan(ln_Sase)) = 0; 
    ln_Sase(isinf(ln_Sase)) = 0;

    Tc = 0.015; % cutoff time for monoexponential regime [s]
    tau_lineID = find(tau > Tc); % tau's to be used for R2' fitting  
    w = 1./tau(1,tau_lineID)'; % weightings for lscov
    p = zeros(x,y,z,2); 

    % Loop through all voxels, solving the linear model for each one
    for xID = 1:x
        for yID = 1:y
            for zID = 1:z
                
                % LSCOV: fit linear regime
                X = [ones(length(tau(1,tau_lineID)'),1) tau(1,tau_lineID)'];
                Y = squeeze(ln_Sase(xID,yID,zID,tau_lineID));

                p(xID,yID,zID,:) = flipud(lscov(X,Y,w));        
            end
        end
    end

    % Calculate Physiological Parameters
    s0_id = find(tau == 0);
    r2p = -p(:,:,:,1); 
    c = p(:,:,:,2);
    dbv = c - ln_Sase(:,:,:,s0_id);
    oef = r2p./(dbv.*gamma.*(4./3).*pi.*dChi0.*Hct.*B0);


    % Output parameter niftis

    save_avw(r2p, strcat('sub-0',num2str(sub_num),'_qbold-linear_param_R2p'), 'f', scales);
    save_avw(dbv, strcat('sub-0',num2str(sub_num),'_qbold-linear_param_DBV'), 'f', scales);
    save_avw(oef, strcat('sub-0',num2str(sub_num),'_qbold-linear_param_OEF'), 'f', scales);

end % if (length(tau) ~= v) ... else ...

end % function