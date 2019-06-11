function [S,PARAMS] = qASE_model(TAU,TE,PARAMS)
    % Calculates ASE signal from a single voxel. Usage:
    %
    %       [S,PARAMS] = qASE_model(TAU,TE,PARAMS)
    %
    % Input:  TAU    - A vector of tau values (refocussing pulse offsets)
    %                  in seconds.
    %         TE     - A single value, or a vector of equal length to TAU,
    %                  of echo time values in seconds.
    %         PARAMS - The structure containing all the parameters used in
    %                  the simulation and model inference. This can be
    %                  optionally returned as an output, with extra values
    %                  (such as dw) which are calculated within, returned
    %                  as the output.
    %
    % Output: S      - A vector of signal values of same size as T.
    %         PARAMS - [OPTIONAL] The modified parameter structure.
    %
    % For use in qASE.m, gridSearchBayesian.m, etc. uses the method described
    % by Yablonskiy & Haacke (1994), and He & Yablonskiy (2007), with additional
    % components as described by Berman & Pike (2017) and Simon et al. (2016)
    %
    % 
    %       Copyright (C) University of Oxford, 2016-2019
    %
    % 
    % Created by MT Cherukara, 16 May 2016
    %
    % CHANGELOG:
    %
    % 2019-05-27 (MTC). Added alternative models for calculating intravascular
    %       signal, especially the 'powder model' (Sukstanskii, 2001), and the
    %       'linear' model (Simon et al., 2016). Removed other compartments that
    %       aren't included in the paper.
    %
    % 2018-10-30 (MTC). Added in a parameter for short-tau DBV offset (referred
    %       to in my notes as beta), and for linear R2' scaling factor SR, and
    %       for [dHb] exponent beta (a.k.a. gamma from Ogawa)
    %
    % 2018-10-23 (MTC). Re-wrote the Asymptotic model such that it depends on
    %       R2', rather than dw, in its actual calculations. This makes no
    %       actual difference to the calculation, but makes it more intiutive
    %       when seeing what R2' actually does. This will also allow for
    %       arbitrary scaling of R2', as part of adapting the model to suit
    %       simulated data.
    %
    % 2018-10-10 (MTC). Added the option to specify which model to use in the
    %       PARAMS structure, so that one doesn't have to keep coming into this
    %       code and commenting section in and out. Also added an option to
    %       include the T1 (steady-state-magnetization) contribution, and 
    %       whether to include T2 effects.
    %
    % 2018-10-01 (MTC). Added the option for calculating dw as a function of dHb
    %       concentration, instead of using OEF and Hct. Requires a scaling
    %       parameter kappa (0.03), and an assumtpion of [Hb] and/or Hct.
    %
    % 2018-09-13 (MTC). Condensing the various compartment functions into the
    %       one script, and bringing nomenclature into line with fmriphysiology
    %       github.
    %
    % 2018-05-16 (MTC). Added the T1-based compartment weightings
    %
    % 2018-01-12 (MTC). R2_blood and R2_blood_star were the wrong way
    %       around. Now we have fixed that.
    %
    % 2017-10-10 (MTC). Added the option to supply a vector of TE values
    %       outside the PARAMS struct, and made changes to the called model
    %       functions MTC_ASE_tissue.m, MTC_ASE_extra.m, MTC_ASE_blood.m,
    %       in order to allow for inference on data with any arbitrary set
    %       of TE-TAU pairs.
    %
    % 2017-08-07 (MTC). Added a control in the form of PARAMS.noDW to  
    %       change the way things are done when we're trying to infer on 
    %       R2'. Added PARAMS as an optional output. 
    %       
    % 2017-07-10 (MTC). Commented out the PARAMS.dw recalculation in order
    %       to use this for the R2' and zeta estimation (alongside
    %       MTC_qASE_model_long.m). This should be put back to normal if
    %       using this for any other grid-search (e.g. the classic OEF-DBV)
    %
    % 2016-05-27 (MTC). Updated the way the compartments were summed up.
   
    
%% Check Parameters
% Check whether recently added parameters (e.g. incIV and incT2) are specified
% or not, if not, add them in 
if ~isfield(PARAMS,'incIV')
    PARAMS.incIV = 1;   % included by default
end

if ~isfield(PARAMS,'incT2')
    PARAMS.incT2 = 1;   % included by default
end

if ~isfield(PARAMS,'incT1')
    PARAMS.incT1 = 1;   % included by default
end
    
%% Recalculate dw (or OEF) depending on what is being used to generate the contrast

ctr = lower(PARAMS.contr);

switch ctr
    
    case 'oef'
        PARAMS.dw   = (4/3)*pi*PARAMS.gam*PARAMS.B0*PARAMS.dChi*PARAMS.Hct*PARAMS.OEF;
        PARAMS.R2p  = PARAMS.dw .* PARAMS.zeta;
        
    case 'r2p'
        PARAMS.dw = PARAMS.R2p ./ PARAMS.zeta;
        PARAMS.OEF = PARAMS.dw ./ ( (4/3)*pi*PARAMS.gam*PARAMS.dChi*PARAMS.Hct*PARAMS.B0 );
        
    otherwise
        warning('No contrast source specified, using OEF');
        PARAMS.dw   = (4/3)*pi*PARAMS.gam*PARAMS.B0*PARAMS.dChi*PARAMS.Hct*PARAMS.OEF; 
        PARAMS.R2p  = PARAMS.dw .* PARAMS.zeta;
    
end % switch ctr


%% Calculate important parameters

% relaxation rate constant of blood
PARAMS.R2b  =  4.5 + 16.4*PARAMS.Hct + (165.2*PARAMS.Hct + 55.7)*PARAMS.OEF^2;
PARAMS.R2bp = 10.2 -  1.5*PARAMS.Hct + (136.9*PARAMS.Hct - 13.9)*PARAMS.OEF^2;

% compartment weightings
if PARAMS.incT1
    
    % If we're including T1, we definitely need to include T2 as well
    PARAMS.incT2 = 1;
    
    % spin densities
    nb = 0.775;    
    
    % calculate compartment steady-state magnetization
    m_bld = calcMagnetization(TAU,TE,PARAMS.TR,PARAMS.T1b,1./PARAMS.R2b,PARAMS.TI);
    
    % pull out parameters
    Vb = PARAMS.zeta;
    
    % calculate compartment weightings
    w_bld = m_bld.*nb.*Vb;
    w_tis = 1 - w_bld;
    
else
    
    % Extract compartment weightings
    w_bld = PARAMS.zeta;
    w_tis = 1 - w_bld;
    
end % if PARAMS.incT1 ... else ...


%% Calculate the model

mdl = lower(PARAMS.model);

% Tissue Compartment - based on PARAMS.model choice
switch mdl
    
    case 'full'
        S_tis = calcTissueCompartment(TAU,TE,PARAMS);
    case 'asymp'
        S_tis = calcTissueAsymp(TAU,TE,PARAMS);
    otherwise
        warning('No model specified, using Asymptotic qBOLD');
        S_tis = calcTissueAsymp(TAU,TE,PARAMS);
        
end % switch mdl

% Other compartments
if PARAMS.incIV
    
    % Include and weight the blood compartments
    S_bld = w_bld.*calcBloodCompartment(TAU,TE,PARAMS);
%     S_bld = w_bld.*calcBloodPowder(TAU,TE,PARAMS);

    
    % add it all together:
    S = PARAMS.S0.*( (w_tis.*S_tis) + S_bld);
else
    
    S = PARAMS.S0.*S_tis;
    
end % if PARAMS.incIV
      


end % MAIN


%% calcMagnetization function
function M = calcMagnetization(tau,TE,TR,T1,T2,TI)
    % Calculate steady-state magnetization in an ASE sequence

    % put it all together
    M = 1 - ( 2 - exp(-(TR-TI)./T1)) .* exp(-TI/T1);

end


%% calcTissueCompartment function
function ST = calcTissueCompartment(TAU,TE,PARAMS)
    % Calculate tissue compartment ASE qBOLD signal using the static dephasing
    % model (Yablonskiy & Haacke, 1994)

    % pull out constants
    dw   = PARAMS.dw;
    zeta = PARAMS.zeta;

    % check whether one TE, or a vector, is supplied
    if length(TE) ~= length(TAU)
        TE(2:length(TAU)) = TE(1);
    end

    t0 = abs(TAU.*dw);              % predefine tau
    fint = zeros(1,length(TAU));    % pre-allocate

    for ii = 1:length(TAU)

        % integrate
        fnc0 = @(u) (2+u).*sqrt(1-u).*(1-besselj(0,1.5*t0(ii).*u))./(u.^2);
        fint(ii) = integral(fnc0,0,1);

    end

    ST = exp(-zeta.*fint./3);
    
    % add T2 effect
    if PARAMS.incT2
        ST = ST .* exp(-TE.*PARAMS.R2t);
    end
    
end % function ST = calcTissueCompartment(TAU,TE,PARAMS)


%% calcBloodCompartment function
function SB = calcBloodCompartment(TAU,TE,PARAMS)
    % Calculate intravascular contribution to ASE qBOLD signal using the
    % motional narrowing model (Berman & Pike, 2017)

    % pull out constants
    gam = PARAMS.gam;
    Hct = PARAMS.Hct;
    OEF = PARAMS.OEF;
    B0  = PARAMS.B0;
    Dx0 = PARAMS.dChi;

    % assign constants
    td = 0.0045067;     % diffusion time

    % calculate parameters
    dChi = Dx0.*OEF + 0.14e-6;
    G0   = (4/45)*Hct*(1-Hct)*((dChi*B0)^2);
    kk   = 0.5*(gam^2)*G0*(td^2);


    % check whether one TE, or a vector, is supplied
    if length(TE) ~= length(TAU)
        TE(2:length(TAU)) = TE(1);
    end

    % calculate model
    SB  = exp(-kk.*( (TE./td) + sqrt(0.25 + (TE./td)) + 1.5 - ...
                    (2.*sqrt( 0.25 + ( ((TE + TAU).^2) ./ td ) ) ) - ...
                    (2.*sqrt( 0.25 + ( ((TE - TAU).^2) ./ td ) ) ) ) );

    % add T2 effect
    if PARAMS.incT2
        SB = SB.*exp(-PARAMS.R2b.*TE);
    end
    
end % function SB = calcBloodCompartment(TAU,TE,PARAMS)


%% calcBloodPowder function
function SB = calcBloodPowder(TAU,TE,PARAMS)
    % Calculate intravascular contribution to ASE qBOLD signal using the powder
    % model (summarized in Yablonskiy, 2013).
    
    % pull out parameters
    dw0 = 1.5*PARAMS.dw;
    eta = sqrt((2.*dw0.*abs(TAU))./pi);
    
    SI = exp(1i.*dw0.*(TAU)./3) .* ( (fresnelc(eta) - 1i.*fresnels(eta)) ./eta);
    
    SB = abs(SI);
    
    % add T2 effect
    if PARAMS.incT2
        SB = SB.*exp(-PARAMS.R2b.*TE);
    end
    
end % function SB = calcBloodPowder(TAU,TE,PARAMS)


%% calcTissueAsymp function
function ST = calcTissueAsymp(TAU,TE,PARAMS)
    % Calculate the tissue contribution to ASE qBOLD signal using the asymptotic
    % version of the Yablonskiy (1994) static dephasing model
    
    % pull out constants
    dw   = PARAMS.dw;
    R2p  = PARAMS.R2p;
    zeta = PARAMS.zeta;
    
    % define the regime boundary
    if PARAMS.tc_man
        tc = PARAMS.tc_val;
    else
        tc = 1.76/dw;
    end

    % pre-allocate
    ST = zeros(1,length(TAU)); 
    
    % loop through tau values
    for ii = 1:length(TAU)

        if abs(TAU(ii)) < tc
            % short tau regime
            ST(ii) = exp(-(0.3*(R2p.*TAU(ii)).^2)./(zeta));

        else
            % long tau regime
            ST(ii) = exp((zeta)-(R2p*abs(TAU(ii))));

        end
    end
    
    % add T2 effect
    if PARAMS.incT2
        ST = ST .* exp(-PARAMS.R2t.*TE);
    end
    
end % function ST = calcTissueAsymp(TAU,TE,PARAMS)



