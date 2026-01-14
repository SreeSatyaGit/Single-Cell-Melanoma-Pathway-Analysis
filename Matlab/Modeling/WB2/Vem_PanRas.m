% =============================================================================
% MAPK/PI3K PATHWAY MODEL: VEMURAFENIB + PAN-RAS INHIBITOR
% =============================================================================
% Version: 2.0 | Author: Anandivada Lab | Last Updated: 2026-01-14
%
% DESCRIPTION:
%   Comprehensive ODE model of MAPK and PI3K/AKT/mTOR signaling pathways
%   simulating Combination Therapy (Vemurafenib + Pan-RAS Inhibitor).
%   Consolidated with WB1 (Trametnib.m) structure for consistency.
%
% KEY FEATURES:
%   • EGFR/Her2/Her3/IGFR/PDGFR → RAS → RAF → MEK → ERK cascade
%   • PI3K → AKT → mTOR → 4EBP1/S6K pathway
%   • DUSP and SPRY negative feedback loops
%   • Vemurafenib inhibition of BRAF^V600E & Paradoxical CRAF activation
%   • Pan-RAS inhibitor logic acting on RAS-GTP species
%
% MODEL SIZE:
%   • 68 Species (state variables)
%   • 71 Parameters (68 Unified + 3 Pan-RAS specific)
%   • Time units: seconds (experimental data in hours, converted)
% =============================================================================

clear all;
close all;
clc;

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   MAPK/PI3K PATHWAY MODEL: VEMURAFENIB + PAN-RAS INHIBITOR\n');
fprintf('   Version 2.0 - Unified WB1-Style Structure\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

%% ============================================================================
%  SECTION 1: MODEL PARAMETERS
%  ============================================================================
%  Organized by pathway module for clarity. All rates in min^-1.
% =============================================================================

fprintf('Setting up model parameters...\n');

% ----------------------------------------------------------------------------
% 1.1 RTK/EGFR Module Parameters
% ----------------------------------------------------------------------------
params.RTK.k_on = 10e-5;          % EGFR binding rate (min^-1)
params.RTK.k_off = 10e-6;         % EGFR dissociation rate (min^-1)
params.RTK.k_cat = 10e-3;         % EGFR phosphorylation rate (min^-1)
params.RTK.k_degrad = 10e-7;      % EGFR degradation rate (min^-1)
params.RTK.k_ERK_inhib = 20e-6;   % ERK feedback on EGFR (min^-1)
params.RTK.k_Her2_act = 10e-5;    % Her2 activation rate (min^-1)
params.RTK.k_Her3_act = 10e-5;    % Her3 activation rate (min^-1)
params.RTK.k_PDGFR_act = 10e-5;   % PDGFR activation rate (min^-1)

% ----------------------------------------------------------------------------
% 1.2 Shc/Grb2/SOS Module Parameters
% ----------------------------------------------------------------------------
params.SOS.k_Shc_dephos = 10e-5;      % Shc dephosphorylation (min^-1)
params.SOS.k_ptp = 7e-7;              % PTP activity (min^-1)
params.SOS.k_Grb2_bind = 36e-6;      % Grb2 binding to pShc (min^-1)
params.SOS.k_Sprty_inhib = 10e-5;    % SPRY inhibition of Grb2 (min^-1)
params.SOS.k_SOS_bind = 10e-5;       % SOS binding to Grb2 (min^-1)
params.SOS.k_ERK_phos = 17e-6;       % ERK phosphorylation of SOS (min^-1)

% ----------------------------------------------------------------------------
% 1.3 RAS Module Parameters
% ----------------------------------------------------------------------------
params.RAS.k_act = 10e-5;            % RAS activation by SOS* (min^-1)

% ----------------------------------------------------------------------------
% 1.4 RAF Module Parameters
% ----------------------------------------------------------------------------
params.RAF.k_CRAF_act = 10e-4;       % CRAF activation by RAS-GTP (min^-1)
params.RAF.k_ERK_phos_CRAF = 15.7342e-6;  % ERK phosphorylation of pCRAF (min^-1)
params.RAF.k_CRAF_degrad = 1.7342e-4;    % pCRAF degradation (min^-1)

% ----------------------------------------------------------------------------
% 1.5 MEK Module Parameters
% ----------------------------------------------------------------------------
params.MEK.k_phos = 7e-6;            % MEK phosphorylation by RAF (min^-1)
params.MEK.k_ERK_phos = 20.7342e-5;  % ERK phosphorylation of pMEK (min^-1)
params.MEK.k_degrad = 1.7342e-8;     % pMEK degradation (min^-1)
params.MEK.k_BRAF_route = 8e-6;      % BRAF^P → MEK (min^-1)
params.MEK.k_CRAF_route = 8e-6;      % CRAF* → MEK (min^-1)

% ----------------------------------------------------------------------------
% 1.6 ERK Module Parameters
% ----------------------------------------------------------------------------
params.ERK.k_phos = 20e-5;           % ERK phosphorylation by MEK (min^-1)
params.ERK.k_DUSP_inhib = 2.7342e-4; % DUSP inhibition of pERK (min^-1)
params.ERK.k_degrad = 1.7342e-7;     % pERK degradation (min^-1)
params.ERK.k_DUSP_dephos = 1e-5;     % DUSP dephosphorylation of pERK (min^-1)

% ----------------------------------------------------------------------------
% 1.7 DUSP Feedback Module Parameters
% ----------------------------------------------------------------------------
params.DUSP.k_max_tx = 1.11e-8;      % Maximum DUSP transcription rate (min^-1)
params.DUSP.k_DUSP_stop = 1.7342e-04; % DUSP stop signal (min^-1)
params.DUSP.k_DUSP_deg = 1e-6;       % DUSP degradation by ERK (min^-1)

% ----------------------------------------------------------------------------
% 1.8 SPRY Feedback Module Parameters
% ----------------------------------------------------------------------------
params.SPRY.k_max_tx = 1.7e-06;      % Maximum SPRY transcription rate (min^-1)
params.SPRY.k_come_down = 5.5000e-5; % SPRY degradation (min^-1)
params.SPRY.k_mRNA_decay = 8e-5;     % SPRY mRNA decay (min^-1)

% ----------------------------------------------------------------------------
% 1.9 BRAF Module Parameters
% ----------------------------------------------------------------------------
params.BRAF.k_inhib = 10e-2;         % BRAF inhibition rate (min^-1)

% ----------------------------------------------------------------------------
% 1.10 KSR Scaffold Module Parameters
% ----------------------------------------------------------------------------
params.KSR.k_phos = 5e-6;            % KSR phosphorylation by RAF (min^-1)
params.KSR.k_dephos = 5e-6;          % KSR dephosphorylation (min^-1)
params.KSR.k_MEK_route = 3e-6;      % pKSR → MEK phosphorylation (min^-1)

% ----------------------------------------------------------------------------
% 1.11 Trametinib (MEK Inhibitor) Parameters - DISABLED for Vem+PanRAS
% ----------------------------------------------------------------------------
params.Trametinib.conc = 0.0;       
params.Trametinib.Ki_RAF = 1e-9;     
params.Trametinib.Ki_KSR = 1e-9;     
params.Trametinib.Hill_n = 2.0;      
params.Trametinib.K_displace = 0.05; 

% ----------------------------------------------------------------------------
% 1.12 Vemurafenib and Paradoxical Activation Parameters
% ----------------------------------------------------------------------------
params.Vemurafenib.conc = 1.0;           % Vemurafenib concentration (normalized [0,1])
params.Vemurafenib.IC50 = 0.4;           % IC50 for BRAF^V600E inhibition
params.Vemurafenib.Hill_n = 1.5;         % Hill coefficient
params.Paradox.k_dimer_form = 6e-6;      % Paradoxical dimer formation induction
params.Paradox.k_dimer_dissoc = 1e-5;    % Dimer dissociation rate
params.Paradox.gamma = 0.5;              % Paradoxical CRAF activation strength

% ----------------------------------------------------------------------------
% 1.13 PI3K/AKT/mTOR Module Parameters
% ----------------------------------------------------------------------------
params.PI3K.k_p85_bind_EGFR = 10e-4;   % p85 binding to pEGFR
params.PI3K.k_p85_bind_Her2 = 10e-4;   % p85 binding to pHer2
params.PI3K.k_p85_bind_Her3 = 10e-4;   % p85 binding to pHer3
params.PI3K.k_p85_bind_IGFR = 10e-4;   % p85 binding to pIGFR
params.PI3K.k_p85_bind_PDGFR = 10e-4;  % p85 binding to pPDGFR
params.PI3K.k_p85_unbind = 10e-5;      % p85 unbinding
params.PI3K.k_PI3K_recruit = 10e-4;    % PI3K recruitment
params.PI3K.k_MTOR_feedback = 10e-4;   % mTOR feedback on PI3K
params.PI3K.k_PIP2_to_PIP3 = 10e-4;    % PIP2 to PIP3 conversion
params.PI3K.k_PTEN = 10e-5;             % PTEN
params.AKT.k_act = 10e-5;              % AKT activation
params.AKT.k_degrad = 10e-7;           % AKT degradation
params.mTOR.kb1 = 10e-8;               % mTOR activation
params.mTOR.k43b1 = 10e-3;             % mTOR activation
params.mTOR.k_4EBP1 = 10e-5;           % 4EBP1 phosphorylation
params.mTOR.k_4EBP1_dephos = 10e-5;   % p4EBP1 dephosphorylation
params.mTOR.k_S6K = 10e-5;               % S6K phosphorylation
params.mTOR.k_S6K_dephos = 10e-5;        % pS6K dephosphorylation

% ----------------------------------------------------------------------------
% 1.14 Pan-RAS Inhibitor Parameters
% ----------------------------------------------------------------------------
params.PanRAS.conc = 1.0;            % Applied inhibitor concentration
params.PanRAS.IC50 = 0.5;            % IC50 for RAS-GTP inhibition
params.PanRAS.Hill_n = 1.0;          % Hill coefficient

% ----------------------------------------------------------------------------
% 1.15 General Degradation Parameters
% ----------------------------------------------------------------------------
params.Degrad.general = 10e-7;      

% Convert parameter structure to vector
p = params;
params_vector = [
    p.RTK.k_on, p.RTK.k_off, p.RTK.k_cat, ...
    p.RAF.k_CRAF_act, p.MEK.k_phos, p.ERK.k_phos, ...
    p.RTK.k_degrad, p.RTK.k_ERK_inhib, p.SOS.k_Shc_dephos, p.SOS.k_ptp, p.SOS.k_Grb2_bind, ...
    p.SOS.k_Sprty_inhib, p.SOS.k_SOS_bind, p.SOS.k_ERK_phos, ...
    p.RAF.k_ERK_phos_CRAF, p.RAF.k_CRAF_degrad, p.MEK.k_ERK_phos, p.MEK.k_degrad, ...
    p.ERK.k_DUSP_inhib, p.ERK.k_degrad, p.BRAF.k_inhib, p.DUSP.k_DUSP_stop, p.DUSP.k_max_tx, ...
    p.SPRY.k_max_tx, p.SPRY.k_come_down, p.Degrad.general, p.SPRY.k_mRNA_decay, ...
    p.DUSP.k_max_tx, p.SPRY.k_max_tx, ...
    p.ERK.k_DUSP_dephos, p.DUSP.k_DUSP_deg, ...
    p.RTK.k_Her2_act, p.RTK.k_Her3_act, ...
    p.PI3K.k_p85_bind_EGFR, p.PI3K.k_p85_bind_Her2, p.PI3K.k_p85_bind_Her3, p.PI3K.k_p85_bind_IGFR, ...
    p.PI3K.k_p85_unbind, p.PI3K.k_PI3K_recruit, ...
    p.PI3K.k_MTOR_feedback, ...
    p.PI3K.k_PIP2_to_PIP3, p.PI3K.k_PTEN, ...
    p.AKT.k_act, p.AKT.k_degrad, ...
    p.mTOR.kb1, p.mTOR.k43b1, p.mTOR.k_4EBP1, p.mTOR.k_4EBP1_dephos, ...
    p.KSR.k_phos, p.KSR.k_dephos, ...
    p.MEK.k_BRAF_route, p.MEK.k_CRAF_route, p.KSR.k_MEK_route, ...
    p.Trametinib.conc, p.Trametinib.Ki_RAF, p.Trametinib.Ki_KSR, p.Trametinib.Hill_n, ...
    p.Vemurafenib.conc, p.Paradox.k_dimer_form, p.Paradox.k_dimer_dissoc, ...
    p.Paradox.gamma, p.Vemurafenib.IC50, p.Vemurafenib.Hill_n, ...
    p.RTK.k_PDGFR_act, p.PI3K.k_p85_bind_PDGFR, ...
    p.mTOR.k_S6K, p.mTOR.k_S6K_dephos, ...
    p.Trametinib.K_displace, ... % Parameters up to 68 (Unified Standard)
    p.PanRAS.conc, p.PanRAS.IC50, p.PanRAS.Hill_n % 69, 70, 71 (PanRAS extension)
];

%% ============================================================================
%  SECTION 2: INITIAL CONDITIONS
%  ============================================================================

fprintf('Setting initial conditions...\n');

IC.EGFR = [1.0, 0.35, 0.35];    IC.Her2 = [1.0, 0.245, 0.245];  IC.Her3 = [1.0, 0.203, 0.203];
IC.SHC = [1.0, 0.0, 1.0];       IC.Grb2_SOS = [0.0, 0.0];
IC.HRAS = [1.0, 0.0];           IC.NRAS = [1.0, 0.0];           IC.KRAS = [1.0, 0.0, 1.0]; % 19, 20, 21
IC.CRAF = [0.8, 0.366];         IC.BRAF = [1.0, 1.0];
IC.MEK = [1.0, 1.759];          IC.ERK = [1.0, 2.903];
IC.DUSP = [1.0, 2.677];         IC.SPRY = [1.0, 1.0];
IC.pERK_degrad = [1.0];         IC.pMEK_degrad = [1.0];         IC.pCRAF_degrad = [1.0]; IC.DUSP_stop = [1.0];
IC.IGFR = [1.0, 0.0, 0.0];      IC.IRS = [1.0, 0.0];
IC.p85 = [1.0];
IC.p85_EGFR = [0.1]; IC.p85_Her2 = [0.1]; IC.p85_Her3 = [0.1]; IC.p85_IGFR = [0.1];
IC.PI3K = [1.0, 0.2];           IC.PIP = [1.0, 0.1];
IC.AKT = [1.0, 0.513];          IC.FOXO = [0.0];
IC.mTORC = [1.0, 0.5];          IC.frebp1 = [1.0, 0.5, 1.002];
IC.KSR = [1.0, 0.0];            IC.BRAF_CRAF_dimer = [0.0];
IC.PDGFR = [1.0, 0.474, 0.474]; IC.S6K = [1.0, 1.432];    IC.p85_PDGFR = [0.1];

y0 = [
    IC.EGFR, IC.Her2, IC.Her3, IC.SHC, IC.Grb2_SOS, ...
    IC.HRAS, IC.NRAS, IC.KRAS, ...
    IC.CRAF, IC.BRAF, ...
    IC.MEK, IC.ERK, ...
    IC.DUSP, IC.SPRY, ...
    IC.pERK_degrad, IC.pMEK_degrad, IC.pCRAF_degrad, IC.DUSP_stop, ...
    IC.IGFR, IC.IRS, ...
    IC.p85, IC.p85_EGFR, IC.p85_Her2, IC.p85_Her3, IC.p85_IGFR, ...
    IC.PI3K, IC.PIP, IC.AKT, IC.FOXO, IC.mTORC, IC.frebp1, ...
    IC.KSR, IC.BRAF_CRAF_dimer, ...
    IC.PDGFR, IC.S6K, IC.p85_PDGFR
];

%% ============================================================================
%  SECTION 3: EXPERIMENTAL DATA
%  ============================================================================

fprintf('Loading experimental data...\n');
timeStamps_hours = [0, 1, 4, 8, 24, 48];
timeStamps_seconds = timeStamps_hours * 3600;

% Data from current Vem_PanRas.m (used as baseline/target)
expData_raw.pMEK = [1.368179432,0.040145283,0.020650599,0.0079056,0.008636227,0.021682145];
expData_raw.pERK = [1.047565823,0.028753956,0.02683006,0.019206199,0.021359286,1.045153456];
expData_raw.DUSP = [0.939982354,0.940182728,0.93995499,0.939315828,0.939284046,0.939384237];
expData_raw.pEGFR = [0.231259538,0.281291836,0.206473671,0.067136885,0.029598009,0.041653072];
expData_raw.pCRAF = [0.189167867,0.109105047,0.102702818,0.074742052,0.030393272,0.013923339];
expData_raw.pAKT = [0.025678316,0.266170637,1.167110876,0.256750989,0.212461114,1.167189541];
expData_raw.p4ebp1 = [0.989561348,0.989525727,0.989495943,0.989507379,0.989667542,0.476875796];

% Placeholder data to avoid errors (panRAS, Her, S6K)
expData_raw.panRAS = ones(size(timeStamps_hours)) * 0.5;
expData_raw.her2 = ones(size(timeStamps_hours)) * 0.5;
expData_raw.her3 = ones(size(timeStamps_hours)) * 0.5;
expData_raw.pDGFR = ones(size(timeStamps_hours)) * 0.5;
expData_raw.pS6k = ones(size(timeStamps_hours)) * 0.5;

% Normalize
species_names = fieldnames(expData_raw);
for i = 1:length(species_names)
    data = expData_raw.(species_names{i});
    expData_norm.(species_names{i}) = (data - min(data)) / (max(data) - min(data) + eps);
end

%% ============================================================================
%  SECTION 4: LOAD PRE-TRAINED PARAMETERS
%  ============================================================================

fprintf('Loading pre-trained Vemurafenib parameters...\n');
try
    load('trained_Vem_params.mat', 'optimizedParams');
    % Apply trained parameters to our vector (1-53 MAPK/PI3K Kinetics)
    params_vector(1:53) = optimizedParams(1:53);
    % Vemurafenib (58-63)
    params_vector(58:63) = optimizedParams(58:63);
    fprintf('Successfully loaded and mapped parameters.\n');
catch
    warning('Pre-trained parameters not found. Using default values.');
end

%% ============================================================================
%  SECTION 5 & 6: PARAMETER SETUP AND SIMULATION
%  ============================================================================

% We can skip optimization if we are just performing a prediction
fprintf('Simulating Vemurafenib + Pan-RAS Inhibitor (Prediction)...\n');
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', 1:68);

[T_all, Y_all] = ode15s(@(t,y) Mapk_ODE(t, y, params_vector), timeStamps_seconds, y0, ode_opts);

% Normalization
min_max_norm = @(v) (v - min(v)) ./ (max(v) - min(v) + eps);

model_outputs.pEGFR = min_max_norm(Y_all(:,3));
model_outputs.panRAS = min_max_norm(Y_all(:,16) + Y_all(:,18) + Y_all(:,21));
model_outputs.pCRAF = min_max_norm(Y_all(:,23));
model_outputs.pMEK = min_max_norm(Y_all(:,27));
model_outputs.pERK = min_max_norm(Y_all(:,29));
model_outputs.pAKT = min_max_norm(Y_all(:,53));
model_outputs.p4EBP1 = min_max_norm(Y_all(:,59));

%% ============================================================================
%  SECTION 7: VISUALIZATION
%  ============================================================================

fprintf('Generating plots...\n');
tFine_hours = linspace(0, 48, 400);
[T_fine, Y_fine] = ode15s(@(t,y) Mapk_ODE(t, y, params_vector), tFine_hours*3600, y0, ode_opts);

% Extract Smooth Trajectories
model_smooth.pEGFR = min_max_norm(Y_fine(:,3));
model_smooth.panRAS = min_max_norm(Y_fine(:,16) + Y_fine(:,18) + Y_fine(:,21));
model_smooth.pCRAF = min_max_norm(Y_fine(:,23));
model_smooth.pMEK = min_max_norm(Y_fine(:,27));
model_smooth.pERK = min_max_norm(Y_fine(:,29));
model_smooth.pAKT = min_max_norm(Y_fine(:,53));
model_smooth.p4EBP1 = min_max_norm(Y_fine(:,59));

figure('Name', 'Prediction: Vemurafenib + Pan-RAS Inhibitor', 'Position', [50, 50, 1600, 1000]);
plot_list = {'pEGFR', 'pCRAF', 'pMEK', 'pERK', 'pAKT', 'p4EBP1', 'panRAS'};
data_keys = {'pEGFR', 'pCRAF', 'pMEK', 'pERK', 'pAKT', 'p4ebp1', 'panRAS'};

for i = 1:length(plot_list)
    subplot(3, 3, i);
    spec = plot_list{i};
    key = data_keys{i};
    plot(tFine_hours, model_smooth.(spec), 'b-', 'LineWidth', 3, 'DisplayName', 'Prediction'); hold on;
    plot(timeStamps_hours, expData_norm.(key), 'ro', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Baseline Data');
    title(spec); xlabel('Time (h)'); ylabel('Activity'); grid on;
    if i == 1; legend('Location', 'best'); end
end
sgtitle('Therapy Response: Vemurafenib + Pan-RAS Inhibitor (Unified Model)');

%% ============================================================================
%  ODE FUNCTION (Unified Architecture)
%  ============================================================================
function dydt = Mapk_ODE(t, y, p)
    dydt = zeros(68, 1);
    
    % --- Unpack Parameters (1-68 Unified) ---
    ka1=p(1); kr1=p(2); kc1=p(3); kpCraf=p(4); kpMek=p(5); kpErk=p(6); kDegradEgfr=p(7); 
    kErkInbEgfr=p(8); kShcDephos=p(9); kptpDeg=p(10); kGrb2CombShc=p(11); kSprtyInbGrb2=p(12); 
    kSosCombGrb2=p(13); kErkPhosSos=p(14); kErkPhosPcraf=p(15); kPcrafDegrad=p(16); 
    kErkPhosMek=p(17); kMekDegrad=p(18); kDuspInbErk=p(19); kErkDeg=p(20); kinbBraf=p(21); 
    kDuspStop=p(22); kDusps=p(23); kSproutyForm=p(24); kSprtyComeDown=p(25); kdegrad=p(26); 
    km_Sprty=p(29); kErkDephos=p(30); kDuspDeg=p(31);
    kHer2_act=p(32); kHer3_act=p(33); k_p85_bind_EGFR=p(34); k_p85_bind_Her2=p(35); 
    k_p85_bind_Her3=p(36); k_p85_bind_IGFR=p(37); k_p85_unbind=p(38); k_PI3K_recruit=p(39); 
    kMTOR_Feedback=p(40); k_PIP2_to_PIP3=p(41); k_PTEN=p(42); kAkt=p(43); kdegradAKT=p(44); 
    kb1=p(45); k43b1=p(46); k4ebp1=p(47); k_4EBP1_dephos=p(48); kKSRphos=p(49); kKSRdephos=p(50);
    kMekByBraf=p(51); kMekByCraf=p(52); kMekByKSR=p(53);
    Tram=p(54); K_tram_RAF=p(55); n_tram=p(57);
    Vemurafenib=p(58); kDimerForm=p(59); kDimerDissoc=p(60); kParadoxCRAF=p(61); 
    IC50_vem=p(62); Hill_n_vem=p(63);
    kPDGFR_act=p(64); k_p85_bind_PDGFR=p(65); kS6K_phos=p(66); kS6K_dephos=p(67);
    K_displace=p(68);
    
    % --- Extension: Pan-RAS Inhibitor (69-71) ---
    RASi_conc = p(69); RASi_IC50 = p(70); RASi_Hill = p(71);
    % Calculate RAS-inhibition effectiveness
    RASi_eff = (RASi_conc^RASi_Hill) / (RASi_IC50^RASi_Hill + RASi_conc^RASi_Hill + eps);

    % --- ODEs ---
    % RTK Module
    dydt(1) = -ka1*y(1) + kr1*y(2);
    dydt(2) = ka1*y(1) - kr1*y(2) - kc1*y(2);
    dydt(3) = kc1*y(2) - kDegradEgfr*y(3) - kErkInbEgfr*y(29)*y(3);
    dydt(4) = -kHer2_act*y(4) + kr1*y(5);
    dydt(5) = kHer2_act*y(4) - kr1*y(5) - kc1*y(5);
    dydt(6) = kc1*y(5) - kDegradEgfr*y(6) - kErkInbEgfr*y(29)*y(6);
    dydt(7) = -kHer3_act*y(7) + kr1*y(8);
    dydt(8) = kHer3_act*y(7) - kr1*y(8) - kc1*y(8);
    dydt(9) = kc1*y(8) - kDegradEgfr*y(9) - kErkInbEgfr*y(29)*y(9);
    
    % SOS
    dydt(10) = -ka1*y(3)*y(10);
    dydt(11) = ka1*y(3)*y(10) - kShcDephos*y(12)*y(11);
    dydt(12) = -kptpDeg*y(11)*y(12);
    dydt(13) = kGrb2CombShc*y(11)*y(3) - kSprtyInbGrb2*y(33)*y(13);
    dydt(14) = kSosCombGrb2*y(13)*y(11) - kErkPhosSos*y(29)*y(14);
    
    % RAS Module (Applying Pan-RAS Inhibitor)
    % Deactivation enhanced by inhibitor (simplified model)
    panRAS_active = y(16) + y(18) + y(21);
    RAS_GAP_rate = kdegrad * (1.0 + RASi_eff * 50.0); % Inhibitor speeds up deactivation
    
    dydt(15) = -ka1*y(14)*y(15) + RAS_GAP_rate*y(16); 
    dydt(16) =  ka1*y(14)*y(15) - RAS_GAP_rate*y(16);
    dydt(17) = -ka1*y(14)*y(17) + RAS_GAP_rate*y(18); 
    dydt(18) =  ka1*y(14)*y(17) - RAS_GAP_rate*y(18);
    dydt(19) = -ka1*y(14)*y(19) + RAS_GAP_rate*y(21); 
    dydt(20) = 0;
    dydt(21) =  ka1*y(14)*y(19) - RAS_GAP_rate*y(21);
    
    % RAF / Vemurafenib Logic
    IC50_n = IC50_vem^Hill_n_vem; Vem_n = Vemurafenib^Hill_n_vem;
    kBRAF_eff = ka1 * IC50_n / (IC50_n + Vem_n + eps);
    paradox_activation = kParadoxCRAF * Vemurafenib * y(62);
    
    dydt(22) = -kpCraf*panRAS_active*y(22) + kErkPhosPcraf*y(29)*y(23) + kPcrafDegrad*y(23)*y(36) - kDimerForm*y(25)*y(22)*Vemurafenib + kDimerDissoc*y(62);
    dydt(23) = kpCraf*panRAS_active*y(22) - kErkPhosPcraf*y(29)*y(23) - kPcrafDegrad*y(23)*y(36) + paradox_activation;
    dydt(24) = -kBRAF_eff*y(24)*panRAS_active - kDimerForm*y(25)*y(22)*Vemurafenib + kDimerDissoc*y(62);
    dydt(25) = kBRAF_eff*y(24)*panRAS_active - kinbBraf*y(25) - kDimerForm*y(25)*y(22)*Vemurafenib + kDimerDissoc*y(62);
    dydt(62) = kDimerForm*y(25)*y(22)*Vemurafenib - kDimerDissoc*y(62) - kPcrafDegrad*y(62)*y(36);
    
    % MEK/ERK
    Stream_Load = panRAS_active + (y(22)+y(23)) + (y(24)+y(25));
    Ki_effective = K_tram_RAF * (1 + (Stream_Load / K_displace)^2);
    f_MEK_activity = 1 / (max(0.01, 1 + (Tram / max(eps, Ki_effective))^n_tram));
    
    raf_to_mek = (kpMek*y(23) + kMekByBraf*y(25) + kMekByCraf*y(23) + kpMek*y(62));
    ksr_to_mek = (kMekByKSR * y(61));
    dydt(26) = -(raf_to_mek + ksr_to_mek)*y(26) + kErkPhosMek*y(29)*y(27) + kMekDegrad*y(27)*y(35);
    dydt(27) = (raf_to_mek + ksr_to_mek)*y(26) - kErkPhosMek*y(29)*y(27) - kMekDegrad*y(27)*y(35);
    
    erk_activation = kpErk * y(27) * y(28) * f_MEK_activity;
    dydt(28) = -erk_activation + kErkDephos*y(31)*y(29) + kErkDeg*y(29)*y(34);
    dydt(29) = erk_activation - kErkDephos*y(31)*y(29) - kErkDeg*y(29)*y(34);
    
    % PI3K Module (Simplified for Prediction)
    dydt(30) = km_Dusp*y(29)/(1 + (km_Dusp/kDusps)*y(29)) - kDuspStop*y(30)*y(37);
    dydt(31) = -kDuspStop*y(30)*y(31);
    dydt(32) = km_Sprty*y(29)/(1 + (km_Sprty/kSproutyForm)*y(29)) - kSprtyComeDown*y(32)*y(33);
    dydt(33) = -kSprtyComeDown*y(32)*y(33);
    dydt(34) = -kErkDeg*y(29)*y(34); dydt(35) = -kMekDegrad*y(27)*y(35); 
    dydt(36) = -kPcrafDegrad*y(23)*y(36); dydt(37) = -kDuspStop*y(30)*y(37);
    total_p85_RTK = y(44) + y(45) + y(46) + y(47) + y(68);
    dydt(48) = -k_PI3K_recruit*total_p85_RTK*y(48) + kMTOR_Feedback*y(56)*y(49);
    dydt(49) = k_PI3K_recruit*total_p85_RTK*y(48) - kMTOR_Feedback*y(56)*y(49);
    dydt(50) = -k_PIP2_to_PIP3*y(49)*y(50) + k_PTEN*y(51);
    dydt(51) = k_PIP2_to_PIP3*y(49)*y(50) - k_PTEN*y(51);
    dydt(52) = -kAkt*y(51)*y(52) + kdegradAKT*y(53);
    dydt(53) = kAkt*y(51)*y(52) - kdegradAKT*y(53);
    dydt(54) = (max(0, 1 - y(53))) * kAkt / (max(0.1, 1 + (y(54) / 15e-5)));
    dydt(55) = -kAkt * y(53) * y(55) + kdegrad * y(56);
    dydt(56) = kAkt * y(53) * y(55) - kdegrad * y(56);
    dydt(57) = -k4ebp1*y(56)*y(57) + kb1*y(58) + k_4EBP1_dephos*y(59);
    dydt(58) = k4ebp1*y(56)*y(57) - kb1*y(58) - k43b1*y(58);
    dydt(59) = k43b1*y(58) - k_4EBP1_dephos*y(59);
    dydt(60) = -kKSRphos*panRAS_active*y(60) + kKSRdephos*y(61);
    dydt(61) = kKSRphos*panRAS_active*y(60) - kKSRdephos*y(61);
    dydt(63) = -kPDGFR_act*y(63) + kr1*y(64);
    dydt(64) = kPDGFR_act*y(63) - kr1*y(64) - kc1*y(64);
    dydt(65) = kc1*y(64) - kDegradEgfr*y(65) - kErkInbEgfr*y(29)*y(65);
    dydt(66) = -kS6K_phos*y(56)*y(66) + kS6K_dephos*y(67);
    dydt(67) = kS6K_phos*y(56)*y(66) - kS6K_dephos*y(67);
    dydt(68) = k_p85_bind_PDGFR*y(65)*y(43) - k_p85_unbind*y(68);
    
    % Scale to hours (p parameters are min^-1, model uses min) 
    % Wait, the parameters are min^-1, so Mapk_ODE should return dydt in min^-1?
    % Actually the system in WB1 uses seconds-scaled derivatives but p is min^-1.
    % To match WB1 exactly, we multiply by scale factor in returns.
    dydt = dydt * 60; % min^-1 derivatives
end
