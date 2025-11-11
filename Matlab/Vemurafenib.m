% =============================================================================
% MAPK/PI3K PATHWAY MODEL WITH VEMURAFENIB AND PARADOXICAL ACTIVATION
% =============================================================================
%
% Model Description:
% This script implements a comprehensive ODE model of the MAPK and PI3K/AKT/mTOR
% signaling pathways with vemurafenib treatment and paradoxical activation.
%
% Key Features:
% - EGFR → RAS → RAF → MEK → ERK cascade
% - PI3K → AKT → mTOR pathway
% - DUSP and SPRY negative feedback loops
% - Vemurafenib inhibition of BRAF^V600E (Hill equation)
% - Paradoxical activation: vemurafenib-induced CRAF activation via dimerization
% - KSR scaffold protein for MEK phosphorylation
% - Parameter optimization to fit experimental data
%
% Model Structure:
% - 49 species (48 original + 1 dimer complex)
% - 52 parameters (46 original + 6 paradoxical activation)
% - Time units: seconds (experimental data in hours, converted)
%
% =============================================================================

clear all;
close all;
clc;

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   MAPK/PI3K PATHWAY MODEL WITH VEMURAFENIB\n');
fprintf('   Paradoxical Activation Included\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

%% ============================================================================
% SECTION 1: MODEL PARAMETERS
% ============================================================================
% Organized by pathway/module for clarity

fprintf('Setting up model parameters...\n');

% ----------------------------------------------------------------------------
% 1.1 RTK/EGFR Module Parameters
% ----------------------------------------------------------------------------
params.RTK.k_on = 10e-5;          % EGFR binding rate (min^-1)
params.RTK.k_off = 10e-6;         % EGFR dissociation rate (min^-1)
params.RTK.k_cat = 10e-3;         % EGFR phosphorylation rate (min^-1)
params.RTK.k_degrad = 10e-7;      % EGFR degradation rate (min^-1)
params.RTK.k_ERK_inhib = 20e-6;   % ERK feedback on EGFR (min^-1)

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
% Note: RAS isoforms (HRAS, NRAS, KRAS) share the same activation rate

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
params.MEK.k_degrad = 1.7342e-4;     % pMEK degradation (min^-1)
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
params.DUSP.k_max_tx = 1.11e-4;      % Maximum DUSP transcription rate (min^-1)
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
% 1.11 Trametinib (MEK Inhibitor) Parameters
% ----------------------------------------------------------------------------
params.Trametinib.conc = 1e-6;       % Trametinib concentration
params.Trametinib.Ki_RAF = 5e-1;     % Ki for RAF→MEK inhibition
params.Trametinib.Ki_KSR = 5e-12;    % Ki for KSR→MEK inhibition
params.Trametinib.Hill_n = 2;        % Hill coefficient

% ----------------------------------------------------------------------------
% 1.12 Vemurafenib and Paradoxical Activation Parameters
% ----------------------------------------------------------------------------
params.Vemurafenib.conc = 1.0;           % Vemurafenib concentration (normalized [0,1])
params.Vemurafenib.IC50 = 0.4;           % IC50 for BRAF^V600E inhibition (normalized [0,1])
params.Vemurafenib.Hill_n = 1.5;         % Hill coefficient for inhibition
params.Paradox.k_dimer_form = 6e-6;      % Dimer formation rate (BRAF-WT* + CRAF → dimer)
params.Paradox.k_dimer_dissoc = 1e-5;    % Dimer dissociation rate (min^-1)
params.Paradox.gamma = 0.5;              % Paradoxical CRAF activation strength

% ----------------------------------------------------------------------------
% 1.13 PI3K/AKT/mTOR Module Parameters
% ----------------------------------------------------------------------------
params.PI3K.k_EGFR_act = 10e-4;     % PI3K activation by EGFR (min^-1)
params.PI3K.k_MTOR_feedback = 10e-4; % mTOR feedback on PI3K (min^-1)
params.AKT.k_act = 10e-5;           % AKT activation rate (min^-1)
params.AKT.k_degrad = 10e-7;        % AKT degradation (min^-1)
params.mTOR.kb1 = 10e-8;            % mTOR activation parameter
params.mTOR.k43b1 = 10e-3;          % mTOR activation parameter

% ----------------------------------------------------------------------------
% 1.14 General Degradation Parameters
% ----------------------------------------------------------------------------
params.Degrad.general = 10e-7;      % General degradation rate (min^-1)

% Convert parameter structure to vector for optimization
% Order must match objectiveFunction_all parameter unpacking
p = params;
params_vector = [
    p.RTK.k_on, p.RTK.k_off, p.RTK.k_cat, ...
    p.RAF.k_CRAF_act, p.MEK.k_phos, p.ERK.k_phos, ...
    p.RTK.k_degrad, p.RTK.k_ERK_inhib, p.SOS.k_Shc_dephos, p.SOS.k_ptp, p.SOS.k_Grb2_bind, ...
    p.SOS.k_Sprty_inhib, p.SOS.k_SOS_bind, p.SOS.k_ERK_phos, ...
    p.RAF.k_ERK_phos_CRAF, p.RAF.k_CRAF_degrad, p.MEK.k_ERK_phos, p.MEK.k_degrad, ...
    p.ERK.k_DUSP_inhib, p.ERK.k_degrad, p.BRAF.k_inhib, p.DUSP.k_DUSP_stop, p.DUSP.k_max_tx, ...
    p.SPRY.k_max_tx, p.SPRY.k_come_down, p.Degrad.general, p.SPRY.k_mRNA_decay, p.DUSP.k_max_tx, p.SPRY.k_max_tx, ...
    p.ERK.k_DUSP_dephos, p.DUSP.k_DUSP_deg, ...
    p.PI3K.k_EGFR_act, p.PI3K.k_MTOR_feedback, ...
    p.AKT.k_act, p.AKT.k_degrad, ...
    p.mTOR.kb1, p.mTOR.k43b1, ...
    p.KSR.k_phos, p.KSR.k_dephos, ...
    p.MEK.k_BRAF_route, p.MEK.k_CRAF_route, p.KSR.k_MEK_route, ...
    p.Trametinib.conc, p.Trametinib.Ki_RAF, p.Trametinib.Ki_KSR, p.Trametinib.Hill_n, ...
    p.Vemurafenib.conc, p.Paradox.k_dimer_form, p.Paradox.k_dimer_dissoc, ...
    p.Paradox.gamma, p.Vemurafenib.IC50, p.Vemurafenib.Hill_n
];

fprintf('Parameters organized by pathway modules.\n\n');

%% ============================================================================
% SECTION 2: INITIAL CONDITIONS
% ============================================================================

fprintf('Setting initial conditions...\n');

% Define initial conditions for each species
% Format: [inactive, intermediate, active] or [inactive, active]

% RTK/EGFR Module
IC.EGFR = [1.0, 0.0, 0.0];          % EGFR, EGFR:ligand, pEGFR
IC.SHC = [1.0, 0.0, 1.0];           % Shc, pEGFR:Shc, pShc
IC.Grb2_SOS = [0.0, 0.0];           % Grb2, pShc:Grb2:SOS

% RAS Module
IC.HRAS = [0.0, 0.0];               % HRAS-GDP, HRAS-GTP
IC.NRAS = [0.0, 0.0];               % NRAS-GDP, NRAS-GTP
IC.KRAS = [1.0, 0.0, 1.0];          % KRAS-GDP, intermediate, KRAS-GTP

% RAF Module
IC.CRAF = [0.8, 0.2];               % CRAF, pCRAF
IC.BRAF = [0.0, 1.0];               % BRAF, BRAF^P (BRAF^V600E*)

% MEK/ERK Module
IC.MEK = [1.0, 1.0];                % MEK, pMEK
IC.ERK = [1.0, 0.8];                % ERK, pERK

% Feedback Module
IC.DUSP = [1.0, 1.0];               % DUSP mRNA, DUSP protein
IC.SPRY = [1.0, 1.0];               % SPRY mRNA, SPRY protein

% Degradation Trackers
IC.pERK_degrad = [1.0];
IC.pMEK_degrad = [1.0];
IC.pCRAF_degrad = [1.0];
IC.DUSP_stop = [1.0];

% PI3K/AKT/mTOR Module
IC.IGFR = [1.0, 0.0, 0.0];          % IGFR, IGFR:ligand, pIGFR
IC.IRS = [1.0, 0.0];                % IRS, pIRS
IC.PI3K = [1.0, 0.0];               % PI3K, PI3K_active
IC.AKT = [1.0, 0.0];                % AKT, pAKT
IC.FOXO = [0.0];
IC.mTORC = [1.0, 0.0];              % mTORC, mTORC_active
IC.frebp1 = [1.0, 0.0, 0.0];        % 4EBP1, intermediate, p4EBP1

% KSR Module
IC.KSR = [1.0, 0.0];                % KSR, pKSR

% Paradoxical Activation
IC.BRAF_CRAF_dimer = [0.0];         % BRAF-WT:CRAF dimer

% Combine all initial conditions into state vector
y0 = [
    IC.EGFR, IC.SHC, IC.Grb2_SOS, ...
    IC.HRAS, IC.NRAS, IC.KRAS, ...
    IC.CRAF, IC.BRAF, ...
    IC.MEK, IC.ERK, ...
    IC.DUSP, IC.SPRY, ...
    IC.pERK_degrad, IC.pMEK_degrad, IC.pCRAF_degrad, IC.DUSP_stop, ...
    IC.IGFR, IC.IRS, IC.PI3K, IC.AKT, IC.FOXO, IC.mTORC, IC.frebp1, ...
    IC.KSR, IC.BRAF_CRAF_dimer
];

fprintf('Initial conditions set. Total species: %d\n\n', length(y0));

%% ============================================================================
% SECTION 3: EXPERIMENTAL DATA
% ============================================================================

fprintf('Loading experimental data...\n');

% Experimental time points (hours)
timeStamps_hours = [0, 1, 4, 8, 24, 48];
timeStamps_seconds = timeStamps_hours * 3600;  % Convert to seconds

% Raw experimental data (not normalized)
expData_raw.RAS_GTP = [0.558156831, 0.61832384, 0.614585492, 0.708019641, 0.999240675, 1.0];
expData_raw.pMEK = [1.75938884, 0.170160085, 0.095112609, 0.201000276, 0.219207054, 0.502831668];
expData_raw.pERK = [2.903209735, 0.207867788, 0.303586121, 0.805254439, 1.408362153, 1.847606441];
expData_raw.DUSP = [2.677161325, 2.782754577, 1.130758062, 0.395642757, 0.828575853, 0.916618219];
expData_raw.pEGFR = [0.291928893, 0.392400458, 0.265016688, 0.394238749, 0.006158316, 0.008115099];
expData_raw.pCRAF = [0.366397596, 0.537106733, 0.465541704, 0.586732657, 1.102322681, 0.269181259];
expData_raw.pAKT = [0.513544148, 0.613178403, 1.03451863, 1.113391047, 0.535242724, 0.538273551];

% Note: RAS_GTP data is the same as RASGtpValsExp_all from original code

% Normalize experimental data to [0, 1] using min-max normalization
species_names = fieldnames(expData_raw);
expData_norm = struct();
for i = 1:length(species_names)
    data = expData_raw.(species_names{i});
    data_min = min(data);
    data_max = max(data);
    expData_norm.(species_names{i}) = (data - data_min) / (data_max - data_min + eps);
end

fprintf('Experimental data normalized. Time points: %d\n\n', length(timeStamps_hours));

%% ============================================================================
% SECTION 4: PARAMETER OPTIMIZATION SETUP
% ============================================================================

fprintf('Setting up parameter optimization...\n');

% Initial parameter guess (from params_vector)
params0 = params_vector;

% Parameter bounds
num_params = length(params0);
lb = 1e-8 * ones(num_params, 1);
ub = 1e-4 * ones(num_params, 1);

% Enhanced bounds for specific parameters
% DUSP-ERK interaction (parameter 30)
lb(30) = 1e-6;  % kErkDephos
ub(30) = 1e-3;

% DUSP degradation (parameter 31)
lb(31) = 1e-7;  % kDuspDeg
ub(31) = 1e-4;

% Paradoxical activation parameters (parameters 47-52)
lb(47) = 0.0;      % Vemurafenib concentration [0, 1]
ub(47) = 1.0;
lb(48) = 1e-7;     % kDimerForm
ub(48) = 1e-3;
lb(49) = 1e-6;     % kDimerDissoc
ub(49) = 1e-3;
lb(50) = 0.1;      % kParadoxCRAF
ub(50) = 2.0;
lb(51) = 0.01;     % IC50_vem
ub(51) = 0.9;
lb(52) = 0.5;      % Hill_n_vem
ub(52) = 3.0;

% Optimization options
opts = optimoptions(@fmincon, ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxIterations', 150, ...
    'FunctionTolerance', 1e-6, ...
    'StepTolerance', 1e-8);

fprintf('Optimization setup complete. Parameters to optimize: %d\n\n', num_params);

%% ============================================================================
% SECTION 5: RUN PARAMETER OPTIMIZATION
% ============================================================================

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   RUNNING PARAMETER OPTIMIZATION\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

tic;
[optimizedParams, errorOpt] = fmincon( ...
    @(p) objectiveFunction_all(p, timeStamps_seconds, expData_norm, y0), ...
    params0, [], [], [], [], lb, ub, [], opts);
optimization_time = toc;

fprintf('\nOptimization completed in %.2f seconds.\n', optimization_time);
fprintf('Optimized Parameters:\n');
disp(optimizedParams);
fprintf('Minimum Fit Error: %.6e\n\n', errorOpt);

%% ============================================================================
% SECTION 6: SIMULATE WITH OPTIMIZED PARAMETERS
% ============================================================================

fprintf('Simulating with optimized parameters...\n');

% Simulate at experimental time points
[T_all, Y_all] = ode45(@(t,y) Mapk_ODE(t, y, optimizedParams), ...
                       timeStamps_seconds, y0);

% Normalization helper function
normit = @(v, mode) v ./ max(eps, ...
    (strcmpi(mode,'first'))*v(1) + ...
    (strcmpi(mode,'last')) *v(end) + ...
    (strcmpi(mode,'max'))  *max(v));

% Extract and normalize model outputs
model_outputs.pEGFR = normit(Y_all(:,3), 'max');
model_outputs.RAS_GTP = normit(Y_all(:,14), 'last');
model_outputs.pCRAF = normit(Y_all(:,17), 'max');
model_outputs.pMEK = normit(Y_all(:,21), 'first');
model_outputs.pERK = normit(Y_all(:,23), 'first');
model_outputs.DUSP = normit(Y_all(:,24), 'max');
model_outputs.pAKT = normit(Y_all(:,40), 'max');

% Calculate fit errors
fit_errors = struct();
fit_errors.pEGFR = abs(model_outputs.pEGFR - expData_norm.pEGFR(:));
fit_errors.RAS = abs(model_outputs.RAS_GTP - expData_norm.RAS_GTP(:));
fit_errors.pCRAF = abs(model_outputs.pCRAF - expData_norm.pCRAF(:));
fit_errors.pMEK = abs(model_outputs.pMEK - expData_norm.pMEK(:));
fit_errors.pERK = abs(model_outputs.pERK - expData_norm.pERK(:));
fit_errors.DUSP = abs(model_outputs.DUSP - expData_norm.DUSP(:));
fit_errors.pAKT = abs(model_outputs.pAKT - expData_norm.pAKT(:));

% Display fit results
fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   MODEL FIT RESULTS\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('Total Fit Error: %.4f\n', errorOpt);
fprintf('Individual Protein Errors (sum across all timepoints):\n');
fprintf('  pEGFR: %.4f\n', sum(fit_errors.pEGFR));
fprintf('  RAS:   %.4f\n', sum(fit_errors.RAS));
fprintf('  pCRAF: %.4f\n', sum(fit_errors.pCRAF));
fprintf('  pMEK:  %.4f\n', sum(fit_errors.pMEK));
fprintf('  pERK:  %.4f\n', sum(fit_errors.pERK));
fprintf('  DUSP:  %.4f\n', sum(fit_errors.DUSP));
fprintf('  pAKT:  %.4f\n', sum(fit_errors.pAKT));
fprintf('\n');

%% ============================================================================
% SECTION 7: VISUALIZATION
% ============================================================================

fprintf('Generating visualization plots...\n');

% Generate smooth time course for plotting
tFine_hours = linspace(0, 48, 400);
tFine_seconds = tFine_hours * 3600;

[T_fine, Y_fine] = ode45(@(t,y) Mapk_ODE(t, y, optimizedParams), ...
                         tFine_seconds, y0);

% Normalize smooth model outputs
model_smooth.pEGFR = normit(Y_fine(:,3), 'max');
model_smooth.RAS_GTP = normit(Y_fine(:,14), 'last');
model_smooth.pCRAF = normit(Y_fine(:,17), 'max');
model_smooth.pMEK = normit(Y_fine(:,21), 'first');
model_smooth.pERK = normit(Y_fine(:,23), 'first');
model_smooth.DUSP = normit(Y_fine(:,24), 'max');
model_smooth.pAKT = normit(Y_fine(:,40), 'max');

% Plotting colors
data_color = [0.8, 0.2, 0.2];    % Red for experimental data
model_color = [0.2, 0.6, 0.2];   % Green for model fit

% Create figure for all species
figure('Name', 'Model Fit - All Species', 'Position', [50, 50, 1600, 1000]);

species_to_plot = {'pEGFR', 'RAS_GTP', 'pCRAF', 'pMEK', 'pERK', 'DUSP', 'pAKT'};
for i = 1:length(species_to_plot)
    subplot(3, 3, i);
    species = species_to_plot{i};
    
    % Plot model fit
    plot(tFine_hours, model_smooth.(species), '-', 'Color', model_color, ...
         'LineWidth', 3, 'DisplayName', 'Model Fit');
    hold on;
    
    % Plot experimental data
    plot(timeStamps_hours, expData_norm.(species), 'o', 'Color', data_color, ...
         'MarkerSize', 10, 'MarkerFaceColor', data_color, 'LineWidth', 2, ...
         'DisplayName', 'Experimental');
    
    xlabel('Time (hours)', 'FontSize', 12);
    ylabel(species, 'FontSize', 12);
    title(sprintf('%s Model Fit', species), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    grid on;
    xlim([0, 50]);
    set(gca, 'FontSize', 11);
    hold off;
end

% Fit error summary plot
subplot(3, 3, 8);
total_errors = [sum(fit_errors.pEGFR), sum(fit_errors.RAS), sum(fit_errors.pCRAF), ...
                sum(fit_errors.pMEK), sum(fit_errors.pERK), sum(fit_errors.DUSP), ...
                sum(fit_errors.pAKT)];
protein_names = {'pEGFR', 'RAS', 'pCRAF', 'pMEK', 'pERK', 'DUSP', 'pAKT'};
bar(1:length(protein_names), total_errors, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k');
set(gca, 'XTick', 1:length(protein_names), 'XTickLabel', protein_names);
ylabel('Total Fit Error', 'FontSize', 12);
title('Fit Errors by Protein', 'FontSize', 14, 'FontWeight', 'bold');
xtickangle(45);
grid on;
set(gca, 'FontSize', 11);

sgtitle('MAPK/PI3K Pathway Model Fit Results', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('Visualization complete.\n\n');

%% ============================================================================
% SECTION 8: ODE FUNCTION
% ============================================================================

function dydt = Mapk_ODE(t, y, p)
    % MAPK/PI3K Pathway ODE System with Vemurafenib and Paradoxical Activation
    %
    % Inputs:
    %   t - time (seconds)
    %   y - state vector (49 species)
    %   p - parameter vector (52 parameters)
    %
    % Outputs:
    %   dydt - derivatives of all species
    %
    % State Vector (y):
    %   [1-3]   EGFR module: EGFR, EGFR:ligand, pEGFR
    %   [4-6]   Shc module: Shc, pEGFR:Shc, pShc
    %   [7-8]   Grb2/SOS: Grb2, pShc:Grb2:SOS
    %   [9-10]  HRAS: HRAS-GDP, HRAS-GTP
    %   [11-12] NRAS: NRAS-GDP, NRAS-GTP
    %   [13-15] KRAS: KRAS-GDP, intermediate, KRAS-GTP
    %   [16-17] CRAF: CRAF, pCRAF
    %   [18-19] BRAF: BRAF, BRAF^P (BRAF^V600E*)
    %   [20-21] MEK: MEK, pMEK
    %   [22-23] ERK: ERK, pERK
    %   [24-25] DUSP: mDUSP, DUSP
    %   [26-27] SPRY: mSPRY, SPRY
    %   [28-31] Degradation trackers
    %   [32-34] IGFR module
    %   [35-36] PI3K: PI3K, PI3K_active
    %   [37-38] PIP: PIP2, PIP3
    %   [39-40] AKT: AKT, pAKT
    %   [41-43] FOXO, mTORC
    %   [44-46] 4EBP1 module
    %   [47-48] KSR: KSR, pKSR
    %   [49]    BRAF-WT:CRAF dimer
    
    % Unpack parameters
    ka1 = p(1);  kr1 = p(2);  kc1 = p(3);
    kpCraf = p(4); kpMek = p(5); kpErk = p(6);
    kDegradEgfr = p(7); kErkInbEgfr = p(8); kShcDephos = p(9); kptpDeg = p(10); kGrb2CombShc = p(11);
    kSprtyInbGrb2 = p(12); kSosCombGrb2 = p(13); kErkPhosSos = p(14);
    kErkPhosPcraf = p(15); kPcrafDegrad = p(16); kErkPhosMek = p(17); kMekDegrad = p(18);
    kDuspInbErk = p(19); kErkDeg = p(20); kinbBraf = p(21); kDuspStop = p(22); kDusps = p(23);
    kSproutyForm = p(24); kSprtyComeDown = p(25); kdegrad = p(26); km_Sprty_decay = p(27); km_Dusp = p(28); km_Sprty = p(29);
    kErkDephos = p(30); kDuspDeg = p(31);
    kEGFRPI3k = p(32); kMTOR_Feedback = p(33);
    kAkt = p(34); kdegradAKT = p(35);
    kb1 = p(36); k43b1 = p(37);
    kKSRphos = p(38); kKSRdephos = p(39);
    kMekByBraf = p(40); kMekByCraf = p(41); kMekByKSR = p(42);
    Tram = p(43); K_tram_RAF = p(44); K_tram_KSR = p(45); n_tram = p(46);
    Vemurafenib = p(47); kDimerForm = p(48); kDimerDissoc = p(49);
    kParadoxCRAF = p(50); IC50_vem = p(51); Hill_n_vem = p(52);
    
    % Initialize derivatives
    dydt = zeros(49, 1);
    
    % ========================================================================
    % MODULE 1: RTK/EGFR SIGNALING
    % ========================================================================
    dydt(1) = -ka1*y(1) + kr1*y(2);                                    % EGFR
    dydt(2) = ka1*y(1) - kr1*y(2) - kc1*y(2);                          % EGFR:ligand
    dydt(3) = kc1*y(2) - kDegradEgfr*y(3) - kErkInbEgfr*y(23)*y(3);    % pEGFR (with ERK feedback)
    
    % ========================================================================
    % MODULE 2: Shc/Grb2/SOS SIGNALING
    % ========================================================================
    dydt(4) = -ka1*y(3)*y(4);                                          % Shc
    dydt(5) = ka1*y(3)*y(4) - kShcDephos*y(6)*y(5);                    % pEGFR:Shc
    dydt(6) = -kptpDeg*y(5)*y(6);                                      % pShc
    
    dydt(7) = kGrb2CombShc*y(5)*y(3) - kSprtyInbGrb2*y(21)*y(7);       % Grb2 (with SPRY feedback)
    dydt(8) = kSosCombGrb2*y(7)*y(5) - kErkPhosSos*y(19)*y(8);         % pShc:Grb2:SOS (with ERK feedback)
    
    % ========================================================================
    % MODULE 3: RAS ACTIVATION
    % ========================================================================
    dydt(9) = -ka1*y(8)*y(9);                                          % HRAS-GDP
    dydt(10) = ka1*y(8)*y(9);                                          % HRAS-GTP
    dydt(11) = -ka1*y(8)*y(11);                                        % NRAS-GDP
    dydt(12) = ka1*y(8)*y(11);                                         % NRAS-GTP
    dydt(13) = -ka1*y(8)*y(13);                                        % KRAS-GDP
    dydt(14) = ka1*y(8)*y(13) - ka1*y(14)*y(15);                       % KRAS-GTP
    dydt(15) = -ka1*y(14)*y(15);                                       % KRAS intermediate
    
    % ========================================================================
    % MODULE 4: RAF SIGNALING WITH PARADOXICAL ACTIVATION
    % ========================================================================
    
    % Vemurafenib inhibition of BRAF^V600E (Hill equation)
    IC50_n = IC50_vem^Hill_n_vem;
    Vem_n = Vemurafenib^Hill_n_vem;
    kBRAF_eff = ka1 * IC50_n / (IC50_n + Vem_n + eps);  % Effective BRAF activation (inhibited)
    
    % Paradoxical activation: vemurafenib promotes BRAF-WT*:CRAF dimer → CRAF* activation
    paradox_activation = kParadoxCRAF * Vemurafenib * y(49);  % y(49) = BRAF-WT:CRAF dimer
    
    % CRAF dynamics
    dydt(16) = -kpCraf*y(14)*y(16) + kErkPhosPcraf*y(23)*y(17) + kPcrafDegrad*y(17)*y(30) ...
               - kDimerForm*y(19)*y(16)*Vemurafenib + kDimerDissoc*y(49);  % CRAF (used in dimer)
    dydt(17) = kpCraf*y(14)*y(16) - kErkPhosPcraf*y(23)*y(17) - kPcrafDegrad*y(17)*y(30) ...
               + paradox_activation;  % pCRAF (with paradoxical activation)
    
    % BRAF dynamics (with vemurafenib inhibition)
    dydt(18) = -kBRAF_eff*y(18) - kDimerForm*y(19)*y(16)*Vemurafenib + kDimerDissoc*y(49);
    dydt(19) = kBRAF_eff*y(18) - kinbBraf*y(19) - kDimerForm*y(19)*y(16)*Vemurafenib + kDimerDissoc*y(49);
    
    % BRAF-WT:CRAF dimer dynamics (y(49) = dimer complex)
    dydt(49) = kDimerForm*y(19)*y(16)*Vemurafenib - kDimerDissoc*y(49) - kPcrafDegrad*y(49)*y(30);
    
    % ========================================================================
    % MODULE 5: MEK PHOSPHORYLATION
    % ========================================================================
    % MEK can be phosphorylated by CRAF, BRAF, or KSR
    raf_to_mek = (kpMek*y(17) + kMekByBraf*y(19) + kMekByCraf*y(17));
    ksr_to_mek = (kMekByKSR * y(48));
    
    dydt(20) = -(raf_to_mek + ksr_to_mek)*y(20) + kErkPhosMek*y(23)*y(21) + kMekDegrad*y(21)*y(29);
    dydt(21) = (raf_to_mek + ksr_to_mek)*y(20) - kErkPhosMek*y(23)*y(21) - kMekDegrad*y(21)*y(29);
    
    % ========================================================================
    % MODULE 6: ERK PHOSPHORYLATION
    % ========================================================================
    dydt(22) = -kpErk*y(21)*y(22) + kDuspInbErk*y(24)*y(23) + kErkDeg*y(23)*y(28) + kErkDephos*y(24)*y(23);
    dydt(23) = kpErk*y(21)*y(22) - kDuspInbErk*y(24)*y(23) - kErkDeg*y(23)*y(28) - kErkDephos*y(24)*y(23);
    
    % ========================================================================
    % MODULE 7: DUSP FEEDBACK
    % ========================================================================
    dydt(24) = km_Dusp*y(23)/(1 + (km_Dusp/kDusps)*y(23)) - kDuspStop*y(24)*y(31) - kDuspDeg*y(24)*y(23);
    dydt(25) = -kDuspStop*y(24)*y(25);
    
    % ========================================================================
    % MODULE 8: SPRY FEEDBACK
    % ========================================================================
    dydt(26) = km_Sprty*y(23)/(1 + (km_Sprty/kSproutyForm)*y(23)) - kSprtyComeDown*y(26)*y(27);
    dydt(27) = -kSprtyComeDown*y(26)*y(27);
    
    % ========================================================================
    % MODULE 9: DEGRADATION TRACKERS
    % ========================================================================
    dydt(28) = -kErkDeg*y(23)*y(28);      % pERK degradation tracker
    dydt(29) = -kMekDegrad*y(21)*y(29);   % pMEK degradation tracker
    dydt(30) = -kPcrafDegrad*y(17)*y(30); % pCRAF degradation tracker
    dydt(31) = -kDuspStop*y(24)*y(31);    % DUSP stop tracker
    
    % ========================================================================
    % MODULE 10: PI3K/AKT/mTOR SIGNALING
    % ========================================================================
    % IGFR module
    dydt(32) = -ka1*y(32) + kr1*y(33);
    dydt(33) = ka1*y(32) - kr1*y(33) - kc1*y(33);
    dydt(34) = kc1*y(33) - kErkInbEgfr*y(23)*y(34);
    
    % PI3K module
    dydt(35) = -ka1*y(3)*y(35) - kEGFRPI3k*y(34)*y(35);
    dydt(36) = ka1*y(3)*y(35) + kEGFRPI3k*y(34)*y(35) - kMTOR_Feedback*y(43)*y(36);
    
    % PIP module
    dydt(37) = -(ka1*y(36)*y(37) + ka1*y(14)*y(37));
    dydt(38) = ka1*y(36)*y(37) + ka1*y(14)*y(37) - kdegrad*y(38);
    
    % AKT module
    dydt(39) = -kAkt*y(38)*y(39);
    dydt(40) = kAkt*y(38)*y(39) - kdegradAKT*y(40);
    
    % FOXO module
    dydt(41) = (1-y(40))*ka1/(1 + (y(41)/15e-5));
    
    % mTORC module
    dydt(42) = -ka1*y(40)*y(42) + kdegrad*y(43);
    dydt(43) = ka1*y(40)*y(42) - kdegrad*y(43);
    
    % 4EBP1 module
    dydt(44) = -ka1*y(43)*y(44) + kb1*y(44);
    dydt(45) = ka1*y(43)*y(44) - kb1*y(44) - k43b1*y(45);
    dydt(46) = k43b1*y(45);
    
    % ========================================================================
    % MODULE 11: KSR SCAFFOLD
    % ========================================================================
    dydt(47) = -kKSRphos*(y(17)+y(19))*y(47) + kKSRdephos*y(48);
    dydt(48) = kKSRphos*(y(17)+y(19))*y(47) - kKSRdephos*y(48);
    
    % Note: dydt(49) for BRAF-WT:CRAF dimer is defined above in RAF module
end

%% ============================================================================
% SECTION 9: OBJECTIVE FUNCTION FOR OPTIMIZATION
% ============================================================================

function err = objectiveFunction_all(p, timeStamps_seconds, expData_norm, y0)
    % Objective function for parameter optimization
    %
    % Inputs:
    %   p - parameter vector (52 parameters)
    %   timeStamps_seconds - experimental time points (seconds)
    %   expData_norm - normalized experimental data (structure)
    %   y0 - initial conditions
    %
    % Outputs:
    %   err - sum of squared residuals
    
    % Integrate ODE system at experimental time points
    [T, Y] = ode45(@(t,y) Mapk_ODE(t, y, p), timeStamps_seconds, y0);
    
    % Extract model outputs
    m_pEGFR = Y(:,3);
    m_pCRAF = Y(:,17);
    m_pMEK = Y(:,21);
    m_pERK = Y(:,23);
    m_DUSP = Y(:,24);
    m_pAKT = Y(:,40);
    m_RAS_GTP = Y(:,14);
    
    % Normalization function
    norm = @(v, mode) v ./ max(eps, (mode=="first")*v(1) + (mode=="last")*v(end) + (mode=="max")*max(v));
    
    % Normalize model outputs
    m_pEGFR = norm(m_pEGFR, 'max');
    m_pMEK = norm(m_pMEK, 'first');
    m_pERK = norm(m_pERK, 'first');
    m_DUSP = norm(m_DUSP, 'max');
    m_pCRAF = norm(m_pCRAF, 'max');
    m_pAKT = norm(m_pAKT, 'max');
    m_RAS_GTP = norm(m_RAS_GTP, 'last');
    
    % Calculate weighted sum of squared residuals
    w = struct('EGFR', 1, 'MEK', 1, 'ERK', 3, 'DUSP', 3, 'CRAF', 1, 'AKT', 1, 'RAS', 1);
    
    err = 0;
    err = err + w.EGFR * sum((m_pEGFR - expData_norm.pEGFR(:)).^2);
    err = err + w.MEK * sum((m_pMEK - expData_norm.pMEK(:)).^2);
    err = err + w.ERK * sum((m_pERK - expData_norm.pERK(:)).^2);
    err = err + w.DUSP * sum((m_DUSP - expData_norm.DUSP(:)).^2);
    err = err + w.CRAF * sum((m_pCRAF - expData_norm.pCRAF(:)).^2);
    err = err + w.AKT * sum((m_pAKT - expData_norm.pAKT(:)).^2);
    err = err + w.RAS * sum((m_RAS_GTP - expData_norm.RAS_GTP(:)).^2);
end

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   SIMULATION COMPLETED SUCCESSFULLY!\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');
