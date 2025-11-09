% =============================================================================
% RAS-RAF-MEK-ERK PATHWAY WITH BRAF^V600E MUTATION AND pERK FEEDBACK
% MAPK Signaling with Constitutive BRAF Mutant and Negative Feedback
% =============================================================================
%
% Model Description:
% This script models the MAPK signaling cascade including:
% - RAS activation by SOS*
% - Wild-type RAF activation (RAS-dependent)
% - BRAF^V600E constitutive activity (RAS-independent)
% - MEK and ERK sequential phosphorylation
% - pERK (ERK-P) negative feedback on upstream components
%
% Key Features:
% - BRAF^V600E is constitutively active and bypasses RAS requirement
% - Both RAF* and BRAF^V600E* phosphorylate MEK
% - pERK feedback affects RAS-dependent pathway more than BRAF^V600E
%
% =============================================================================

clear all;
close all;
clc;

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   RAS-RAF-MEK-ERK PATHWAY WITH BRAF^V600E MUTATION\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

%% ============================================================================
% INPUT SIGNAL: SOS*(t)
% ============================================================================

% Define SOS* as a function of time
% Option 1: Transient pulse (exponential decay)
SOS_star_func = @(t) 1.0 * exp(-0.05 * t) .* (t >= 0);

% Option 2: Step function (constant)
% SOS_star_func = @(t) 0.8 * (t >= 0);

% Option 3: Load from previous simulation
% if exist('sos_output.mat', 'file')
%     load('sos_output.mat', 't', 'SOS_star');
%     t_sos = t;
%     SOS_star_data = SOS_star;
%     SOS_star_func = @(t) interp1(t_sos, SOS_star_data, t, 'linear', SOS_star_data(end));
% end

fprintf('Using SOS*(t) as input signal.\n\n');

%% ============================================================================
% PARAMETERS
% ============================================================================

fprintf('Setting up parameters...\n');

% ============================================================================
% RAS Module Parameters
% ============================================================================

% Reaction 1: RAS-GDP + SOS* → RAS-GTP + SOS*
k_SOS = 0.5;        % min^-1 - RAS activation rate by SOS* (base rate)

% Reaction 2: RAS-GTP → RAS-GDP (GTP hydrolysis)
k_GTPase = 0.1;     % min^-1 - RAS intrinsic GTP hydrolysis rate (base rate)

% ============================================================================
% RAF Module Parameters
% ============================================================================

% Reaction 3: RAS-GTP + RAF → RAF* (Wild-type RAF activation)
k_RAF_act = 0.3;    % min^-1 - RAF activation rate by RAS-GTP (base rate)

% Reaction 4: RAF* → RAF (dephosphorylation)
k_RAF_deact = 0.05; % min^-1 - RAF* deactivation rate

% Reaction 5: BRAF^V600E → BRAF^V600E* (constitutive activation)
k_BRAF_mut = 0.5;   % min^-1 - BRAF^V600E constitutive activation rate

% Reaction 6: BRAF^V600E* → BRAF^V600E (deactivation, slower than wild-type)
k_BRAF_mut_deact = 0.02; % min^-1 - BRAF^V600E* deactivation rate (lower than wild-type)

% ============================================================================
% MEK Module Parameters
% ============================================================================

% Reaction 7: RAF* + MEK → RAF* + MEK-P
k_MEK_wt = 0.5;     % min^-1 - MEK phosphorylation rate by RAF*

% Reaction 8: BRAF^V600E* + MEK → BRAF^V600E* + MEK-P
k_MEK_mut = 0.5;    % min^-1 - MEK phosphorylation rate by BRAF^V600E*

% Reaction 9: MEK-P → MEK (dephosphorylation)
k_MEK_deact = 0.05; % min^-1 - MEK-P dephosphorylation rate

% ============================================================================
% ERK Module Parameters
% ============================================================================

% Reaction 10: MEK-P + ERK → MEK-P + ERK-P
k_ERK = 0.5;        % min^-1 - ERK phosphorylation rate by MEK-P

% Reaction 11: ERK-P → ERK (dephosphorylation)
k_ERK_deact = 0.05; % min^-1 - ERK-P dephosphorylation rate

% ============================================================================
% Negative Feedback Parameters
% ============================================================================

alpha1 = 0.5;       % Feedback strength: ERK-P inhibition of Shc binding (upstream)
alpha2 = 0.5;       % Feedback strength: ERK-P inhibition of SOS activation
alpha3 = 0.5;       % Feedback strength: ERK-P inhibition of RAF activation
alpha4 = 0.2;       % Feedback strength: ERK-P enhancement of RAS GTPase

% Feedback on BRAF^V600E (reduced, as it's RAS-independent)
alpha_BRAF_mut = 0.1; % Feedback strength: ERK-P on BRAF^V600E (much weaker)

% Store parameters in structure
params.k_SOS = k_SOS;
params.k_GTPase = k_GTPase;
params.k_RAF_act = k_RAF_act;
params.k_RAF_deact = k_RAF_deact;
params.k_BRAF_mut = k_BRAF_mut;
params.k_BRAF_mut_deact = k_BRAF_mut_deact;
params.k_MEK_wt = k_MEK_wt;
params.k_MEK_mut = k_MEK_mut;
params.k_MEK_deact = k_MEK_deact;
params.k_ERK = k_ERK;
params.k_ERK_deact = k_ERK_deact;
params.alpha1 = alpha1;
params.alpha2 = alpha2;
params.alpha3 = alpha3;
params.alpha4 = alpha4;
params.alpha_BRAF_mut = alpha_BRAF_mut;
params.SOS_star_func = SOS_star_func;

fprintf('Parameters configured.\n\n');

%% ============================================================================
% INITIAL CONDITIONS
% ============================================================================

fprintf('Setting initial conditions...\n');

% RAS species
RAS_GDP0 = 10.0;    % nM - inactive RAS (GDP-bound)
RAS_GTP0 = 0.0;     % nM - active RAS (GTP-bound)

% Wild-type RAF species
RAF0 = 5.0;         % nM - inactive wild-type RAF
RAF_star0 = 0.0;    % nM - active wild-type RAF

% BRAF^V600E mutant species
BRAF_mut0 = 1.0;    % nM - inactive BRAF^V600E (constitutive mutant)
BRAF_mut_star0 = 0.0; % nM - active BRAF^V600E (initially zero, but activates quickly)

% MEK species
MEK0 = 5.0;         % nM - inactive MEK
MEK_P0 = 0.0;       % nM - phosphorylated MEK

% ERK species
ERK0 = 5.0;         % nM - inactive ERK
ERK_P0 = 0.0;       % nM - phosphorylated ERK (pERK)

% Initial state vector: [RAS_GDP, RAS_GTP, RAF, RAF_star, BRAF_mut, 
%                        BRAF_mut_star, MEK, MEK_P, ERK, ERK_P]
y0 = [RAS_GDP0; RAS_GTP0; RAF0; RAF_star0; BRAF_mut0; BRAF_mut_star0; ...
      MEK0; MEK_P0; ERK0; ERK_P0];

fprintf('Initial conditions set.\n\n');

%% ============================================================================
% TIME SPAN
% ============================================================================

tspan = [0, 200];  % minutes

%% ============================================================================
% SOLVE ODE SYSTEM
% ============================================================================

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   SOLVING ODE SYSTEM\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

fprintf('Solving ODE system (this may take a moment)...\n');
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, y] = ode45(@(t,y) braf_mutant_odes(t, y, params), tspan, y0, options);

% Extract species concentrations
RAS_GDP = y(:, 1);
RAS_GTP = y(:, 2);
RAF = y(:, 3);
RAF_star = y(:, 4);
BRAF_mut = y(:, 5);
BRAF_mut_star = y(:, 6);
MEK = y(:, 7);
MEK_P = y(:, 8);
ERK = y(:, 9);
ERK_P = y(:, 10);  % pERK - feedback signal

% Calculate SOS_star input signal over time
SOS_star_signal = arrayfun(SOS_star_func, t);

fprintf('Simulation completed successfully!\n\n');

%% ============================================================================
% PLOTTING
% ============================================================================

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   GENERATING VISUALIZATIONS\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

% Figure 1: Key Species Time Courses
figure('Position', [100, 100, 1400, 900]);

% Plot 1: RAS species
subplot(2, 3, 1);
plot(t, RAS_GDP, 'b-', 'LineWidth', 2, 'DisplayName', 'RAS-GDP');
hold on;
plot(t, RAS_GTP, 'r-', 'LineWidth', 2, 'DisplayName', 'RAS-GTP');
xlabel('Time (minutes)', 'FontSize', 11);
ylabel('Concentration (nM)', 'FontSize', 11);
title('RAS Species', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, 200]);
hold off;

% Plot 2: Wild-type RAF species
subplot(2, 3, 2);
plot(t, RAF, 'b-', 'LineWidth', 2, 'DisplayName', 'RAF (WT)');
hold on;
plot(t, RAF_star, 'r-', 'LineWidth', 2, 'DisplayName', 'RAF* (WT)');
xlabel('Time (minutes)', 'FontSize', 11);
ylabel('Concentration (nM)', 'FontSize', 11);
title('Wild-type RAF', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, 200]);
hold off;

% Plot 3: BRAF^V600E mutant species
subplot(2, 3, 3);
plot(t, BRAF_mut, 'b-', 'LineWidth', 2, 'DisplayName', 'BRAF^{V600E}');
hold on;
plot(t, BRAF_mut_star, 'r-', 'LineWidth', 2, 'DisplayName', 'BRAF^{V600E}*');
xlabel('Time (minutes)', 'FontSize', 11);
ylabel('Concentration (nM)', 'FontSize', 11);
title('BRAF^{V600E} Mutant', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, 200]);
hold off;

% Plot 4: MEK species
subplot(2, 3, 4);
plot(t, MEK, 'b-', 'LineWidth', 2, 'DisplayName', 'MEK');
hold on;
plot(t, MEK_P, 'r-', 'LineWidth', 2, 'DisplayName', 'MEK-P');
xlabel('Time (minutes)', 'FontSize', 11);
ylabel('Concentration (nM)', 'FontSize', 11);
title('MEK Species', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, 200]);
hold off;

% Plot 5: ERK species
subplot(2, 3, 5);
plot(t, ERK, 'b-', 'LineWidth', 2, 'DisplayName', 'ERK');
hold on;
plot(t, ERK_P, 'r-', 'LineWidth', 2.5, 'DisplayName', 'ERK-P (pERK)');
xlabel('Time (minutes)', 'FontSize', 11);
ylabel('Concentration (nM)', 'FontSize', 11);
title('ERK Species (Feedback Signal)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, 200]);
hold off;

% Plot 6: Input signal
subplot(2, 3, 6);
plot(t, SOS_star_signal, 'k-', 'LineWidth', 2);
xlabel('Time (minutes)', 'FontSize', 11);
ylabel('Concentration (nM)', 'FontSize', 11);
title('Input Signal: SOS*', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
xlim([0, 200]);

sgtitle('RAS-RAF-MEK-ERK Pathway with BRAF^{V600E} Mutation', ...
        'FontSize', 15, 'FontWeight', 'bold');

% Figure 2: Comparison of Wild-type vs Mutant Pathways
figure('Position', [150, 150, 1400, 800]);

% Plot 1: RAF activation comparison
subplot(2, 2, 1);
plot(t, RAF_star, 'b-', 'LineWidth', 2, 'DisplayName', 'RAF* (WT, RAS-dependent)');
hold on;
plot(t, BRAF_mut_star, 'r-', 'LineWidth', 2, 'DisplayName', 'BRAF^{V600E}* (constitutive)');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Concentration (nM)', 'FontSize', 12);
title('RAF Activation: Wild-type vs Mutant', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, 200]);
hold off;

% Plot 2: Cascade progression
subplot(2, 2, 2);
plot(t, RAS_GTP, 'c-', 'LineWidth', 2, 'DisplayName', 'RAS-GTP');
hold on;
plot(t, RAF_star, 'b-', 'LineWidth', 2, 'DisplayName', 'RAF* (WT)');
plot(t, BRAF_mut_star, 'r-', 'LineWidth', 2, 'DisplayName', 'BRAF^{V600E}*');
plot(t, MEK_P, 'm-', 'LineWidth', 2, 'DisplayName', 'MEK-P');
plot(t, ERK_P, 'k-', 'LineWidth', 2.5, 'DisplayName', 'ERK-P');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Concentration (nM)', 'FontSize', 12);
title('Complete Cascade Progression', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, 200]);
hold off;

% Plot 3: Feedback effects on rate constants
subplot(2, 2, 3);
% Calculate effective rates (with feedback) over time
k_SOS_eff = k_SOS ./ (1 + alpha2 * ERK_P);
k_RAF_act_eff = k_RAF_act ./ (1 + alpha3 * ERK_P);
k_GTPase_eff = k_GTPase + alpha4 * ERK_P;
k_BRAF_mut_eff = k_BRAF_mut ./ (1 + alpha_BRAF_mut * ERK_P);

plot(t, k_SOS_eff, 'g-', 'LineWidth', 2, 'DisplayName', 'k_{SOS} (RAS activation)');
hold on;
plot(t, k_RAF_act_eff, 'b-', 'LineWidth', 2, 'DisplayName', 'k_{RAF,act} (WT RAF)');
plot(t, k_GTPase_eff, 'c-', 'LineWidth', 2, 'DisplayName', 'k_{GTPase}');
plot(t, k_BRAF_mut_eff, 'r-', 'LineWidth', 2, 'DisplayName', 'k_{BRAF^{V600E}}');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Effective Rate Constant', 'FontSize', 12);
title('Feedback Effects on Rate Constants', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, 200]);
hold off;

% Plot 4: Contribution to MEK-P from each RAF species
subplot(2, 2, 4);
% Calculate contributions (approximate from rates)
contribution_wt = k_MEK_wt * RAF_star .* MEK;
contribution_mut = k_MEK_mut * BRAF_mut_star .* MEK;

plot(t, contribution_wt, 'b-', 'LineWidth', 2, 'DisplayName', 'Contribution from RAF* (WT)');
hold on;
plot(t, contribution_mut, 'r-', 'LineWidth', 2, 'DisplayName', 'Contribution from BRAF^{V600E}*');
plot(t, MEK_P, 'k--', 'LineWidth', 2, 'DisplayName', 'Total MEK-P');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Rate / Concentration (nM/min)', 'FontSize', 12);
title('MEK-P Activation: Wild-type vs Mutant Contribution', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, 200]);
hold off;

sgtitle('Wild-type vs BRAF^{V600E} Pathway Comparison', 'FontSize', 16, 'FontWeight', 'bold');

%% ============================================================================
% SUMMARY STATISTICS
% ============================================================================

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   SUMMARY STATISTICS\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

fprintf('KEY SPECIES PEAK VALUES:\n');
[max_RAS_GTP, idx_RAS] = max(RAS_GTP);
[max_RAF_star, idx_RAF] = max(RAF_star);
[max_BRAF_mut_star, idx_BRAF] = max(BRAF_mut_star);
[max_MEK_P, idx_MEK] = max(MEK_P);
[max_ERK_P, idx_ERK] = max(ERK_P);

fprintf('  RAS-GTP:            %.4f nM at t = %.2f min\n', max_RAS_GTP, t(idx_RAS));
fprintf('  RAF* (WT):           %.4f nM at t = %.2f min\n', max_RAF_star, t(idx_RAF));
fprintf('  BRAF^{V600E}*:       %.4f nM at t = %.2f min\n', max_BRAF_mut_star, t(idx_BRAF));
fprintf('  MEK-P:              %.4f nM at t = %.2f min\n', max_MEK_P, t(idx_MEK));
fprintf('  ERK-P (pERK):       %.4f nM at t = %.2f min\n', max_ERK_P, t(idx_ERK));

fprintf('\nSTEADY-STATE VALUES (at t = 200 min):\n');
fprintf('  RAS-GTP:            %.4f nM\n', RAS_GTP(end));
fprintf('  RAF* (WT):          %.4f nM\n', RAF_star(end));
fprintf('  BRAF^{V600E}*:      %.4f nM\n', BRAF_mut_star(end));
fprintf('  MEK-P:              %.4f nM\n', MEK_P(end));
fprintf('  ERK-P (pERK):       %.4f nM\n', ERK_P(end));

fprintf('\nPATHWAY COMPARISON:\n');
fprintf('  Peak RAF* (WT) / Peak BRAF^{V600E}*: %.2f\n', max_RAF_star / max_BRAF_mut_star);
fprintf('  Steady-state RAF* (WT) / BRAF^{V600E}*: %.2f\n', RAF_star(end) / BRAF_mut_star(end));

fprintf('\nFEEDBACK STRENGTHS:\n');
fprintf('  alpha1 (Shc binding):     %.2f\n', alpha1);
fprintf('  alpha2 (SOS activation):  %.2f\n', alpha2);
fprintf('  alpha3 (RAF activation):  %.2f\n', alpha3);
fprintf('  alpha4 (RAS GTPase):      %.2f\n', alpha4);
fprintf('  alpha_BRAF_mut:           %.2f (reduced feedback on mutant)\n', alpha_BRAF_mut);

fprintf('\n═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   SIMULATION COMPLETED SUCCESSFULLY!\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

%% ============================================================================
% SAVE OUTPUTS
% ============================================================================

save('braf_mutant_output.mat', 't', 'RAS_GTP', 'RAF_star', 'BRAF_mut_star', ...
     'MEK_P', 'ERK_P', 'SOS_star_signal');
fprintf('Saved outputs to braf_mutant_output.mat\n');

%% ============================================================================
% ODE FUNCTION
% ============================================================================

function dydt = braf_mutant_odes(t, y, p)
    % ODE system for RAS-RAF-MEK-ERK pathway with BRAF^V600E mutation and pERK feedback
    % y = [RAS_GDP, RAS_GTP, RAF, RAF_star, BRAF_mut, BRAF_mut_star, 
    %      MEK, MEK_P, ERK, ERK_P]
    % p = parameters structure
    
    % Extract species
    RAS_GDP = y(1);
    RAS_GTP = y(2);
    RAF = y(3);
    RAF_star = y(4);
    BRAF_mut = y(5);
    BRAF_mut_star = y(6);
    MEK = y(7);
    MEK_P = y(8);
    ERK = y(9);
    ERK_P = y(10);  % pERK - feedback signal
    
    % Get SOS_star input signal at current time
    SOS_star = p.SOS_star_func(t);
    
    % ========================================================================
    % CALCULATE FEEDBACK-MODULATED RATE CONSTANTS
    % ========================================================================
    
    % Feedback 1: ERK-P reduces SOS activation (affects RAS activation)
    k_SOS_eff = p.k_SOS / (1 + p.alpha2 * ERK_P);
    
    % Feedback 2: ERK-P increases RAS GTPase activity
    k_GTPase_eff = p.k_GTPase + p.alpha4 * ERK_P;
    
    % Feedback 3: ERK-P reduces wild-type RAF activation
    k_RAF_act_eff = p.k_RAF_act / (1 + p.alpha3 * ERK_P);
    
    % Feedback 4: ERK-P has reduced effect on BRAF^V600E (RAS-independent)
    k_BRAF_mut_eff = p.k_BRAF_mut / (1 + p.alpha_BRAF_mut * ERK_P);
    
    % ========================================================================
    % RAS MODULE REACTIONS
    % ========================================================================
    
    % Reaction 1: RAS-GDP + SOS* → RAS-GTP + SOS*
    r1 = k_SOS_eff * SOS_star * RAS_GDP;  % Feedback-modulated
    
    % Reaction 2: RAS-GTP → RAS-GDP (GTP hydrolysis)
    r2 = k_GTPase_eff * RAS_GTP;  % Feedback-modulated
    
    % ========================================================================
    % RAF MODULE REACTIONS
    % ========================================================================
    
    % Reaction 3: RAS-GTP + RAF → RAF* (Wild-type RAF activation)
    r3 = k_RAF_act_eff * RAS_GTP * RAF;  % Feedback-modulated, RAS-dependent
    
    % Reaction 4: RAF* → RAF (dephosphorylation)
    r4 = p.k_RAF_deact * RAF_star;
    
    % Reaction 5: BRAF^V600E → BRAF^V600E* (constitutive activation)
    r5 = k_BRAF_mut_eff * BRAF_mut;  % Feedback-modulated (weakly), RAS-independent
    
    % Reaction 6: BRAF^V600E* → BRAF^V600E (deactivation)
    r6 = p.k_BRAF_mut_deact * BRAF_mut_star;
    
    % ========================================================================
    % MEK MODULE REACTIONS
    % ========================================================================
    
    % Reaction 7: RAF* + MEK → RAF* + MEK-P (Wild-type pathway)
    r7 = p.k_MEK_wt * RAF_star * MEK;
    
    % Reaction 8: BRAF^V600E* + MEK → BRAF^V600E* + MEK-P (Mutant pathway)
    r8 = p.k_MEK_mut * BRAF_mut_star * MEK;
    
    % Reaction 9: MEK-P → MEK (dephosphorylation)
    r9 = p.k_MEK_deact * MEK_P;
    
    % ========================================================================
    % ERK MODULE REACTIONS
    % ========================================================================
    
    % Reaction 10: MEK-P + ERK → MEK-P + ERK-P
    r10 = p.k_ERK * MEK_P * ERK;
    
    % Reaction 11: ERK-P → ERK (dephosphorylation)
    r11 = p.k_ERK_deact * ERK_P;
    
    % ========================================================================
    % ODEs FOR EACH SPECIES
    % ========================================================================
    
    % RAS module
    dRAS_GDP_dt = -r1 + r2;
    dRAS_GTP_dt = r1 - r2 - r3;
    
    % Wild-type RAF module
    dRAF_dt = -r3 + r4;
    dRAF_star_dt = r3 - r4;
    
    % BRAF^V600E mutant module
    dBRAF_mut_dt = -r5 + r6;
    dBRAF_mut_star_dt = r5 - r6;
    
    % MEK module (activated by both RAF* and BRAF^V600E*)
    dMEK_dt = -r7 - r8 + r9;
    dMEK_P_dt = r7 + r8 - r9;
    
    % ERK module
    dERK_dt = -r10 + r11;
    dERK_P_dt = r10 - r11;
    
    % Return derivatives
    dydt = [dRAS_GDP_dt; dRAS_GTP_dt; dRAF_dt; dRAF_star_dt; ...
            dBRAF_mut_dt; dBRAF_mut_star_dt; dMEK_dt; dMEK_P_dt; ...
            dERK_dt; dERK_P_dt];
end

