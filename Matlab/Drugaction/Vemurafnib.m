% =============================================================================
% VEMURAFENIB EFFECTS ON MAPK AND PI3K/AKT/mTOR PATHWAYS
% Drug Inhibition of BRAF^V600E with Adaptive Signaling and Cross-Talk
% =============================================================================
%
% IMPORTANT: All concentrations and species values are normalized to [0, 1]
% - Initial conditions: Total pools = 1.0, distributed between active/inactive forms
% - Input signals: pRTK and vemurafenib normalized to [0, 1]
% - IC50: Normalized value (0.1 = 10% of maximum drug concentration)
%
% Model Description:
% This script models the effects of vemurafenib (BRAF inhibitor) on:
% - MAPK pathway with BRAF^V600E mutation
% - PI3K/AKT/mTOR pathway
% - Cross-talk mechanisms and adaptive responses
%
% PATHWAY 1: EGFR → MAPK (RAS-RAF-MEK-ERK)
%   - Shc recruitment and phosphorylation by pEGFR
%   - Grb2 binding and SOS activation
%   - RAS activation by SOS*
%   - Wild-type RAF activation (RAS-dependent, AKT-inhibited)
%   - BRAF^V600E constitutive activation (inhibited by vemurafenib)
%   - Sequential MEK and ERK phosphorylation
%   - ERK-P negative feedback (relieved when vemurafenib reduces ERK-P)
%
% PATHWAY 2: EGFR → PI3K/AKT/mTOR
%   - PI3K activation by pEGFR/pHER3 (and RAS-GTP cross-talk)
%   - PIP2 to PIP3 conversion
%   - AKT dual phosphorylation (Thr308, Ser473)
%   - TSC1/TSC2 inhibition by AKT and ERK-P (cross-talk)
%   - Rheb-GTP accumulation and mTORC1 activation
%   - S6K and 4EBP1 phosphorylation
%   - pS6K negative feedback on PI3K
%
% VEMURAFENIB EFFECTS:
% 1. Direct inhibition: Vemurafenib inhibits BRAF^V600E kinase activity using Hill equation
%    - Hill equation: k_eff = k_base * IC50^n / (IC50^n + [Drug]^n)
%    - IC50: half-maximal inhibitory concentration
%    - Hill coefficient n: cooperativity parameter (typically 1-2)
% 2. Paradoxical activation: Vemurafenib can trans-activate wild-type RAF in dimers
% 3. Feedback relief: Reduced ERK-P relieves negative feedback on SOS/Shc
% 4. Adaptive signaling: RAS-GTP rebound → PI3K activation rebound
% 5. Cross-talk changes: Reduced ERK-P → reduced ERK→TSC2 enhancement
%
% CROSS-TALK MECHANISMS:
% 1. RAS → PI3K: RAS-GTP activates PI3K (rebound after vemurafenib)
% 2. ERK → TSC2/mTORC1: ERK-P phosphorylates TSC2 (reduced after vemurafenib)
% 3. AKT → RAF: AKT inhibits wild-type RAF (not BRAF^V600E)
% 4. pS6K → PI3K: Negative feedback on PI3K activation
% 5. PTEN: Degrades PIP3, limiting AKT activation
%
% =============================================================================

clear all;
close all;
clc;

% Experimental values normalized with Vinculin, Experimental time is hours

EGFR = [0.123025235,0.107275504,0.269197631,0.307347132,0.415349555,0.330053274]
pEGFR = [0.291928893,0.392400458,0.265016688,0.394238749,0.006158316,0.008115099]
her3 = [0.203233765,0.194358998,0.303475212,0.674083831,0.89702403,0.459831389]
her2 = [0.245236744,0.177917339,0.239075259,0.306884773,1.066654783,1.005085151]
IGF1R = [0.96024414,1.070624757,1.126551098,1.157072001,0.734617533,0.457922622]
PDGFR = [0.474174188,0.492132953,0.743620725,1.266460499,2.514722273,2.482761079]
PanRAS = [0.839133594,0.833259289,0.919508516,1.240888235,1.582734859,1.468310571]
pCRAF = [0.366397596,0.537106733,0.465541704,0.586732657,1.102322681,0.269181259]
pMEK = [1.75938884,0.170160085,0.095112609,0.201000276,0.219207054,0.502831668]
pERK = [2.903209735,0.207867788,0.303586121,0.805254439,1.408362153,1.847606441]
DUSP6 = [2.677161325,2.782754577,1.130758062,0.395642757,0.828575853,0.916618219]
DUSP4 = [2.035822288,1.560630445,0.235971416,0.3076027,0.143555376,0.120326106]
pAKT308 = [0.407894525,0.552421756,1.005350576,0.860214813,0.444584869,0.339957286]
pAKT473 = [0.513544148,0.613178403,1.03451863,1.113391047,0.535242724,0.538273551]
pS6K = [1.432459522,1.520433646,1.542177411,1.248505245,0.109963216,0.013374136]
p4EBP1 = [1.002468056,1.276793699,1.252681407,1.707504483,1.271216967,0.61389625]
experimentalTime = [0,1,4,8,24,48];  % Time in hours
experimentalTime_min = experimentalTime * 60;  % Convert to minutes

% Normalize experimental values to [0, 1] using min-max normalization
% Function to normalize a vector to [0, 1]
normalize_0_1 = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);  % Add eps to avoid division by zero for constant vectors

% Normalize all experimental datasets
EGFR = normalize_0_1(EGFR);
pEGFR = normalize_0_1(pEGFR);
her3 = normalize_0_1(her3);
her2 = normalize_0_1(her2);
IGF1R = normalize_0_1(IGF1R);
PDGFR = normalize_0_1(PDGFR);
PanRAS = normalize_0_1(PanRAS);
pCRAF = normalize_0_1(pCRAF);
pMEK = normalize_0_1(pMEK);
pERK = normalize_0_1(pERK);
DUSP6 = normalize_0_1(DUSP6);
DUSP4 = normalize_0_1(DUSP4);
pAKT308 = normalize_0_1(pAKT308);
pAKT473 = normalize_0_1(pAKT473);
pS6K = normalize_0_1(pS6K);
p4EBP1 = normalize_0_1(p4EBP1);

fprintf('Experimental data normalized to [0, 1] range.\n');

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   VEMURAFENIB EFFECTS ON MAPK AND PI3K/AKT/mTOR PATHWAYS\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

%% ============================================================================
% INPUT SIGNALS
% ============================================================================

% Define pRTK (pEGFR or pHER3) as a function of time (normalized to [0, 1])
% Option 1: Transient pulse (exponential decay, normalized)
pRTK_func = @(t) 1.0 * exp(-0.05 * t) .* (t >= 0);  % Normalized: starts at 1.0, decays

% Option 2: Step function (constant, normalized)
% pRTK_func = @(t) 1.0 * (t >= 0);  % Normalized: constant at 1.0

% Option 3: Load from previous simulation
% if exist('egfr_output.mat', 'file')
%     load('egfr_output.mat', 't', 'EEp', 'EHp');
%     t_rtk = t;
%     pRTK_data = EEp + 0.5 * EHp;
%     pRTK_func = @(t) interp1(t_rtk, pRTK_data, t, 'linear', pRTK_data(end));
% end

% Define vemurafenib concentration as a function of time (normalized to [0, 1])
% Option 1: Step function (constant drug concentration, normalized)
% Drug is present from t=0 (experimental data includes drug from start)
vemurafenib_func = @(t) 1.0 * (t >= 0);  % Normalized: drug = 1.0 (maximum) from t = 0

% Option 2: Ramp function (gradual increase)
% vemurafenib_func = @(t) 1.0 * max(0, (t - 50) / 10) .* (t >= 50) .* (t <= 60) + 1.0 * (t > 60);

% Option 3: Time-dependent (decay)
% vemurafenib_func = @(t) 1.0 * exp(-0.02 * (t - 50)) .* (t >= 50);

fprintf('Using pRTK (pEGFR/pHER3) as input signal.\n');
fprintf('Using vemurafenib(t) as drug input signal.\n\n');

%% ============================================================================
% PARAMETERS
% ============================================================================

fprintf('Setting up parameters for both pathways with vemurafenib effects...\n');
fprintf('NOTE: All rate constants have been reduced by a factor of 10 to slow down reactions.\n\n');

% ============================================================================
% MAPK PATHWAY PARAMETERS
% ============================================================================

% Shc–Grb2–SOS Module (slowed down by factor of 10)
k_on_Shc = 0.001;      % nM^-1 min^-1 - Shc binding to pEGFR (base rate)
k_off_Shc = 0.001;     % min^-1 - Dissociation of pEGFR:Shc
k_cat1 = 0.05;         % min^-1 - Shc phosphorylation rate
k_cat_SOS = 0.05;      % min^-1 - SOS activation rate (base rate)
k_on2 = 0.005;         % nM^-1 min^-1 - Grb2 binding to pShc
k_off2 = 0.005;        % min^-1 - Dissociation of pShc:Grb2
k_on3 = 0.005;         % nM^-1 min^-1 - SOS binding to pShc:Grb2
k_off3 = 0.005;        % min^-1 - Dissociation of ternary complex
k_deg4 = 0.005;        % min^-1 - SOS* deactivation rate

% RAS–RAF–MEK–ERK Module (slowed down by factor of 10)
k_SOS = 0.05;          % min^-1 - RAS activation rate by SOS*
k_GTPase = 0.01;       % min^-1 - RAS GTP hydrolysis rate (base rate)
k_RAF_act = 0.03;      % min^-1 - Wild-type RAF activation rate (base rate)
k_RAF_deact = 0.005;   % min^-1 - RAF* deactivation rate

% BRAF^V600E Mutation Parameters (slowed down by factor of 10)
k_BRAF_mut = 0.05;     % min^-1 - BRAF^V600E constitutive activation rate
k_BRAF_mut_deact = 0.002; % min^-1 - BRAF^V600E* deactivation rate

% MEK and ERK Module (slowed down by factor of 10)
k_MEK = 0.05;          % min^-1 - MEK phosphorylation rate
k_MEK_deact = 0.005;   % min^-1 - MEK-P dephosphorylation rate
k_ERK = 0.05;          % min^-1 - ERK phosphorylation rate
k_ERK_deact = 0.005;   % min^-1 - ERK-P dephosphorylation rate

% MAPK Feedback Parameters (ERK-P feedback)
alpha_ERK_Shc = 0.5;  % ERK-P inhibition of Shc binding
alpha_ERK_SOS = 0.5;  % ERK-P inhibition of SOS activation
alpha_ERK_RAF = 0.5;  % ERK-P inhibition of wild-type RAF activation
alpha_ERK_BRAF_mut = 0.1; % ERK-P on BRAF^V600E (much weaker feedback)
alpha_ERK_GTPase = 0.2; % ERK-P enhancement of RAS GTPase

% ============================================================================
% PI3K/AKT/mTOR PATHWAY PARAMETERS
% ============================================================================

% PI3K Module (slowed down by factor of 10)
k_PI3K = 0.05;         % min^-1 - PI3K activation rate by pRTK (base rate)
k_PI3K_deact = 0.005;  % min^-1 - PI3K deactivation rate
k_PIP2_to_PIP3 = 0.05; % min^-1 - PIP2 to PIP3 conversion rate
k_PTEN = 0.01;         % min^-1 - PTEN-mediated PIP3 dephosphorylation

% Cross-talk: RAS → PI3K (slowed down by factor of 10)
k_PI3K_RAS = 0.02;     % min^-1 - PI3K activation rate by RAS-GTP

% AKT Module (slowed down by factor of 10)
k_AKT_recruit = 0.03;  % min^-1 - AKT recruitment rate by PIP3
k_AKT_Thr308 = 0.05;   % min^-1 - AKT Thr308 phosphorylation rate
k_AKT_Ser473 = 0.03;   % min^-1 - AKT Ser473 phosphorylation rate
k_AKT_deact = 0.005;   % min^-1 - AKT dephosphorylation rate

% TSC1/TSC2 and Rheb Module (slowed down by factor of 10)
k_TSC2_inh = 0.05;     % min^-1 - TSC2 inhibition rate by pAKT (base rate)
k_TSC2_base = 0.01;    % min^-1 - Base TSC2 inhibition rate
k_Rheb_act = 0.05;     % min^-1 - Rheb-GTP activation rate
k_Rheb_GTPase = 0.01;  % min^-1 - Rheb GTPase activity

% Cross-talk: ERK → TSC2/mTORC1
alpha3 = 0.2;         % ERK-P enhancement of TSC2 inhibition

% mTORC1 Module (slowed down by factor of 10)
k_mTORC1 = 0.05;       % min^-1 - mTORC1 activation rate
k_mTORC1_deact = 0.005; % min^-1 - mTORC1 deactivation rate

% S6K and 4EBP1 Module (slowed down by factor of 10)
k_S6K = 0.05;          % min^-1 - S6K phosphorylation rate
k_S6K_deact = 0.005;   % min^-1 - pS6K dephosphorylation rate
k_4EBP1 = 0.05;        % min^-1 - 4EBP1 phosphorylation rate
k_4EBP1_deact = 0.005; % min^-1 - p4EBP1 dephosphorylation rate

% PI3K Feedback Parameters
alpha1 = 0.5;         % pS6K inhibition of PI3K activation

% Cross-talk: AKT → RAF (wild-type only)
alpha2 = 0.3;         % pAKT inhibition of wild-type RAF activation

% ============================================================================
% VEMURAFENIB PARAMETERS (Hill Equation)
% ============================================================================

% Hill equation parameters for vemurafenib inhibition of BRAF^V600E (normalized)
IC50_vem = 0.1;       % Normalized [0,1] - Half-maximal inhibitory concentration (IC50)
                     % At 10% of max drug concentration, 50% inhibition occurs
Hill_n = 1.5;         % Hill coefficient (cooperativity, typically 1-2 for drugs)
% Hill equation inhibition: k_BRAF_mut_eff = k_BRAF_mut * IC50^n / (IC50^n + [vemurafenib]^n)
% This gives: 0% inhibition at [Drug]=0, 50% at [Drug]=IC50, ~100% at [Drug]>>IC50
% All concentrations are normalized to [0, 1]

% Paradoxical activation parameter (vemurafenib can activate wild-type RAF) (slowed down by factor of 10)
k_RAF_paradox = 0.01;  % min^-1 - Paradoxical RAF activation rate by vemurafenib

% Store all parameters in structure
params.k_on_Shc = k_on_Shc;
params.k_off_Shc = k_off_Shc;
params.k_cat1 = k_cat1;
params.k_cat_SOS = k_cat_SOS;
params.k_on2 = k_on2;
params.k_off2 = k_off2;
params.k_on3 = k_on3;
params.k_off3 = k_off3;
params.k_deg4 = k_deg4;
params.k_SOS = k_SOS;
params.k_GTPase = k_GTPase;
params.k_RAF_act = k_RAF_act;
params.k_RAF_deact = k_RAF_deact;
params.k_BRAF_mut = k_BRAF_mut;
params.k_BRAF_mut_deact = k_BRAF_mut_deact;
params.k_MEK = k_MEK;
params.k_MEK_deact = k_MEK_deact;
params.k_ERK = k_ERK;
params.k_ERK_deact = k_ERK_deact;
params.alpha_ERK_Shc = alpha_ERK_Shc;
params.alpha_ERK_SOS = alpha_ERK_SOS;
params.alpha_ERK_RAF = alpha_ERK_RAF;
params.alpha_ERK_BRAF_mut = alpha_ERK_BRAF_mut;
params.alpha_ERK_GTPase = alpha_ERK_GTPase;
params.k_PI3K = k_PI3K;
params.k_PI3K_deact = k_PI3K_deact;
params.k_PI3K_RAS = k_PI3K_RAS;
params.k_PIP2_to_PIP3 = k_PIP2_to_PIP3;
params.k_PTEN = k_PTEN;
params.k_AKT_recruit = k_AKT_recruit;
params.k_AKT_Thr308 = k_AKT_Thr308;
params.k_AKT_Ser473 = k_AKT_Ser473;
params.k_AKT_deact = k_AKT_deact;
params.k_TSC2_inh = k_TSC2_inh;
params.k_TSC2_base = k_TSC2_base;
params.k_Rheb_act = k_Rheb_act;
params.k_Rheb_GTPase = k_Rheb_GTPase;
params.alpha3 = alpha3;
params.k_mTORC1 = k_mTORC1;
params.k_mTORC1_deact = k_mTORC1_deact;
params.k_S6K = k_S6K;
params.k_S6K_deact = k_S6K_deact;
params.k_4EBP1 = k_4EBP1;
params.k_4EBP1_deact = k_4EBP1_deact;
params.alpha1 = alpha1;
params.alpha2 = alpha2;
params.IC50_vem = IC50_vem;
params.Hill_n = Hill_n;
params.k_RAF_paradox = k_RAF_paradox;
params.pRTK_func = pRTK_func;
params.vemurafenib_func = vemurafenib_func;

fprintf('Parameters configured for both pathways with vemurafenib effects.\n\n');

%% ============================================================================
% INITIAL CONDITIONS
% ============================================================================

fprintf('Setting initial conditions (normalized to [0, 1])...\n');

% MAPK Pathway species (normalized: total pools = 1.0)
% Shc pool: Shc + pEGFR:Shc + pShc = 1.0
Shc0 = 1.0;           % Normalized: inactive Shc
pEGFR_Shc0 = 0.0;     % Normalized: pEGFR:Shc complex
pShc0 = 0.0;          % Normalized: phosphorylated Shc

% Grb2 pool: Grb2 + pShc:Grb2 = 1.0
Grb2_0 = 1.0;         % Normalized: free Grb2
pShc_Grb2_0 = 0.0;    % Normalized: pShc:Grb2 complex

% SOS pool: SOS + pShc:Grb2:SOS + SOS* = 1.0
SOS_0 = 1.0;          % Normalized: inactive SOS
pShc_Grb2_SOS_0 = 0.0; % Normalized: ternary complex
SOS_star0 = 0.0;      % Normalized: active SOS*

% RAS pool: RAS_GDP + RAS_GTP = 1.0
RAS_GDP0 = 1.0;       % Normalized: inactive RAS-GDP
RAS_GTP0 = 0.0;       % Normalized: active RAS-GTP

% RAF pool: RAF + RAF* = 1.0
RAF0 = 1.0;           % Normalized: inactive RAF
RAF_star0 = 0.0;      % Normalized: active RAF*

% BRAF^V600E pool: BRAF_mut + BRAF_mut* = 1.0
BRAF_mut0 = 0.0;      % Normalized: inactive BRAF^V600E
BRAF_mut_star0 = 1.0; % Normalized: active BRAF^V600E* (constitutively active)

% MEK pool: MEK + MEK_P = 1.0
MEK0 = 1.0;           % Normalized: inactive MEK
MEK_P0 = 0.0;         % Normalized: phosphorylated MEK

% ERK pool: ERK + ERK_P = 1.0
ERK0 = 1.0;           % Normalized: inactive ERK
ERK_P0 = 0.0;         % Normalized: phosphorylated ERK

% PI3K/AKT/mTOR Pathway species (normalized: total pools = 1.0)
% PI3K pool: PI3K + PI3K_active = 1.0
PI3K0 = 1.0;          % Normalized: inactive PI3K
PI3K_active0 = 0.0;   % Normalized: active PI3K

% PIP pool: PIP2 + PIP3 = 1.0
PIP2_0 = 1.0;         % Normalized: PIP2
PIP3_0 = 0.0;         % Normalized: PIP3

% AKT pool: AKT + pAKT_Thr308 + pAKT_Ser473 + pAKT_full = 1.0
AKT0 = 1.0;           % Normalized: inactive AKT
pAKT_Thr308_0 = 0.0;  % Normalized: AKT phosphorylated at Thr308
pAKT_Ser473_0 = 0.0;  % Normalized: AKT phosphorylated at Ser473
pAKT_full_0 = 0.0;    % Normalized: AKT phosphorylated at both sites

% TSC2 pool: TSC2_active = 1.0 (simplified, active form)
TSC2_active0 = 1.0;   % Normalized: active TSC2

% Rheb pool: Rheb_GDP + Rheb_GTP = 1.0
Rheb_GTP0 = 0.0;      % Normalized: active Rheb-GTP

% mTORC1 pool: mTORC1 + mTORC1_active = 1.0
mTORC1_0 = 1.0;       % Normalized: inactive mTORC1
mTORC1_active0 = 0.0; % Normalized: active mTORC1

% S6K pool: S6K + pS6K = 1.0
S6K0 = 1.0;           % Normalized: inactive S6K
pS6K0 = 0.0;          % Normalized: phosphorylated S6K

% 4EBP1 pool: 4EBP1 + p4EBP1 = 1.0
EBP1_0 = 1.0;         % Normalized: inactive 4EBP1
p4EBP1_0 = 0.0;       % Normalized: phosphorylated 4EBP1

% Combined initial state vector: [MAPK species (18), PI3K species (16)]
y0 = [Shc0; pEGFR_Shc0; pShc0; Grb2_0; pShc_Grb2_0; SOS_0; ...
      pShc_Grb2_SOS_0; SOS_star0; RAS_GDP0; RAS_GTP0; RAF0; ...
      RAF_star0; BRAF_mut0; BRAF_mut_star0; MEK0; MEK_P0; ERK0; ERK_P0; ...
      PI3K0; PI3K_active0; PIP2_0; PIP3_0; AKT0; pAKT_Thr308_0; ...
      pAKT_Ser473_0; pAKT_full_0; TSC2_active0; Rheb_GTP0; ...
      mTORC1_0; mTORC1_active0; S6K0; pS6K0; EBP1_0; p4EBP1_0];

fprintf('Initial conditions set. BRAF^{V600E}* initially active.\n\n');

%% ============================================================================
% TIME SPAN
% ============================================================================

tspan = [0, 200];  % minutes (0 to 50 hours)

%% ============================================================================
% SOLVE ODE SYSTEM
% ============================================================================

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   SOLVING ODE SYSTEM WITH VEMURAFENIB INHIBITION\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

fprintf('Solving ODE system (this may take a moment)...\n');
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, y] = ode45(@(t,y) vemurafenib_pathways_odes(t, y, params), tspan, y0, options);

% Extract MAPK pathway species
Shc = y(:, 1);
pEGFR_Shc = y(:, 2);
pShc = y(:, 3);
Grb2 = y(:, 4);
pShc_Grb2 = y(:, 5);
SOS = y(:, 6);
pShc_Grb2_SOS = y(:, 7);
SOS_star = y(:, 8);
RAS_GDP = y(:, 9);
RAS_GTP = y(:, 10);
RAF = y(:, 11);
RAF_star = y(:, 12);
BRAF_mut = y(:, 13);
BRAF_mut_star = y(:, 14);
MEK = y(:, 15);
MEK_P = y(:, 16);
ERK = y(:, 17);
ERK_P = y(:, 18);

% Extract PI3K/AKT/mTOR pathway species
PI3K = y(:, 19);
PI3K_active = y(:, 20);
PIP2 = y(:, 21);
PIP3 = y(:, 22);
AKT = y(:, 23);
pAKT_Thr308 = y(:, 24);
pAKT_Ser473 = y(:, 25);
pAKT_full = y(:, 26);
TSC2_active = y(:, 27);
Rheb_GTP = y(:, 28);
mTORC1 = y(:, 29);
mTORC1_active = y(:, 30);
S6K = y(:, 31);
pS6K = y(:, 32);
EBP1 = y(:, 33);
p4EBP1 = y(:, 34);

% Calculate input signals over time
pRTK_signal = arrayfun(pRTK_func, t);
vemurafenib_signal = arrayfun(vemurafenib_func, t);

fprintf('Simulation completed successfully!\n\n');

%% ============================================================================
% PLOTTING
% ============================================================================

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   GENERATING VISUALIZATIONS\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

% Figure 1: Vemurafenib Effects on MAPK Pathway
figure('Position', [50, 50, 1600, 1000]);

% MAPK Pathway - Left column
subplot(3, 3, 1);
plot(t, pShc, 'b-', 'LineWidth', 2, 'DisplayName', 'pShc');
hold on;
plot(t, SOS_star, 'g-', 'LineWidth', 2, 'DisplayName', 'SOS*');
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('MAPK: Adaptor Module', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, tspan(2)]);
hold off;

subplot(3, 3, 2);
plot(t, RAS_GTP, 'c-', 'LineWidth', 2, 'DisplayName', 'RAS-GTP');
hold on;
plot(t, RAF_star, 'm-', 'LineWidth', 2, 'DisplayName', 'RAF* (WT)');
plot(t, BRAF_mut_star, 'r-', 'LineWidth', 2, 'DisplayName', 'BRAF^{V600E}*');
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('MAPK: RAS & RAF Species', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, tspan(2)]);
hold off;

subplot(3, 3, 3);
plot(t, MEK_P, 'y-', 'LineWidth', 2, 'DisplayName', 'MEK-P');
hold on;
plot(t, ERK_P, 'r-', 'LineWidth', 2.5, 'DisplayName', 'ERK-P');
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('MAPK: MEK & ERK', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, tspan(2)]);
hold off;

% PI3K Pathway - Middle column
subplot(3, 3, 4);
plot(t, PI3K_active, 'b-', 'LineWidth', 2, 'DisplayName', 'PI3K_active');
hold on;
plot(t, PIP3, 'r-', 'LineWidth', 2, 'DisplayName', 'PIP3');
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('PI3K: PI3K & PIP3', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, tspan(2)]);
hold off;

subplot(3, 3, 5);
plot(t, pAKT_full, 'b-', 'LineWidth', 2.5, 'DisplayName', 'pAKT(Full)');
hold on;
plot(t, Rheb_GTP, 'g-', 'LineWidth', 2, 'DisplayName', 'Rheb-GTP');
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('PI3K: AKT & Rheb', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, tspan(2)]);
hold off;

subplot(3, 3, 6);
plot(t, mTORC1_active, 'm-', 'LineWidth', 2, 'DisplayName', 'mTORC1_active');
hold on;
plot(t, pS6K, 'r-', 'LineWidth', 2, 'DisplayName', 'pS6K');
plot(t, p4EBP1, 'b-', 'LineWidth', 2, 'DisplayName', 'p4EBP1');
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('PI3K: mTORC1 Outputs', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, tspan(2)]);
hold off;

% Drug and comparison plots - Right column
subplot(3, 3, 7);
plot(t, vemurafenib_signal, 'k-', 'LineWidth', 2);
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('Vemurafenib Concentration', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, tspan(2)]);
yline(0, 'k--', 'LineWidth', 1);

subplot(3, 3, 8);
plot(t, BRAF_mut_star, 'r-', 'LineWidth', 2.5, 'DisplayName', 'BRAF^{V600E}*');
hold on;
plot(t, ERK_P, 'b-', 'LineWidth', 2, 'DisplayName', 'ERK-P');
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('Vemurafenib Effect: BRAF^{V600E}* and ERK-P', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, tspan(2)]);
hold off;

subplot(3, 3, 9);
plot(t, ERK_P, 'r-', 'LineWidth', 2, 'DisplayName', 'ERK-P');
hold on;
plot(t, pS6K, 'g-', 'LineWidth', 2, 'DisplayName', 'pS6K');
plot(t, p4EBP1, 'b-', 'LineWidth', 2, 'DisplayName', 'p4EBP1');
xlabel('Time (minutes)', 'FontSize', 10);
ylabel('Normalized Concentration [0,1]', 'FontSize', 10);
title('Final Pathway Outputs', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
xlim([0, tspan(2)]);
hold off;

sgtitle('Vemurafenib Effects on MAPK and PI3K/AKT/mTOR Pathways', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Figure 2: Vemurafenib Effects and Adaptive Responses
figure('Position', [100, 100, 1400, 900]);

% Vemurafenib inhibition of BRAF^V600E
subplot(2, 3, 1);
yyaxis left;
plot(t, vemurafenib_signal, 'k-', 'LineWidth', 2);
ylabel('Normalized Vemurafenib [0,1]', 'FontSize', 12);
yyaxis right;
plot(t, BRAF_mut_star, 'r-', 'LineWidth', 2.5, 'DisplayName', 'BRAF^{V600E}*');
hold on;
plot(t, ERK_P, 'b-', 'LineWidth', 2, 'DisplayName', 'ERK-P');
ylabel('Normalized Concentration [0,1]', 'FontSize', 12);
xlabel('Time (minutes)', 'FontSize', 12);
title('Vemurafenib Inhibition of BRAF^{V600E}', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, tspan(2)]);
hold off;

% Adaptive response: RAS-GTP rebound
subplot(2, 3, 2);
plot(t, ERK_P, 'r-', 'LineWidth', 2, 'DisplayName', 'ERK-P');
hold on;
plot(t, RAS_GTP, 'c-', 'LineWidth', 2, 'DisplayName', 'RAS-GTP');
plot(t, RAF_star, 'm-', 'LineWidth', 2, 'DisplayName', 'RAF* (WT)');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Normalized Concentration [0,1]', 'FontSize', 12);
title('Adaptive Response: RAS-GTP Rebound', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, tspan(2)]);
hold off;

% PI3K pathway rebound
subplot(2, 3, 3);
plot(t, RAS_GTP, 'c-', 'LineWidth', 2, 'DisplayName', 'RAS-GTP');
hold on;
plot(t, PI3K_active, 'b-', 'LineWidth', 2, 'DisplayName', 'PI3K_active');
plot(t, PIP3, 'r-', 'LineWidth', 2, 'DisplayName', 'PIP3');
plot(t, pAKT_full, 'g-', 'LineWidth', 2, 'DisplayName', 'pAKT(Full)');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Normalized Concentration [0,1]', 'FontSize', 12);
title('PI3K Pathway Rebound (RAS → PI3K)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, tspan(2)]);
hold off;

% mTORC1 response to ERK-P reduction
subplot(2, 3, 4);
plot(t, ERK_P, 'r-', 'LineWidth', 2, 'DisplayName', 'ERK-P');
hold on;
plot(t, TSC2_active, 'g-', 'LineWidth', 2, 'DisplayName', 'TSC2_active');
plot(t, mTORC1_active, 'm-', 'LineWidth', 2, 'DisplayName', 'mTORC1_active');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Normalized Concentration [0,1]', 'FontSize', 12);
title('mTORC1 Response to ERK-P Reduction', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, tspan(2)]);
hold off;

% Effective rate constants with vemurafenib (using Hill equation)
subplot(2, 3, 5);
% Hill equation: k_eff = k_base * IC50^n / (IC50^n + [Drug]^n)
IC50_n = IC50_vem^Hill_n;
vem_n = vemurafenib_signal.^Hill_n;
% Handle zero drug concentrations
k_BRAF_mut_eff = k_BRAF_mut * IC50_n ./ (IC50_n + vem_n);
k_BRAF_mut_eff(vemurafenib_signal <= 0) = k_BRAF_mut;  % No inhibition when drug = 0
k_RAF_act_eff = k_RAF_act ./ (1 + alpha_ERK_RAF * ERK_P) ./ (1 + alpha2 * pAKT_full);
k_PI3K_eff = k_PI3K ./ (1 + alpha1 * pS6K);
k_TSC2_ERK = k_TSC2_base + alpha3 * ERK_P;
plot(t, k_BRAF_mut_eff, 'r-', 'LineWidth', 2, 'DisplayName', 'k_{BRAF^{V600E}} (inhibited)');
hold on;
plot(t, k_RAF_act_eff, 'm-', 'LineWidth', 2, 'DisplayName', 'k_{RAF} (WT)');
plot(t, k_PI3K_eff, 'b-', 'LineWidth', 2, 'DisplayName', 'k_{PI3K}');
plot(t, k_TSC2_ERK, 'g-', 'LineWidth', 2, 'DisplayName', 'k_{TSC2} (ERK-enhanced)');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Effective Rate Constant', 'FontSize', 12);
title('Rate Constants: Vemurafenib Effects', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, tspan(2)]);
hold off;

% Pathway outputs comparison
subplot(2, 3, 6);
plot(t, ERK_P, 'r-', 'LineWidth', 2.5, 'DisplayName', 'ERK-P (MAPK)');
hold on;
plot(t, pAKT_full, 'b-', 'LineWidth', 2, 'DisplayName', 'pAKT (PI3K)');
plot(t, pS6K, 'g-', 'LineWidth', 2, 'DisplayName', 'pS6K (PI3K)');
plot(t, p4EBP1, 'm-', 'LineWidth', 2, 'DisplayName', 'p4EBP1 (PI3K)');
xlabel('Time (minutes)', 'FontSize', 12);
ylabel('Normalized Concentration [0,1]', 'FontSize', 12);
title('All Pathway Outputs', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, tspan(2)]);
hold off;

sgtitle('Vemurafenib Effects: Adaptive Signaling and Cross-Talk', ...
        'FontSize', 16, 'FontWeight', 'bold');

%% ============================================================================
% SUMMARY STATISTICS
% ============================================================================

fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   SUMMARY STATISTICS\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

% Find drug addition time
drug_time_idx = find(vemurafenib_signal > 0, 1);
if isempty(drug_time_idx)
    drug_time_idx = length(t);
end
drug_time = t(drug_time_idx);

fprintf('VEMURAFENIB TREATMENT:\n');
fprintf('  Drug added at:        t = %.2f minutes\n', drug_time);
fprintf('  Peak drug concentration: %.4f (normalized [0,1])\n', max(vemurafenib_signal));

fprintf('\nMAPK PATHWAY - BEFORE vs AFTER VEMURAFENIB:\n');
% Handle case where drug starts at t=0 (drug_time_idx = 1)
if drug_time_idx > 1
    before_idx = drug_time_idx - 1;
else
    before_idx = 1;  % Use initial condition if drug starts at t=0
end
ERK_P_before = ERK_P(before_idx);
ERK_P_after = ERK_P(end);
BRAF_mut_star_before = BRAF_mut_star(before_idx);
BRAF_mut_star_after = BRAF_mut_star(end);
RAS_GTP_before = RAS_GTP(before_idx);
RAS_GTP_after = max(RAS_GTP(drug_time_idx:end));  % Peak after drug
fprintf('  BRAF^{V600E}* before:  %.4f (normalized [0,1])\n', BRAF_mut_star_before);
fprintf('  BRAF^{V600E}* after:   %.4f (normalized [0,1], %.1f%% reduction)\n', ...
        BRAF_mut_star_after, 100 * (1 - BRAF_mut_star_after / BRAF_mut_star_before));
fprintf('  ERK-P before:          %.4f (normalized [0,1])\n', ERK_P_before);
fprintf('  ERK-P after:            %.4f (normalized [0,1], %.1f%% reduction)\n', ...
        ERK_P_after, 100 * (1 - ERK_P_after / ERK_P_before));
fprintf('  RAS-GTP before:        %.4f (normalized [0,1])\n', RAS_GTP_before);
fprintf('  RAS-GTP peak after:    %.4f (normalized [0,1], %.1f%% increase)\n', ...
        RAS_GTP_after, 100 * (RAS_GTP_after / RAS_GTP_before - 1));

fprintf('\nPI3K/AKT/mTOR PATHWAY - ADAPTIVE RESPONSE:\n');
PI3K_active_before = PI3K_active(before_idx);
PI3K_active_after = max(PI3K_active(drug_time_idx:end));
pAKT_full_before = pAKT_full(before_idx);
pAKT_full_after = max(pAKT_full(drug_time_idx:end));
mTORC1_active_before = mTORC1_active(before_idx);
mTORC1_active_after = mTORC1_active(end);
fprintf('  PI3K_active before:   %.4f (normalized [0,1])\n', PI3K_active_before);
fprintf('  PI3K_active peak after: %.4f (normalized [0,1], %.1f%% change)\n', ...
        PI3K_active_after, 100 * (PI3K_active_after / PI3K_active_before - 1));
fprintf('  pAKT(Full) before:    %.4f (normalized [0,1])\n', pAKT_full_before);
fprintf('  pAKT(Full) peak after: %.4f (normalized [0,1], %.1f%% change)\n', ...
        pAKT_full_after, 100 * (pAKT_full_after / pAKT_full_before - 1));
fprintf('  mTORC1_active before: %.4f (normalized [0,1])\n', mTORC1_active_before);
fprintf('  mTORC1_active after:  %.4f (normalized [0,1], %.1f%% change)\n', ...
        mTORC1_active_after, 100 * (mTORC1_active_after / mTORC1_active_before - 1));

fprintf('\nCROSS-TALK AND FEEDBACK STRENGTHS:\n');
fprintf('  Vemurafenib inhibition (Hill): IC50 = %.4f (normalized [0,1]), Hill n = %.2f\n', IC50_vem, Hill_n);
fprintf('  RAS → PI3K:           k_PI3K_RAS = %.2f\n', k_PI3K_RAS);
fprintf('  ERK → TSC2:           alpha3 = %.2f\n', alpha3);
fprintf('  AKT → RAF (WT only):  alpha2 = %.2f\n', alpha2);
fprintf('  pS6K → PI3K:         alpha1 = %.2f\n', alpha1);

fprintf('\n═══════════════════════════════════════════════════════════════════════════\n');
fprintf('   SIMULATION COMPLETED SUCCESSFULLY!\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n\n');

%% ============================================================================
% SAVE OUTPUTS
% ============================================================================

save('vemurafenib_output.mat', 't', 'ERK_P', 'pAKT_full', 'pS6K', 'p4EBP1', ...
     'BRAF_mut_star', 'RAF_star', 'MEK_P', 'RAS_GTP', 'PIP3', 'Rheb_GTP', ...
     'mTORC1_active', 'pRTK_signal', 'vemurafenib_signal');
fprintf('Saved outputs to vemurafenib_output.mat\n');

%% ============================================================================
% ODE FUNCTION
% ============================================================================

function dydt = vemurafenib_pathways_odes(t, y, p)
    % Integrated ODE system for MAPK and PI3K/AKT/mTOR pathways with vemurafenib
    % y = [MAPK: Shc, pEGFR_Shc, pShc, Grb2, pShc_Grb2, SOS, pShc_Grb2_SOS, SOS_star,
    %      RAS_GDP, RAS_GTP, RAF, RAF_star, BRAF_mut, BRAF_mut_star, MEK, MEK_P, ERK, ERK_P;
    %      PI3K: PI3K, PI3K_active, PIP2, PIP3, AKT, pAKT_Thr308, pAKT_Ser473, pAKT_full,
    %      TSC2_active, Rheb_GTP, mTORC1, mTORC1_active, S6K, pS6K, 4EBP1, p4EBP1]
    
    % Extract MAPK pathway species
    Shc = y(1);
    pEGFR_Shc = y(2);
    pShc = y(3);
    Grb2 = y(4);
    pShc_Grb2 = y(5);
    SOS = y(6);
    pShc_Grb2_SOS = y(7);
    SOS_star = y(8);
    RAS_GDP = y(9);
    RAS_GTP = y(10);
    RAF = y(11);
    RAF_star = y(12);
    BRAF_mut = y(13);
    BRAF_mut_star = y(14);
    MEK = y(15);
    MEK_P = y(16);
    ERK = y(17);
    ERK_P = y(18);
    
    % Extract PI3K/AKT/mTOR pathway species
    PI3K = y(19);
    PI3K_active = y(20);
    PIP2 = y(21);
    PIP3 = y(22);
    AKT = y(23);
    pAKT_Thr308 = y(24);
    pAKT_Ser473 = y(25);
    pAKT_full = y(26);
    TSC2_active = y(27);
    Rheb_GTP = y(28);
    mTORC1 = y(29);
    mTORC1_active = y(30);
    S6K = y(31);
    pS6K = y(32);
    EBP1 = y(33);
    p4EBP1 = y(34);
    
    % Get input signals at current time
    pRTK = p.pRTK_func(t);
    vemurafenib = p.vemurafenib_func(t);
    
    % ========================================================================
    % CALCULATE FEEDBACK-MODULATED AND CROSS-TALK RATE CONSTANTS
    % ========================================================================
    
    % MAPK pathway feedbacks (ERK-P) - feedback is relieved when ERK-P drops
    k_on_Shc_eff = p.k_on_Shc / (1 + p.alpha_ERK_Shc * ERK_P);
    k_cat_SOS_eff = p.k_cat_SOS / (1 + p.alpha_ERK_SOS * ERK_P);
    k_RAF_act_eff = p.k_RAF_act / (1 + p.alpha_ERK_RAF * ERK_P);
    k_GTPase_eff = p.k_GTPase + p.alpha_ERK_GTPase * ERK_P;
    
    % Cross-talk: AKT → RAF (wild-type only, not BRAF^V600E)
    k_RAF_act_eff = k_RAF_act_eff / (1 + p.alpha2 * pAKT_full);
    
    % VEMURAFENIB INHIBITION: BRAF^V600E activity reduced using Hill equation
    % Hill equation: k_eff = k_base * IC50^n / (IC50^n + [Drug]^n)
    % This gives sigmoidal inhibition curve with IC50 as half-maximal concentration
    % When [Drug] = 0: k_eff = k_base (no inhibition)
    % When [Drug] = IC50: k_eff = k_base/2 (50% inhibition)
    % When [Drug] >> IC50: k_eff → 0 (full inhibition)
    if vemurafenib <= 0
        k_BRAF_mut_eff = p.k_BRAF_mut;  % No inhibition when drug concentration is zero
    else
        IC50_n = p.IC50_vem^p.Hill_n;
        vem_n = vemurafenib^p.Hill_n;
        k_BRAF_mut_eff = p.k_BRAF_mut * IC50_n / (IC50_n + vem_n);
    end
    
    % BRAF^V600E feedback (much weaker than wild-type)
    k_BRAF_mut_eff = k_BRAF_mut_eff / (1 + p.alpha_ERK_BRAF_mut * ERK_P);
    
    % Paradoxical activation: Vemurafenib can trans-activate wild-type RAF
    % (simplified: small increase in RAF activation when vemurafenib is present)
    if vemurafenib > 0
        k_RAF_paradox_eff = p.k_RAF_paradox * vemurafenib / (1 + vemurafenib);
    else
        k_RAF_paradox_eff = 0;
    end
    
    % PI3K pathway feedback (pS6K)
    k_PI3K_eff = p.k_PI3K / (1 + p.alpha1 * pS6K);
    
    % Cross-talk: ERK → TSC2 (ERK-P enhances TSC2 inhibition)
    % This is reduced when ERK-P drops after vemurafenib
    k_TSC2_ERK = p.k_TSC2_base + p.alpha3 * ERK_P;
    
    % ========================================================================
    % MAPK PATHWAY REACTIONS
    % ========================================================================
    
    % Shc–Grb2–SOS module
    r1_forward = k_on_Shc_eff * pRTK * Shc;
    r1_reverse = p.k_off_Shc * pEGFR_Shc;
    r1_cat = p.k_cat1 * pEGFR_Shc;
    r2_forward = p.k_on2 * pShc * Grb2;
    r2_reverse = p.k_off2 * pShc_Grb2;
    r3_forward = p.k_on3 * pShc_Grb2 * SOS;
    r3_reverse = p.k_off3 * pShc_Grb2_SOS;
    r3_cat = k_cat_SOS_eff * pShc_Grb2_SOS;
    r4 = p.k_deg4 * SOS_star;
    
    % RAS–RAF–MEK–ERK module
    r5 = p.k_SOS * SOS_star * RAS_GDP;
    r6 = k_GTPase_eff * RAS_GTP;
    
    % Wild-type RAF activation (RAS-dependent, AKT-inhibited, with paradoxical activation)
    r7 = k_RAF_act_eff * RAS_GTP * RAF + k_RAF_paradox_eff * RAF;  % Paradoxical activation
    r8 = p.k_RAF_deact * RAF_star;
    
    % BRAF^V600E constitutive activation (inhibited by vemurafenib)
    r9 = k_BRAF_mut_eff * BRAF_mut;
    r10 = p.k_BRAF_mut_deact * BRAF_mut_star;
    
    % MEK phosphorylation by both RAF* and BRAF^V600E*
    r11 = p.k_MEK * RAF_star * MEK;  % Wild-type RAF contribution
    r12 = p.k_MEK * BRAF_mut_star * MEK;  % BRAF^V600E contribution (reduced by vemurafenib)
    r13 = p.k_MEK_deact * MEK_P;
    
    % ERK phosphorylation
    r14 = p.k_ERK * MEK_P * ERK;
    r15 = p.k_ERK_deact * ERK_P;
    
    % ========================================================================
    % PI3K/AKT/mTOR PATHWAY REACTIONS
    % ========================================================================
    
    % PI3K module
    % PI3K activation by pRTK (with pS6K feedback) AND by RAS-GTP (cross-talk)
    % RAS-GTP rebound after vemurafenib can increase PI3K activation
    r16 = k_PI3K_eff * pRTK * PI3K + p.k_PI3K_RAS * RAS_GTP * PI3K;  % Cross-talk: RAS → PI3K
    r17 = p.k_PI3K_deact * PI3K_active;
    r18 = p.k_PIP2_to_PIP3 * PI3K_active * PIP2;
    r19 = p.k_PTEN * PIP3;
    
    % AKT module
    r20 = p.k_AKT_recruit * PIP3 * AKT;
    r21 = p.k_AKT_Thr308 * PIP3 * AKT;
    r22 = p.k_AKT_Ser473 * pAKT_Thr308;
    r23 = p.k_AKT_Ser473 * PIP3 * AKT;
    r24 = p.k_AKT_Thr308 * pAKT_Ser473;
    r25_Thr308 = p.k_AKT_deact * pAKT_Thr308;
    r25_Ser473 = p.k_AKT_deact * pAKT_Ser473;
    r25_full = p.k_AKT_deact * pAKT_full;
    
    % TSC2/Rheb module
    % TSC2 inhibition by pAKT (base) AND by ERK-P (cross-talk, reduced when ERK-P drops)
    r26 = (p.k_TSC2_inh * pAKT_full + k_TSC2_ERK * ERK_P) * TSC2_active;  % Cross-talk: ERK → TSC2
    Rheb_GDP = 1.0 - Rheb_GTP;  % Normalized: total Rheb pool = 1.0
    r27 = p.k_Rheb_act * Rheb_GDP / (1 + TSC2_active);
    r28 = p.k_Rheb_GTPase * Rheb_GTP * (1 + TSC2_active);  % Updated for normalized values
    
    % mTORC1 module
    r29 = p.k_mTORC1 * Rheb_GTP * mTORC1;
    r30 = p.k_mTORC1_deact * mTORC1_active;
    
    % S6K and 4EBP1 module
    r31 = p.k_S6K * mTORC1_active * S6K;
    r32 = p.k_S6K_deact * pS6K;
    r33 = p.k_4EBP1 * mTORC1_active * EBP1;
    r34 = p.k_4EBP1_deact * p4EBP1;
    
    % ========================================================================
    % ODEs FOR MAPK PATHWAY
    % ========================================================================
    
    dShc_dt = -r1_forward + r1_reverse;
    dpEGFR_Shc_dt = r1_forward - r1_reverse - r1_cat;
    dpShc_dt = r1_cat - r2_forward + r2_reverse;
    dGrb2_dt = -r2_forward + r2_reverse;
    dpShc_Grb2_dt = r2_forward - r2_reverse - r3_forward + r3_reverse + r3_cat;
    dSOS_dt = -r3_forward + r3_reverse + r4;
    dpShc_Grb2_SOS_dt = r3_forward - r3_reverse - r3_cat;
    dSOS_star_dt = r3_cat - r4;
    dRAS_GDP_dt = -r5 + r6;
    dRAS_GTP_dt = r5 - r6 - r7;
    dRAF_dt = -r7 + r8;
    dRAF_star_dt = r7 - r8;
    dBRAF_mut_dt = -r9 + r10;
    dBRAF_mut_star_dt = r9 - r10;
    dMEK_dt = -r11 - r12 + r13;
    dMEK_P_dt = r11 + r12 - r13;
    dERK_dt = -r14 + r15;
    dERK_P_dt = r14 - r15;
    
    % ========================================================================
    % ODEs FOR PI3K/AKT/mTOR PATHWAY
    % ========================================================================
    
    dPI3K_dt = -r16 + r17;
    dPI3K_active_dt = r16 - r17;
    dPIP2_dt = -r18 + r19;
    dPIP3_dt = r18 - r19;
    dAKT_dt = -r20 - r21 - r23 + r25_Thr308 + r25_Ser473 + r25_full;
    dpAKT_Thr308_dt = r20 + r21 - r22 - r25_Thr308;
    dpAKT_Ser473_dt = r23 - r24 - r25_Ser473;
    dpAKT_full_dt = r22 + r24 - r25_full;
    dTSC2_active_dt = -r26;
    dRheb_GTP_dt = r27 - r28;
    dmTORC1_dt = -r29 + r30;
    dmTORC1_active_dt = r29 - r30;
    dS6K_dt = -r31 + r32;
    dpS6K_dt = r31 - r32;
    dEBP1_dt = -r33 + r34;
    dp4EBP1_dt = r33 - r34;
    
    % Return all derivatives
    dydt = [dShc_dt; dpEGFR_Shc_dt; dpShc_dt; dGrb2_dt; dpShc_Grb2_dt; ...
            dSOS_dt; dpShc_Grb2_SOS_dt; dSOS_star_dt; ...
            dRAS_GDP_dt; dRAS_GTP_dt; dRAF_dt; dRAF_star_dt; ...
            dBRAF_mut_dt; dBRAF_mut_star_dt; dMEK_dt; dMEK_P_dt; ...
            dERK_dt; dERK_P_dt; ...
            dPI3K_dt; dPI3K_active_dt; dPIP2_dt; dPIP3_dt; ...
            dAKT_dt; dpAKT_Thr308_dt; dpAKT_Ser473_dt; dpAKT_full_dt; ...
            dTSC2_active_dt; dRheb_GTP_dt; dmTORC1_dt; dmTORC1_active_dt; ...
            dS6K_dt; dpS6K_dt; dEBP1_dt; dp4EBP1_dt];
end

