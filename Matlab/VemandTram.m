ka1 = 10e-5;
kr1 = 10e-6;
kc1 = 10e-3;

kpCraf = 10e-4;
kpMek = 7e-6;
kpErk = 20e-5;

kDegradEgfr = 10e-7;
kErkInbEgfr = 20e-6;
kShcDephos = 10e-5;
kptpDeg = 7e-7;
kGrb2CombShc = 36e-6;
kSprtyInbGrb2 = 10e-5;
kSosCombGrb2 = 10e-5;
kErkPhosSos = 17e-6;
kErkPhosPcraf = 15.7342e-6;
kPcrafDegrad = 1.7342e-4;
kErkPhosMek =    20.7342e-5;
kMekDegrad = 1.7342e-4;
kDuspInbErk =  2.7342e-4;
kErkDeg = 1.7342e-7;
kinbBraf = 10e-2;
kDuspStop =    1.7342e-04;
kDusps =   1.11e-4 ;
kSproutyForm = 1.7e-06;
kSprtyComeDown =    5.5000e-5;
km_Sprty_decay =8e-5;
kdegrad = 10e-7;
km_Dusp = 1.11e-4;
km_Sprty = 1.7e-06;

% Enhanced DUSP-ERK interaction parameters
kErkDephos = 1e-5;  % ERK dephosphorylation by DUSP
kDuspDeg = 1e-6;    % DUSP degradation by ERK

kAkt = 10e-5;
kdegradAKT = 10e-7;
kb1 = 10e-8;
k43b1 = 10e-3;

% --- KSR + MEK phosphorylation route params
kKSRphos      = 5e-6;   % RAF^P-driven phosphorylation of KSR
kKSRdephos    = 5e-6;   % dephosphorylation of pKSR
kMekByBraf    = 8e-6;   % BRAF^P -> MEK phosphorylation (explicit, instead of reusing ka1)
kMekByCraf    = 8e-6;
kMekByKSR     = 3e-6;   % pKSR   -> MEK phosphorylation (weaker than RAF route)

% --- Trametinib (MEK inhibitor) Hill modifiers
Tram          = 1e-6;   % trametinib concentration (use any schedule you like)
K_tram_RAF    = 5e-1;   % Ki-like for inhibiting RAF→MEK (less potent)
K_tram_KSR    = 5e-12;  % Ki-like for inhibiting KSR→MEK (more potent)
n_tram        = 2;      % Hill coefficient

sotrasinb = 0;
kTramInb = 0;
%Pi3k Pathway
kEGFRPI3k = 10e-4;

kKrasActivatePI3k = 10e-5;
kMTOR_Feedback = 10e-4;

kIGFR = 5e-5;


tspan = [0 2000000];

EGFR = [1,0,0];
SHC = [1,0,1];
grb_sos = [0,0];
HRAS = [0,0];
NRAS = [0,0];
KRAS = [1,0,1];
CRAF = [0.8,0.2];
BRAF = [0,1];
MEK = [1,1];
ERK = [1,0.8];
DUSP = [1,1];
SRTY = [1,1];
pERKDegrad = [1];
pMEKDegrad = [1];
pCRAFDegrad = [1];
DUSPStop = [1];

%PI3k
IGFR =[1,0,0];
IRS = [1,0];
pi3k = [1,0];
AKT = [1,0];
FOXO = [0];
MTORC = [1,0];
frebp1 = [1,0,0];
KSR   = [1, 0];  % KSR, pKSR


y0 = [EGFR,SHC,grb_sos,HRAS,NRAS,KRAS,CRAF,BRAF,MEK,ERK,DUSP,SRTY,pERKDegrad,pMEKDegrad,pCRAFDegrad,DUSPStop,IGFR,IRS,pi3k,AKT,FOXO,MTORC,frebp1,KSR]; 


% All experimental timestamps (0, 1, 4, 8, 24, 48 hours)
timeStamps_all = [0,1,4,8,24,48];

% Combined experimental data for all timestamps
RASGtpValsExp_all = [0.558156831,0.61832384,0.614585492,0.708019641,0.866573533,1];
expMEK_all = [1,0.015974647,0.006869778,0.0199368,0.056852335,0.537465586];
PexpErkVals_all = [1,0.024357316,0.003652955,0.002236656,0.01479824,0.308323021];
expDusp_all = [0.953521721,1,0.401783008,0.126565633,0.011764517,0.034619089];
expPEGFR_all = [0.339054463,1,0.991433301,0.92624324,0.039717869,0.016855968];
expCRAF_all = [0.346766996,1,0.86761608,0.684515481,0.997892098,0.409290248];
expAKT_all = [0.514573119,0.791478835,0.944188475,1,0.506364117,0.215157296];




% Initial guess for the parameters
% ----- build parameter vector (now length 46) -----
params0 = [ka1,kr1,kc1, ...
    kpCraf,kpMek,kpErk, ...
    kDegradEgfr,kErkInbEgfr,kShcDephos,kptpDeg,kGrb2CombShc, ...
    kSprtyInbGrb2,kSosCombGrb2,kErkPhosSos, ...
    kErkPhosPcraf,kPcrafDegrad,kErkPhosMek,kMekDegrad, ...
    kDuspInbErk,kErkDeg,kinbBraf,kDuspStop,kDusps, ...
    kSproutyForm,kSprtyComeDown,kdegrad,km_Sprty_decay,km_Dusp,km_Sprty, ...
    kErkDephos,kDuspDeg, ...
    kEGFRPI3k,kMTOR_Feedback, ...
    kAkt,kdegradAKT, ...
    kb1,k43b1, ...
    kKSRphos,kKSRdephos, ...
    kMekByBraf,kMekByCraf,kMekByKSR, ...
    Tram,K_tram_RAF,K_tram_KSR,n_tram];

% bounds must match params0 length (46 parameters)
% Create protein-specific bounds for better optimization
lb = 1e-8 * ones(size(params0));
ub = 1e-4 * ones(size(params0));

% Enhanced bounds for DUSP-ERK interaction parameters
% kErkDephos (parameter 30) - ERK dephosphorylation by DUSP
lb(30) = 1e-6;  % Higher lower bound for meaningful dephosphorylation
ub(30) = 1e-3;  % Higher upper bound for strong dephosphorylation

% kDuspDeg (parameter 31) - DUSP degradation by ERK  
lb(31) = 1e-7;  % Lower bound for DUSP degradation
ub(31) = 1e-4;  % Upper bound for DUSP degradation

% options
opts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','MaxIterations',150);

% fmincon call — optimize on all timestamps (0-48hrs)
[optimizedParams, errorOpt] = fmincon( ...
    @(p) objectiveFunction_all(p, timeStamps_all, RASGtpValsExp_all, expMEK_all, PexpErkVals_all, expDusp_all, expPEGFR_all, expCRAF_all, expAKT_all, y0), ...
    params0, [], [], [], [], lb, ub, [], opts);

disp('Optimized Parameters:');  disp(optimizedParams);
disp('Minimum Fit Error (All Timestamps):');         disp(errorOpt);

% Simulate with optimized parameters for all timepoints
[T_all, Y_all] = ode45(@(t,y) Mapk( ...
    t,y, ...
    optimizedParams(1),optimizedParams(2),optimizedParams(3), ...
    optimizedParams(4),optimizedParams(5),optimizedParams(6), ...
    optimizedParams(7),optimizedParams(8),optimizedParams(9),optimizedParams(10),optimizedParams(11), ...
    optimizedParams(12),optimizedParams(13),optimizedParams(14), ...
    optimizedParams(15),optimizedParams(16),optimizedParams(17),optimizedParams(18), ...
    optimizedParams(19),optimizedParams(20),optimizedParams(21),optimizedParams(22),optimizedParams(23), ...
    optimizedParams(24),optimizedParams(25),optimizedParams(26),optimizedParams(27),optimizedParams(28),optimizedParams(29), ...
    optimizedParams(30),optimizedParams(31), ...
    optimizedParams(32),optimizedParams(33), ...
    optimizedParams(34),optimizedParams(35), ...
    optimizedParams(36),optimizedParams(37), ...
    optimizedParams(38),optimizedParams(39),optimizedParams(40), ...
    optimizedParams(41),optimizedParams(42),optimizedParams(43),optimizedParams(44),optimizedParams(45),optimizedParams(46)), ...
    timeStamps_all*3600, y0);

% Normalization helper
normit = @(v,mode) v ./ max(eps, ...
    (strcmpi(mode,'first'))*v(1) + ...
    (strcmpi(mode,'last')) *v(end) + ...
    (strcmpi(mode,'max'))  *max(v));

% Extract and normalize model outputs at all timepoints
m_pEGFR_all = normit(Y_all(:,3),'max');
m_RASgtp_all = normit(Y_all(:,14),'last');
m_pCRAF_all = normit(Y_all(:,17),'max');
m_pMEK_all = normit(Y_all(:,21),'first');
m_pERK_all = normit(Y_all(:,23),'first');
m_DUSP_all = normit(Y_all(:,24),'max');
m_pAKT_all = normit(Y_all(:,40),'max');

% Calculate fit errors at all timepoints
fit_errors = struct();
fit_errors.pEGFR = abs(m_pEGFR_all - expPEGFR_all(:));
fit_errors.RAS = abs(m_RASgtp_all - RASGtpValsExp_all(:));
fit_errors.pCRAF = abs(m_pCRAF_all - expCRAF_all(:));
fit_errors.pMEK = abs(m_pMEK_all - expMEK_all(:));
fit_errors.pERK = abs(m_pERK_all - PexpErkVals_all(:));
fit_errors.DUSP = abs(m_DUSP_all - expDusp_all(:));
fit_errors.pAKT = abs(m_pAKT_all - expAKT_all(:));

% Display fit results
disp('=== MODEL FIT RESULTS (All Timestamps) ===');
disp(['Total Fit Error: ', num2str(errorOpt, '%.4f')]);
disp('Individual Protein Errors (sum across all timepoints):');
disp(['pEGFR: ', num2str(sum(fit_errors.pEGFR), '%.4f')]);
disp(['RAS: ', num2str(sum(fit_errors.RAS), '%.4f')]);
disp(['pCRAF: ', num2str(sum(fit_errors.pCRAF), '%.4f')]);
disp(['pMEK: ', num2str(sum(fit_errors.pMEK), '%.4f')]);
disp(['pERK: ', num2str(sum(fit_errors.pERK), '%.4f')]);
disp(['DUSP: ', num2str(sum(fit_errors.DUSP), '%.4f')]);
disp(['pAKT: ', num2str(sum(fit_errors.pAKT), '%.4f')]);


function dy = Mapk(t,y, ...
    ka1,kr1,kc1,kpCraf,kpMek,kpErk, ...
    kDegradEgfr,kErkInbEgfr,kShcDephos,kptpDeg,kGrb2CombShc, ...
    kSprtyInbGrb2,kSosCombGrb2,kErkPhosSos, ...
    kErkPhosPcraf,kPcrafDegrad,kErkPhosMek,kMekDegrad, ...
    kDuspInbErk,kErkDeg,kinbBraf,kDuspStop,kDusps, ...
    kSproutyForm,kSprtyComeDown,kdegrad,km_Sprty_decay,km_Dusp,km_Sprty, ...
    kErkDephos,kDuspDeg, ...
    kEGFRPI3k,kMTOR_Feedback, ...
    kAkt,kdegradAKT, ...
    kb1,k43b1, ...
    kKSRphos,kKSRdephos, ...
    kMekByBraf,kMekByCraf,kMekByKSR, ...
    Tram,K_tram_RAF,K_tram_KSR,n_tram)

dy = zeros(48,1);

% Hill-scale trametinib factors (0..1)
f_tram_RAF = 1 / (1 + (Tram / K_tram_RAF)^n_tram);
f_tram_KSR = 1 / (1 + (Tram / K_tram_KSR)^n_tram);

% --- RTK to SOS ---
dy(1) =  -ka1*y(1) + kr1*y(2);
dy(2) =   ka1*y(1) - kr1*y(2) - kc1*y(2);
dy(3) =   kc1*y(2) - kDegradEgfr*y(3) - kErkInbEgfr*y(23)*y(3);

dy(4) =  -ka1*y(3)*y(4);
dy(5) =   ka1*y(3)*y(4) - kShcDephos*y(6)*y(5);
dy(6) =  -kptpDeg*y(5)*y(6);

dy(7) =   kGrb2CombShc*y(5)*y(3) - kSprtyInbGrb2*y(21)*y(7);
dy(8) =   kSosCombGrb2*y(7)*y(5) - kErkPhosSos*y(19)*y(8)  ;

% --- RAS blocks (use KRAS as main activator) ---
dy(9)  = -ka1*y(8)*y(9);
dy(10) =  ka1*y(8)*y(9);
dy(11) = -ka1*y(8)*y(11);
dy(12) =  ka1*y(8)*y(11);
dy(13) = -ka1*y(8)*y(13);
dy(14) =  ka1*y(8)*y(13) - ka1*y(14)*y(15);
dy(15) = -ka1*y(14)*y(15);

% --- RAFs ---
dy(16) = -kpCraf*y(14)*y(16) + kErkPhosPcraf*y(23)*y(17) + kPcrafDegrad*y(17)*y(30);
dy(17) =  kpCraf*y(14)*y(16) - kErkPhosPcraf*y(23)*y(17) - kPcrafDegrad*y(17)*y(30);
dy(18) = -ka1*y(18);
dy(19) =  ka1*y(18) - kinbBraf*y(19);

% --- MEK (RAF routes + KSR route + feedbacks/degrad) ---
raf_to_mek = ( kpMek*y(17) + kMekByBraf*y(19) + kMekByCraf*y(17) ) * f_tram_RAF;
ksr_to_mek = ( kMekByKSR * y(48) ) * f_tram_KSR;

dy(20) = -(raf_to_mek + ksr_to_mek)*y(20) + kErkPhosMek*y(23)*y(21) + kMekDegrad*y(21)*y(29);
dy(21) =  (raf_to_mek + ksr_to_mek)*y(20) - kErkPhosMek*y(23)*y(21) - kMekDegrad*y(21)*y(29);

% --- ERK ---
% Enhanced ERK dynamics with proper dephosphorylation by DUSP
dy(22) = -kpErk*y(21)*y(22) + kDuspInbErk*y(24)*y(23) + kErkDeg*y(23)*y(28) + kErkDephos*y(24)*y(23);
dy(23) =  kpErk*y(21)*y(22) - kDuspInbErk*y(24)*y(23) - kErkDeg*y(23)*y(28) - kErkDephos*y(24)*y(23);

% --- DUSP & Sprouty ---
% Enhanced DUSP dynamics with ERK-dependent degradation
dy(24) =  km_Dusp*y(23)/(1 + (km_Dusp/kDusps)*y(23)) - kDuspStop*y(24)*y(31) - kDuspDeg*y(24)*y(23);
dy(25) = -kDuspStop*y(24)*y(25);
dy(26) =  km_Sprty*y(23)/(1 + (km_Sprty/kSproutyForm)*y(23)) - kSprtyComeDown*y(26)*y(27);
dy(27) = -kSprtyComeDown*y(26)*y(27);

% --- degrad trackers ---
dy(28) = -kErkDeg*y(23)*y(28);
dy(29) = -kMekDegrad*y(21)*y(29);
dy(30) = -kPcrafDegrad*y(17)*y(30);
dy(31) = -kDuspStop*y(24)*y(31);

% --- IGF/PI3K/AKT/MTOR ---
dy(32) = -ka1*y(32) + kr1*y(33);
dy(33) =  ka1*y(32) - kr1*y(33) - kc1*y(33);
dy(34) =  kc1*y(33) - kErkInbEgfr*y(23)*y(34);

dy(35) = -ka1*y(3)*y(35) - kEGFRPI3k*y(34)*y(35);
dy(36) =  ka1*y(3)*y(35) + kEGFRPI3k*y(34)*y(35) - kMTOR_Feedback*y(43)*y(36);

dy(37) = -(ka1*y(36)*y(37) + ka1*y(14)*y(37));
dy(38) =  ka1*y(36)*y(37) + ka1*y(14)*y(37) - kdegrad*y(38);

dy(39) = -kAkt*y(38)*y(39);
dy(40) =  kAkt*y(38)*y(39) - kdegradAKT*y(40);

dy(41) =  (1-y(40))*ka1/(1 + (y(41)/15e-5));
dy(42) = -ka1*y(40)*y(42) + kdegrad*y(43);
dy(43) =  ka1*y(40)*y(42) - kdegrad*y(43);

dy(44) = -ka1*y(43)*y(44) + kb1*y(44);
dy(45) =  ka1*y(43)*y(44) - kb1*y(44) - k43b1*y(45);
dy(46) =  k43b1*y(45);

% --- KSR ---
dy(47) = -kKSRphos*(y(17)+y(19))*y(47) + kKSRdephos*y(48);
dy(48) =  kKSRphos*(y(17)+y(19))*y(47) - kKSRdephos*y(48);
end

% Objective function for fitting on all timestamps
function err = objectiveFunction_all(p, timeStamps_all, RASGtpValsExp_all, expMEK_all, PexpErkVals_all, expDusp_all, expPEGFR_all, expCRAF_all, expAKT_all, y0)

% Unpack params in the SAME order as params0
ka1=p(1);  kr1=p(2);  kc1=p(3);
kpCraf=p(4); kpMek=p(5); kpErk=p(6);
kDegradEgfr=p(7); kErkInbEgfr=p(8); kShcDephos=p(9); kptpDeg=p(10); kGrb2CombShc=p(11);
kSprtyInbGrb2=p(12); kSosCombGrb2=p(13); kErkPhosSos=p(14);
kErkPhosPcraf=p(15); kPcrafDegrad=p(16); kErkPhosMek=p(17); kMekDegrad=p(18);
kDuspInbErk=p(19); kErkDeg=p(20); kinbBraf=p(21); kDuspStop=p(22); kDusps=p(23);
kSproutyForm=p(24); kSprtyComeDown=p(25); kdegrad=p(26); km_Sprty_decay=p(27); km_Dusp=p(28); km_Sprty=p(29);
kErkDephos=p(30); kDuspDeg=p(31);
kEGFRPI3k=p(32); kMTOR_Feedback=p(33);
kAkt=p(34); kdegradAKT=p(35);
kb1=p(36); k43b1=p(37);
kKSRphos=p(38); kKSRdephos=p(39);
kMekByBraf=p(40); kMekByCraf=p(41); kMekByKSR=p(42);
Tram=p(43); K_tram_RAF=p(44); K_tram_KSR=p(45); n_tram=p(46);

% integrate exactly at all experimental times (seconds)
[T, Y] = ode45(@(t,y) Mapk( ...
    t,y, ...
    ka1,kr1,kc1,kpCraf,kpMek,kpErk, ...
    kDegradEgfr,kErkInbEgfr,kShcDephos,kptpDeg,kGrb2CombShc, ...
    kSprtyInbGrb2,kSosCombGrb2,kErkPhosSos, ...
    kErkPhosPcraf,kPcrafDegrad,kErkPhosMek,kMekDegrad, ...
    kDuspInbErk,kErkDeg,kinbBraf,kDuspStop,kDusps, ...
    kSproutyForm,kSprtyComeDown,kdegrad,km_Sprty_decay,km_Dusp,km_Sprty, ...
    kErkDephos,kDuspDeg, ...
    kEGFRPI3k,kMTOR_Feedback, ...
    kAkt,kdegradAKT, ...
    kb1,k43b1, ...
    kKSRphos,kKSRdephos, ...
    kMekByBraf,kMekByCraf,kMekByKSR, ...
    Tram,K_tram_RAF,K_tram_KSR,n_tram), ...
    timeStamps_all*3600, y0);

% pull model outputs at all experimental times
m_pEGFR = Y(:,3);
m_RASgtp = Y(:,14);
m_pCRAF = Y(:,17);
m_pMEK  = Y(:,21);
m_pERK  = Y(:,23);
m_ERK   = Y(:,22);
m_DUSP  = Y(:,24);
m_pAKT  = Y(:,40);

% Normalization (matches original)
norm = @(v,mode) v ./ max(eps, (mode=="first")*v(1) + (mode=="last")*v(end) + (mode=="max")*max(v));
m_pEGFR = norm(m_pEGFR,'max');
m_RASgtp= norm(m_RASgtp,'last');
m_pMEK  = norm(m_pMEK,'first');
m_pERK  = norm(m_pERK,'first');
m_DUSP  = norm(m_DUSP,'max');
m_pCRAF = norm(m_pCRAF,'max');
m_pAKT  = norm(m_pAKT,'max');
m_ERK   = norm(m_ERK,'max');

% squared-error across all experimental data
% Enhanced weights for DUSP-ERK interaction optimization
w = struct('EGFR',1,'RAS',1,'MEK',1,'ERK',3,'DUSP',3,'CRAF',1,'AKT',1);

err = 0;
err = err + w.EGFR * sum((m_pEGFR - expPEGFR_all(:)).^2);
err = err + w.RAS  * sum((m_RASgtp - RASGtpValsExp_all(:)).^2);
err = err + w.MEK  * sum((m_pMEK  - expMEK_all(:)).^2);
err = err + w.ERK  * sum((m_pERK  - PexpErkVals_all(:)).^2);
err = err + w.DUSP * sum((m_DUSP  - expDusp_all(:)).^2);
err = err + w.CRAF * sum((m_pCRAF - expCRAF_all(:)).^2);
err = err + w.AKT  * sum((m_pAKT  - expAKT_all(:)).^2);
end

function err = objectiveFunction(p, timeStamps, RASGtpValsExp, expMEK, PexpErkVals, expDusp, expPEGFR, expCRAF, expAKT, y0)

% Unpack params in the SAME order as params0
ka1=p(1);  kr1=p(2);  kc1=p(3);
kpCraf=p(4); kpMek=p(5); kpErk=p(6);
kDegradEgfr=p(7); kErkInbEgfr=p(8); kShcDephos=p(9); kptpDeg=p(10); kGrb2CombShc=p(11);
kSprtyInbGrb2=p(12); kSosCombGrb2=p(13); kErkPhosSos=p(14);
kErkPhosPcraf=p(15); kPcrafDegrad=p(16); kErkPhosMek=p(17); kMekDegrad=p(18);
kDuspInbErk=p(19); kErkDeg=p(20); kinbBraf=p(21); kDuspStop=p(22); kDusps=p(23);
kSproutyForm=p(24); kSprtyComeDown=p(25); kdegrad=p(26); km_Sprty_decay=p(27); km_Dusp=p(28); km_Sprty=p(29);
kErkDephos=p(30); kDuspDeg=p(31);
kEGFRPI3k=p(32); kMTOR_Feedback=p(33);
kAkt=p(34); kdegradAKT=p(35);
kb1=p(36); k43b1=p(37);
kKSRphos=p(38); kKSRdephos=p(39);
kMekByBraf=p(40); kMekByCraf=p(41); kMekByKSR=p(42);
Tram=p(43); K_tram_RAF=p(44); K_tram_KSR=p(45); n_tram=p(46);

% integrate exactly at the experimental times (seconds)
[T, Y] = ode45(@(t,y) Mapk( ...
    t,y, ...
    ka1,kr1,kc1,kpCraf,kpMek,kpErk, ...
    kDegradEgfr,kErkInbEgfr,kShcDephos,kptpDeg,kGrb2CombShc, ...
    kSprtyInbGrb2,kSosCombGrb2,kErkPhosSos, ...
    kErkPhosPcraf,kPcrafDegrad,kErkPhosMek,kMekDegrad, ...
    kDuspInbErk,kErkDeg,kinbBraf,kDuspStop,kDusps, ...
    kSproutyForm,kSprtyComeDown,kdegrad,km_Sprty_decay,km_Dusp,km_Sprty, ...
    kErkDephos,kDuspDeg, ...
    kEGFRPI3k,kMTOR_Feedback, ...
    kAkt,kdegradAKT, ...
    kb1,k43b1, ...
    kKSRphos,kKSRdephos, ...
    kMekByBraf,kMekByCraf,kMekByKSR, ...
    Tram,K_tram_RAF,K_tram_KSR,n_tram), ...
    timeStamps*3600, y0);

% pull model outputs at those times (already at t = timeStamps*3600)
m_pEGFR = Y(:,3);
m_RASgtp = Y(:,14);
m_pCRAF = Y(:,17);
m_pMEK  = Y(:,21);
m_pERK  = Y(:,23);
m_ERK   = Y(:,22);
m_DUSP  = Y(:,24);
m_pAKT  = Y(:,40);

% (Optional) light normalization if your data are normalized:
norm = @(v,mode) v ./ max(eps, (mode=="first")*v(1) + (mode=="last")*v(end) + (mode=="max")*max(v));
m_pEGFR = norm(m_pEGFR,'max');
m_RASgtp= norm(m_RASgtp,'last');
m_pMEK  = norm(m_pMEK,'first');
m_pERK  = norm(m_pERK,'first');
m_DUSP  = norm(m_DUSP,'max');
m_pCRAF = norm(m_pCRAF,'max');
m_pAKT  = norm(m_pAKT,'max');
m_ERK   = norm(m_ERK,'max');

% squared-error across series (weights can be tuned)
w = struct('EGFR',1,'RAS',1,'MEK',1,'ERK',1,'DUSP',1,'CRAF',1,'AKT',1);

err = 0;
err = err + w.EGFR * sum((m_pEGFR - expPEGFR(:)).^2);
err = err + w.RAS  * sum((m_RASgtp - RASGtpValsExp(:)).^2);
err = err + w.MEK  * sum((m_pMEK  - expMEK(:)).^2);
err = err + w.DUSP * sum((m_DUSP  - expDusp(:)).^2);
err = err + w.CRAF * sum((m_pCRAF - expCRAF(:)).^2);
err = err + w.AKT  * sum((m_pAKT  - expAKT(:)).^2);
end




%% --- Visualization: Model Fit on All Data ---

% Define colors and markers
data_color = [0.8, 0.2, 0.2];   % Red for experimental data
model_color = [0.2, 0.6, 0.2];  % Green for model fit

% All timepoints for smooth model curve
tFine_hr = linspace(0, 48, 400);
tFine_sec = tFine_hr * 3600;

% Simulate with optimized parameters for smooth curve
[T_fine, Y_fine] = ode45(@(t,y) Mapk( ...
    t,y, ...
    optimizedParams(1),optimizedParams(2),optimizedParams(3), ...
    optimizedParams(4),optimizedParams(5),optimizedParams(6), ...
    optimizedParams(7),optimizedParams(8),optimizedParams(9),optimizedParams(10),optimizedParams(11), ...
    optimizedParams(12),optimizedParams(13),optimizedParams(14), ...
    optimizedParams(15),optimizedParams(16),optimizedParams(17),optimizedParams(18), ...
    optimizedParams(19),optimizedParams(20),optimizedParams(21),optimizedParams(22),optimizedParams(23), ...
    optimizedParams(24),optimizedParams(25),optimizedParams(26),optimizedParams(27),optimizedParams(28),optimizedParams(29), ...
    optimizedParams(30),optimizedParams(31), ...
    optimizedParams(32),optimizedParams(33), ...
    optimizedParams(34),optimizedParams(35), ...
    optimizedParams(36),optimizedParams(37), ...
    optimizedParams(38),optimizedParams(39),optimizedParams(40), ...
    optimizedParams(41),optimizedParams(42),optimizedParams(43),optimizedParams(44),optimizedParams(45),optimizedParams(46)), ...
    tFine_sec, y0);

% Normalize smooth model outputs
m_pEGFR_smooth = normit(Y_fine(:,3),'max');
m_RASgtp_smooth = normit(Y_fine(:,14),'last');
m_pCRAF_smooth = normit(Y_fine(:,17),'max');
m_pMEK_smooth = normit(Y_fine(:,21),'first');
m_pERK_smooth = normit(Y_fine(:,23),'first');
m_DUSP_smooth = normit(Y_fine(:,24),'max');
m_pAKT_smooth = normit(Y_fine(:,40),'max');

% Plot each protein species separately

% Plot 1: pEGFR
figure('Name','pEGFR Model Fit');
plot(tFine_hr, m_pEGFR_smooth, '-', 'Color', model_color, 'LineWidth', 3);
hold on;
plot(timeStamps_all, expPEGFR_all, 'o', 'Color', data_color, 'MarkerSize', 10, 'MarkerFaceColor', data_color, 'LineWidth', 2);
xlabel('Time (hrs)', 'FontSize', 12); ylabel('pEGFR', 'FontSize', 12); 
title('pEGFR Model Fit (All Timestamps)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Model Fit', 'Experimental Data'}, 'Location', 'best', 'FontSize', 11);
grid on; xlim([0, 50]);
set(gca, 'FontSize', 11);

% Plot 2: RAS-GTP
figure('Name','RAS-GTP Model Fit');
plot(tFine_hr, m_RASgtp_smooth, '-', 'Color', model_color, 'LineWidth', 3);
hold on;
plot(timeStamps_all, RASGtpValsExp_all, 'o', 'Color', data_color, 'MarkerSize', 10, 'MarkerFaceColor', data_color, 'LineWidth', 2);
xlabel('Time (hrs)', 'FontSize', 12); ylabel('RAS-GTP', 'FontSize', 12); 
title('RAS-GTP Model Fit (All Timestamps)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Model Fit', 'Experimental Data'}, 'Location', 'best', 'FontSize', 11);
grid on; xlim([0, 50]);
set(gca, 'FontSize', 11);

% Plot 3: pCRAF
figure('Name','pCRAF Model Fit');
plot(tFine_hr, m_pCRAF_smooth, '-', 'Color', model_color, 'LineWidth', 3);
hold on;
plot(timeStamps_all, expCRAF_all, 'o', 'Color', data_color, 'MarkerSize', 10, 'MarkerFaceColor', data_color, 'LineWidth', 2);
xlabel('Time (hrs)', 'FontSize', 12); ylabel('pCRAF', 'FontSize', 12); 
title('pCRAF Model Fit (All Timestamps)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Model Fit', 'Experimental Data'}, 'Location', 'best', 'FontSize', 11);
grid on; xlim([0, 50]);
set(gca, 'FontSize', 11);

% Plot 4: pMEK
figure('Name','pMEK Model Fit');
plot(tFine_hr, m_pMEK_smooth, '-', 'Color', model_color, 'LineWidth', 3);
hold on;
plot(timeStamps_all, expMEK_all, 'o', 'Color', data_color, 'MarkerSize', 10, 'MarkerFaceColor', data_color, 'LineWidth', 2);
xlabel('Time (hrs)', 'FontSize', 12); ylabel('pMEK', 'FontSize', 12); 
title('pMEK Model Fit (All Timestamps)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Model Fit', 'Experimental Data'}, 'Location', 'best', 'FontSize', 11);
grid on; xlim([0, 50]);
set(gca, 'FontSize', 11);

% Plot 5: pERK
figure('Name','pERK Model Fit');
plot(tFine_hr, m_pERK_smooth, '-', 'Color', model_color, 'LineWidth', 3);
hold on;
plot(timeStamps_all, PexpErkVals_all, 'o', 'Color', data_color, 'MarkerSize', 10, 'MarkerFaceColor', data_color, 'LineWidth', 2);
xlabel('Time (hrs)', 'FontSize', 12); ylabel('pERK', 'FontSize', 12); 
title('pERK Model Fit (All Timestamps)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Model Fit', 'Experimental Data'}, 'Location', 'best', 'FontSize', 11);
grid on; xlim([0, 50]);
set(gca, 'FontSize', 11);

% Plot 6: DUSP
figure('Name','DUSP Model Fit');
plot(tFine_hr, m_DUSP_smooth, '-', 'Color', model_color, 'LineWidth', 3);
hold on;
plot(timeStamps_all, expDusp_all, 'o', 'Color', data_color, 'MarkerSize', 10, 'MarkerFaceColor', data_color, 'LineWidth', 2);
xlabel('Time (hrs)', 'FontSize', 12); ylabel('DUSP', 'FontSize', 12); 
title('DUSP Model Fit (All Timestamps)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Model Fit', 'Experimental Data'}, 'Location', 'best', 'FontSize', 11);
grid on; xlim([0, 50]);
set(gca, 'FontSize', 11);

% Plot 7: pAKT
figure('Name','pAKT Model Fit');
plot(tFine_hr, m_pAKT_smooth, '-', 'Color', model_color, 'LineWidth', 3);
hold on;
plot(timeStamps_all, expAKT_all, 'o', 'Color', data_color, 'MarkerSize', 10, 'MarkerFaceColor', data_color, 'LineWidth', 2);
xlabel('Time (hrs)', 'FontSize', 12); ylabel('pAKT', 'FontSize', 12); 
title('pAKT Model Fit (All Timestamps)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Model Fit', 'Experimental Data'}, 'Location', 'best', 'FontSize', 11);
grid on; xlim([0, 50]);
set(gca, 'FontSize', 11);

% Plot 8: Fit Error Summary
figure('Name','Model Fit Error Summary');
% Calculate total errors for each protein (sum across all timepoints)
total_errors = [sum(fit_errors.pEGFR), sum(fit_errors.RAS), sum(fit_errors.pCRAF), ...
                sum(fit_errors.pMEK), sum(fit_errors.pERK), sum(fit_errors.DUSP), sum(fit_errors.pAKT)];

protein_names = {'pEGFR', 'RAS', 'pCRAF', 'pMEK', 'pERK', 'DUSP', 'pAKT'};
x_pos = 1:length(protein_names);

% Create bar plot
bar(x_pos, total_errors, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k');
set(gca, 'XTick', x_pos, 'XTickLabel', protein_names);
ylabel('Total Fit Error (sum across all timepoints)', 'FontSize', 12); 
title('Model Fit Errors by Protein (All Timestamps)', 'FontSize', 14, 'FontWeight', 'bold');
xtickangle(45); grid on;
set(gca, 'FontSize', 11);

% Add overall summary text
text(0.02, 0.98, sprintf('Total Fit Error: %.4f\nModel fitted on all timestamps (0-48hrs)', errorOpt), ...
     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', 'white', 'EdgeColor', 'black');
