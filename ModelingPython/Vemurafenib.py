
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, differential_evolution
from numba import jit
import time
import os

# =============================================================================
# MAPK/PI3K PATHWAY MODEL WITH VEMURAFENIB AND PARADOXICAL ACTIVATION
# Converted from MATLAB to Python
# =============================================================================

def print_header():
    print("═══════════════════════════════════════════════════════════════════════════")
    print("   MAPK/PI3K PATHWAY MODEL WITH VEMURAFENIB (PYTHON VERSION)")
    print("   Paradoxical Activation Included")
    print("═══════════════════════════════════════════════════════════════════════════\n")

# =============================================================================
# ODE FUNCTION
# =============================================================================
@jit(nopython=True)
def mapk_ode(t, y, p):
    """
    MAPK/PI3K Pathway ODE System with Vemurafenib and Paradoxical Activation.
    """
    dydt = np.zeros_like(y)
    
    # Unpack parameters (p is 0-indexed in Python)
    ka1 = p[0];  kr1 = p[1];  kc1 = p[2]
    kpCraf = p[3]; kpMek = p[4]; kpErk = p[5]
    kDegradEgfr = p[6]; kErkInbEgfr = p[7]; kShcDephos = p[8]; kptpDeg = p[9]; kGrb2CombShc = p[10]
    kSprtyInbGrb2 = p[11]; kSosCombGrb2 = p[12]; kErkPhosSos = p[13]
    kErkPhosPcraf = p[14]; kPcrafDegrad = p[15]; kErkPhosMek = p[16]; kMekDegrad = p[17]
    kDuspInbErk = p[18]; kErkDeg = p[19]; kinbBraf = p[20]; kDuspStop = p[21]; kDusps = p[22]
    kSproutyForm = p[23]; kSprtyComeDown = p[24]; kdegrad = p[25]; km_Sprty_decay = p[26]; km_Dusp = p[27]; km_Sprty = p[28]
    kErkDephos = p[29]; kDuspDeg = p[30]
    kHer2_act = p[31]; kHer3_act = p[32]
    k_p85_bind_EGFR = p[33]; k_p85_bind_Her2 = p[34]; k_p85_bind_Her3 = p[35]; k_p85_bind_IGFR = p[36]
    k_p85_unbind = p[37]; k_PI3K_recruit = p[38]
    kMTOR_Feedback = p[39]
    k_PIP2_to_PIP3 = p[40]; k_PTEN = p[41]
    kAkt = p[42]; kdegradAKT = p[43]
    kb1 = p[44]; k43b1 = p[45]; k4ebp1 = p[46]; k_4EBP1_dephos = p[47]
    kKSRphos = p[48]; kKSRdephos = p[49]
    kMekByBraf = p[50]; kMekByCraf = p[51]; kMekByKSR = p[52]
    Tram = p[53]; K_tram_RAF = p[54]; K_tram_KSR = p[55]; n_tram = p[56]
    Vemurafenib = p[57]; kDimerForm = p[58]; kDimerDissoc = p[59]
    kParadoxCRAF = p[60]; IC50_vem = p[61]; Hill_n_vem = p[62]
    
    EPS = 1e-16

    # ========================================================================
    # MODULE 1: RTK SIGNALING
    # ========================================================================
    # EGFR (y[0]-y[2])
    dydt[0] = -ka1*y[0] + kr1*y[1]                                    # EGFR
    dydt[1] = ka1*y[0] - kr1*y[1] - kc1*y[1]                          # EGFR:ligand
    dydt[2] = kc1*y[1] - kDegradEgfr*y[2] - kErkInbEgfr*y[28]*y[2]    # pEGFR
    
    # Her2 (y[3]-y[5])
    dydt[3] = -kHer2_act*y[3] + kr1*y[4]                              # Her2
    dydt[4] = kHer2_act*y[3] - kr1*y[4] - kc1*y[4]                    # Her2:ligand
    dydt[5] = kc1*y[4] - kDegradEgfr*y[5] - kErkInbEgfr*y[28]*y[5]    # pHer2
    
    # Her3 (y[6]-y[8])
    dydt[6] = -kHer3_act*y[6] + kr1*y[7]                              # Her3
    dydt[7] = kHer3_act*y[6] - kr1*y[7] - kc1*y[7]                    # Her3:ligand
    dydt[8] = kc1*y[7] - kDegradEgfr*y[8] - kErkInbEgfr*y[28]*y[8]    # pHer3

    # ========================================================================
    # MODULE 2: Shc/Grb2/SOS
    # ========================================================================
    dydt[9]  = -ka1*y[2]*y[9]                                         # Shc
    dydt[10] = ka1*y[2]*y[9] - kShcDephos*y[11]*y[10]                 # pEGFR:Shc
    dydt[11] = -kptpDeg*y[10]*y[11]                                   # pShc
    
    dydt[12] = kGrb2CombShc*y[10]*y[2] - kSprtyInbGrb2*y[26]*y[12]    # Grb2
    dydt[13] = kSosCombGrb2*y[12]*y[10] - kErkPhosSos*y[24]*y[13]     # pShc:Grb2:SOS

    # ========================================================================
    # MODULE 3: RAS ACTIVATION
    # ========================================================================
    dydt[14] = -ka1*y[13]*y[14]                                       # HRAS-GDP
    dydt[15] = ka1*y[13]*y[14]                                        # HRAS-GTP
    dydt[16] = -ka1*y[13]*y[16]                                       # NRAS-GDP
    dydt[17] = ka1*y[13]*y[16]                                        # NRAS-GTP
    dydt[18] = -ka1*y[13]*y[18]                                       # KRAS-GDP
    dydt[19] = ka1*y[13]*y[18] - ka1*y[19]*y[20]                      # KRAS-GTP
    dydt[20] = -ka1*y[19]*y[20]                                       # KRAS intermediate
    
    # ========================================================================
    # MODULE 4: RAF SIGNALING (Paradoxical)
    # ========================================================================
    # Vemurafenib inhibition
    IC50_n = IC50_vem**Hill_n_vem
    Vem_n = Vemurafenib**Hill_n_vem
    kBRAF_eff = ka1 * IC50_n / (IC50_n + Vem_n + EPS)
    
    paradox_activation = kParadoxCRAF * Vemurafenib * y[61]  # y[61] = dimer
    
    # CRAF
    dydt[21] = -kpCraf*y[19]*y[21] + kErkPhosPcraf*y[28]*y[22] + kPcrafDegrad*y[22]*y[35] \
               - kDimerForm*y[24]*y[21]*Vemurafenib + kDimerDissoc*y[61]
    dydt[22] = kpCraf*y[19]*y[21] - kErkPhosPcraf*y[28]*y[22] - kPcrafDegrad*y[22]*y[35] \
               + paradox_activation                                   # pCRAF
    
    # BRAF
    dydt[23] = -kBRAF_eff*y[23] - kDimerForm*y[24]*y[21]*Vemurafenib + kDimerDissoc*y[61]
    dydt[24] = kBRAF_eff*y[23] - kinbBraf*y[24] - kDimerForm*y[24]*y[21]*Vemurafenib + kDimerDissoc*y[61]
    
    # Dimer
    dydt[61] = kDimerForm*y[24]*y[21]*Vemurafenib - kDimerDissoc*y[61] - kPcrafDegrad*y[61]*y[35]

    # ========================================================================
    # MODULE 5: MEK
    # ========================================================================
    raf_to_mek = (kpMek*y[22] + kMekByBraf*y[24] + kMekByCraf*y[22])
    ksr_to_mek = (kMekByKSR * y[60])
    
    dydt[25] = -(raf_to_mek + ksr_to_mek)*y[25] + kErkPhosMek*y[28]*y[26] + kMekDegrad*y[26]*y[34]
    dydt[26] = (raf_to_mek + ksr_to_mek)*y[25] - kErkPhosMek*y[28]*y[26] - kMekDegrad*y[26]*y[34]

    # ========================================================================
    # MODULE 6: ERK
    # ========================================================================
    dydt[27] = -kpErk*y[26]*y[27] + kDuspInbErk*y[30]*y[28] + kErkDeg*y[28]*y[33] + kErkDephos*y[30]*y[28]
    dydt[28] = kpErk*y[26]*y[27] - kDuspInbErk*y[30]*y[28] - kErkDeg*y[28]*y[33] - kErkDephos*y[30]*y[28]

    # ========================================================================
    # FEEDBACK TRACKERS
    # ========================================================================
    dydt[29] = km_Dusp*y[28]/(1 + (km_Dusp/kDusps)*y[28]) - kDuspStop*y[29]*y[36] - kDuspDeg*y[29]*y[28] # mDUSP
    dydt[30] = -kDuspStop*y[29]*y[30] # DUSP
    
    dydt[31] = km_Sprty*y[28]/(1 + (km_Sprty/kSproutyForm)*y[28]) - kSprtyComeDown*y[31]*y[32] # mSPRY
    dydt[32] = -kSprtyComeDown*y[31]*y[32] # SPRY
    
    dydt[33] = -kErkDeg*y[28]*y[33]      # pERK degrad
    dydt[34] = -kMekDegrad*y[26]*y[34]   # pMEK degrad
    dydt[35] = -kPcrafDegrad*y[22]*y[35] # pCRAF degrad
    dydt[36] = -kDuspStop*y[29]*y[36]    # DUSP stop

    # ========================================================================
    # MODULE 10: PI3K/AKT/mTOR
    # ========================================================================
    dydt[37] = -ka1*y[37] + kr1*y[38]
    dydt[38] = ka1*y[37] - kr1*y[38] - kc1*y[38]
    dydt[39] = kc1*y[38] - kErkInbEgfr*y[28]*y[39]
    
    dydt[40] = -ka1*y[2]*y[40]
    dydt[41] = ka1*y[2]*y[40]
    
    dydt[42] = -k_p85_bind_EGFR*y[2]*y[42] \
               - k_p85_bind_Her2*y[5]*y[42] \
               - k_p85_bind_Her3*y[8]*y[42] \
               - k_p85_bind_IGFR*y[39]*y[42] \
               + k_p85_unbind*(y[43] + y[44] + y[45] + y[46])

    dydt[43] = k_p85_bind_EGFR*y[2]*y[42] - k_p85_unbind*y[43]  # p85:EGFR
    dydt[44] = k_p85_bind_Her2*y[5]*y[42] - k_p85_unbind*y[44]  # p85:Her2
    dydt[45] = k_p85_bind_Her3*y[8]*y[42] - k_p85_unbind*y[45]  # p85:Her3
    dydt[46] = k_p85_bind_IGFR*y[39]*y[42] - k_p85_unbind*y[46] # p85:IGFR
    
    total_p85_RTK = y[43] + y[44] + y[45] + y[46]
    
    dydt[47] = -k_PI3K_recruit * total_p85_RTK * y[47] + kMTOR_Feedback * y[55] * y[48] # PI3K
    dydt[48] = k_PI3K_recruit * total_p85_RTK * y[47] - kMTOR_Feedback * y[55] * y[48]  # PI3K*
    
    dydt[49] = -k_PIP2_to_PIP3 * y[48] * y[49] + k_PTEN * y[50] # PIP2
    dydt[50] = k_PIP2_to_PIP3 * y[48] * y[49] - k_PTEN * y[50]  # PIP3
    
    dydt[51] = -kAkt * y[50] * y[51] + kdegradAKT * y[52] # AKT
    dydt[52] = kAkt * y[50] * y[51] - kdegradAKT * y[52]  # pAKT
    
    dydt[53] = (1 - y[52]) * ka1 / (1 + (y[53] / 15e-5))   # FOXO
    
    dydt[54] = -ka1 * y[52] * y[54] + kdegrad * y[55]      # mTORC
    dydt[55] = ka1 * y[52] * y[54] - kdegrad * y[55]       # mTORC*
    
    dydt[56] = -k4ebp1 * y[55] * y[56] + kb1 * y[57] + k_4EBP1_dephos * y[58] # 4EBP1
    dydt[57] = k4ebp1 * y[55] * y[56] - kb1 * y[57] - k43b1 * y[57]           # inter
    dydt[58] = k43b1 * y[57] - k_4EBP1_dephos * y[58]                         # p4EBP1
    
    # ========================================================================
    # KSR
    # ========================================================================
    dydt[59] = -kKSRphos * (y[15] + y[17] + y[20]) * y[59] + kKSRdephos * y[60]
    dydt[60] = kKSRphos * (y[15] + y[17] + y[20]) * y[59] - kKSRdephos * y[60]
    
    return dydt

# =============================================================================
# NORMALIZATION
# =============================================================================
def normit(v, mode):
    arr = np.array(v)
    if mode == 'first':
        denom = arr[0]
    elif mode == 'last':
        denom = arr[-1]
    elif mode == 'max':
        denom = np.max(arr)
    else:
        denom = 1.0
    return arr / max(1e-16, denom)

# =============================================================================
# OBJECTIVE FUNCTION
# =============================================================================
def objective_function(p, time_seconds, exp_data_norm, y0):
    sol = solve_ivp(mapk_ode, [time_seconds[0], time_seconds[-1]], y0, 
                    args=(p,), t_eval=time_seconds, method='Radau', rtol=1e-5, atol=1e-8)
    
    if not sol.success or sol.y.shape[1] != len(time_seconds):
        return 1e9

    Y = sol.y.T 
    m_pEGFR = Y[:, 2]; m_pCRAF = Y[:, 22]; m_pMEK = Y[:, 26]
    m_pERK = Y[:, 28]; m_DUSP = Y[:, 30]; m_pAKT = Y[:, 52]; m_p4EBP1 = Y[:, 58]
    
    m_pEGFR_n = normit(m_pEGFR, 'max')
    m_pMEK_n  = normit(m_pMEK, 'first')
    m_pERK_n  = normit(m_pERK, 'first')
    m_DUSP_n  = normit(m_DUSP, 'max')
    m_pCRAF_n = normit(m_pCRAF, 'max')
    m_pAKT_n  = normit(m_pAKT, 'max')
    m_p4EBP1_n= normit(m_p4EBP1, 'max')
    
    err = 0
    weights = {'EGFR':1, 'MEK':1, 'ERK':3, 'DUSP':3, 'CRAF':1, 'AKT':1, 'p4EBP1':1}
    
    err += weights['EGFR'] * np.sum((m_pEGFR_n - exp_data_norm['pEGFR'])**2)
    err += weights['MEK']  * np.sum((m_pMEK_n  - exp_data_norm['pMEK'])**2)
    err += weights['ERK']  * np.sum((m_pERK_n  - exp_data_norm['pERK'])**2)
    err += weights['DUSP'] * np.sum((m_DUSP_n  - exp_data_norm['DUSP'])**2)
    err += weights['CRAF'] * np.sum((m_pCRAF_n - exp_data_norm['pCRAF'])**2)
    err += weights['AKT']  * np.sum((m_pAKT_n  - exp_data_norm['pAKT'])**2)
    err += weights['p4EBP1']*np.sum((m_p4EBP1_n- exp_data_norm['p4EBP1'])**2)
    
    return err

# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print_header()
    
    # Initial P vector
    p_init = np.array([
        10e-5, 10e-6, 10e-3, 10e-4, 7e-6, 20e-5, 10e-7, 20e-6, 10e-5, 7e-7, 36e-6, 
        10e-5, 10e-5, 17e-6, 15.7e-6, 1.7e-4, 20.7e-5, 1.7e-4, 2.7e-4, 1.7e-7, 10e-2, 1.7e-4, 1.11e-4, 
        1.7e-06, 5.5e-5, 10e-7, 8e-5, 1.11e-4, 1.7e-06, 1e-5, 1e-6, 
        10e-5, 10e-5, 10e-4, 10e-4, 10e-4, 10e-4, 10e-5, 10e-4, 10e-4, 
        10e-4, 10e-5, 10e-5, 10e-7, 10e-8, 10e-3, 10e-5, 10e-5, 5e-6, 5e-6, 
        8e-6, 8e-6, 3e-6, 1e-6, 5e-1, 5e-12, 2, 1.0, 6e-6, 1e-5, 0.5, 0.4, 1.5
    ])

    # Initial Conditions (y0)
    y0 = np.zeros(62)
    y0[0]=1.0; y0[3]=1.0; y0[6]=1.0; y0[9]=1.0; y0[11]=1.0
    y0[14]=0.0; y0[16]=0.0; y0[18]=1.0; y0[20]=1.0
    y0[21]=0.8; y0[22]=0.2; y0[23]=0.0; y0[24]=1.0
    y0[25]=1.0; y0[26]=1.0; y0[27]=1.0; y0[28]=0.8
    y0[29]=1.0; y0[30]=1.0; y0[31]=1.0; y0[32]=1.0
    y0[33]=1.0; y0[34]=1.0; y0[35]=1.0; y0[36]=1.0
    y0[37]=1.0; y0[40]=1.0; y0[42]=1.0
    y0[47]=1.0; y0[49]=1.0; y0[51]=1.0; y0[54]=1.0
    y0[56]=1.0; y0[59]=1.0; y0[61]=0.0

    # Data
    time_hours = np.array([0, 1, 4, 8, 24, 48])
    time_seconds = time_hours * 3600
    exp_data = {
        'pMEK': np.array([1.75938884, 0.170160085, 0.095112609, 0.201000276, 0.219207054, 0.502831668]),
        'pERK': np.array([2.903209735, 0.207867788, 0.303586121, 0.805254439, 1.408362153, 1.847606441]),
        'DUSP': np.array([2.677161325, 2.782754577, 1.130758062, 0.395642757, 0.828575853, 0.916618219]),
        'pEGFR': np.array([0.291928893, 0.392400458, 0.265016688, 0.394238749, 0.006158316, 0.008115099]),
        'pCRAF': np.array([0.366397596, 0.537106733, 0.465541704, 0.586732657, 1.102322681, 0.269181259]),
        'pAKT': np.array([0.513544148, 0.613178403, 1.03451863, 1.113391047, 0.535242724, 0.538273551]),
        'p4EBP1': np.array([1.002468056,1.276793699,1.252681407,1.707504483,1.271216967,0.61389625])
    }
    exp_data_norm = {k: (v-np.min(v))/(np.max(v)-np.min(v)+1e-16) for k,v in exp_data.items()}

    # Bounds
    bnds = [(1e-12, 1e-3) for _ in range(len(p_init))]
    bnds[57] = (0.0, 1.0); bnds[61] = (0.01, 0.99); bnds[62] = (0.5, 5.0)

    print("Starting optimization (Parallel Differential Evolution using JIT)...")
    t0 = time.time()
    
    # Using Differential Evolution for parallel implementation
    # workers=-1 uses all available CPU cores
    # updating='deferred' is best for parallelization
    res = differential_evolution(objective_function, bnds, args=(time_seconds, exp_data_norm, y0),
                                 strategy='best1bin', maxiter=50, popsize=10, 
                                 mutation=(0.5, 1), recombination=0.7,
                                 tol=0.01, workers=-1, updating='deferred', disp=True)
                                 
    print(f"Done in {time.time()-t0:.1f}s. Err: {res.fun:.4f}")
    print("Optimized Parameters:", res.x)

    # Plot
    t_fine = np.linspace(0, 48*3600, 400)
    sol = solve_ivp(mapk_ode, [t_fine[0], t_fine[-1]], y0, args=(res.x,), t_eval=t_fine, method='Radau')
    
    # If the solver failed to reach the end or t_eval was ignored due to failure,
    # sol.t and sol.y will still be consistent with each other.
    # We use sol.t for the x-axis to guarantee dimension match.
    
    plt.figure()
    if sol.success:
        print("Integration successful for plotting.")
    else:
        print(f"Integration failed or truncated during plotting: {sol.message}")
        
    # Using sol.t ensures x and y dimensions match (e.g. if solver stopped early)
    plot_t = sol.t / 3600
    plot_y = normit(sol.y[28], 'first')
    
    plt.plot(plot_t, plot_y, label='Model pERK')
    plt.plot(time_hours, exp_data_norm['pERK'], 'ro', label='Data')
    plt.legend()
    plt.title('pERK Fit')
    plt.xlabel('Time (h)')
    plt.savefig('/Users/bharadwajanandivada/SCMPA/Modelingpython/Model_Fit.png')
