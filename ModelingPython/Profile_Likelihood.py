
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from numba import jit
import time
import os

# =============================================================================
# 1. MODEL DEFINITION (Must match Vemurafenib.py)
# =============================================================================
@jit(nopython=True)
def mapk_ode(t, y, p):
    dydt = np.zeros_like(y)
    
    # Unpack parameters
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

    # MODULE 1: RTK SIGNALING
    dydt[0] = -ka1*y[0] + kr1*y[1]                                    # EGFR
    dydt[1] = ka1*y[0] - kr1*y[1] - kc1*y[1]                          # EGFR:ligand
    dydt[2] = kc1*y[1] - kDegradEgfr*y[2] - kErkInbEgfr*y[28]*y[2]    # pEGFR
    
    dydt[3] = -kHer2_act*y[3] + kr1*y[4]                              # Her2
    dydt[4] = kHer2_act*y[3] - kr1*y[4] - kc1*y[4]                    # Her2:ligand
    dydt[5] = kc1*y[4] - kDegradEgfr*y[5] - kErkInbEgfr*y[28]*y[5]    # pHer2
    
    dydt[6] = -kHer3_act*y[6] + kr1*y[7]                              # Her3
    dydt[7] = kHer3_act*y[6] - kr1*y[7] - kc1*y[7]                    # Her3:ligand
    dydt[8] = kc1*y[7] - kDegradEgfr*y[8] - kErkInbEgfr*y[28]*y[8]    # pHer3

    # MODULE 2: Shc/Grb2/SOS
    dydt[9]  = -ka1*y[2]*y[9]                                         # Shc
    dydt[10] = ka1*y[2]*y[9] - kShcDephos*y[11]*y[10]                 # pEGFR:Shc
    dydt[11] = -kptpDeg*y[10]*y[11]                                   # pShc
    
    dydt[12] = kGrb2CombShc*y[10]*y[2] - kSprtyInbGrb2*y[26]*y[12]    # Grb2
    dydt[13] = kSosCombGrb2*y[12]*y[10] - kErkPhosSos*y[24]*y[13]     # pShc:Grb2:SOS

    # MODULE 3: RAS ACTIVATION
    dydt[14] = -ka1*y[13]*y[14]                                       # HRAS-GDP
    dydt[15] = ka1*y[13]*y[14]                                        # HRAS-GTP
    dydt[16] = -ka1*y[13]*y[16]                                       # NRAS-GDP
    dydt[17] = ka1*y[13]*y[16]                                        # NRAS-GTP
    dydt[18] = -ka1*y[13]*y[18]                                       # KRAS-GDP
    dydt[19] = ka1*y[13]*y[18] - ka1*y[19]*y[20]                      # KRAS-GTP
    dydt[20] = -ka1*y[19]*y[20]                                       # KRAS intermediate
    
    # MODULE 4: RAF SIGNALING
    IC50_n = IC50_vem**Hill_n_vem
    Vem_n = Vemurafenib**Hill_n_vem
    kBRAF_eff = ka1 * IC50_n / (IC50_n + Vem_n + EPS)
    paradox_activation = kParadoxCRAF * Vemurafenib * y[61]
    
    dydt[21] = -kpCraf*y[19]*y[21] + kErkPhosPcraf*y[28]*y[22] + kPcrafDegrad*y[22]*y[35] \
               - kDimerForm*y[24]*y[21]*Vemurafenib + kDimerDissoc*y[61]
    dydt[22] = kpCraf*y[19]*y[21] - kErkPhosPcraf*y[28]*y[22] - kPcrafDegrad*y[22]*y[35] \
               + paradox_activation
    dydt[23] = -kBRAF_eff*y[23] - kDimerForm*y[24]*y[21]*Vemurafenib + kDimerDissoc*y[61]
    dydt[24] = kBRAF_eff*y[23] - kinbBraf*y[24] - kDimerForm*y[24]*y[21]*Vemurafenib + kDimerDissoc*y[61]
    dydt[61] = kDimerForm*y[24]*y[21]*Vemurafenib - kDimerDissoc*y[61] - kPcrafDegrad*y[61]*y[35]

    # MODULE 5: MEK
    raf_to_mek = (kpMek*y[22] + kMekByBraf*y[24] + kMekByCraf*y[22])
    ksr_to_mek = (kMekByKSR * y[60])
    dydt[25] = -(raf_to_mek + ksr_to_mek)*y[25] + kErkPhosMek*y[28]*y[26] + kMekDegrad*y[26]*y[34]
    dydt[26] = (raf_to_mek + ksr_to_mek)*y[25] - kErkPhosMek*y[28]*y[26] - kMekDegrad*y[26]*y[34]

    # MODULE 6: ERK
    dydt[27] = -kpErk*y[26]*y[27] + kDuspInbErk*y[30]*y[28] + kErkDeg*y[28]*y[33] + kErkDephos*y[30]*y[28]
    dydt[28] = kpErk*y[26]*y[27] - kDuspInbErk*y[30]*y[28] - kErkDeg*y[28]*y[33] - kErkDephos*y[30]*y[28]

    # FEEDBACK
    dydt[29] = km_Dusp*y[28]/(1 + (km_Dusp/kDusps)*y[28]) - kDuspStop*y[29]*y[36] - kDuspDeg*y[29]*y[28]
    dydt[30] = -kDuspStop*y[29]*y[30]
    dydt[31] = km_Sprty*y[28]/(1 + (km_Sprty/kSproutyForm)*y[28]) - kSprtyComeDown*y[31]*y[32]
    dydt[32] = -kSprtyComeDown*y[31]*y[32]
    dydt[33] = -kErkDeg*y[28]*y[33]; dydt[34] = -kMekDegrad*y[26]*y[34]
    dydt[35] = -kPcrafDegrad*y[22]*y[35]; dydt[36] = -kDuspStop*y[29]*y[36]

    # MODULE 10: PI3K/AKT/mTOR
    dydt[37] = -ka1*y[37] + kr1*y[38]
    dydt[38] = ka1*y[37] - kr1*y[38] - kc1*y[38]
    dydt[39] = kc1*y[38] - kErkInbEgfr*y[28]*y[39]
    dydt[40] = -ka1*y[2]*y[40]; dydt[41] = ka1*y[2]*y[40]
    
    dydt[42] = -k_p85_bind_EGFR*y[2]*y[42] - k_p85_bind_Her2*y[5]*y[42] \
               - k_p85_bind_Her3*y[8]*y[42] - k_p85_bind_IGFR*y[39]*y[42] \
               + k_p85_unbind*(y[43] + y[44] + y[45] + y[46])

    dydt[43] = k_p85_bind_EGFR*y[2]*y[42] - k_p85_unbind*y[43]
    dydt[44] = k_p85_bind_Her2*y[5]*y[42] - k_p85_unbind*y[44]
    dydt[45] = k_p85_bind_Her3*y[8]*y[42] - k_p85_unbind*y[45]
    dydt[46] = k_p85_bind_IGFR*y[39]*y[42] - k_p85_unbind*y[46]
    
    total_p85_RTK = y[43] + y[44] + y[45] + y[46]
    dydt[47] = -k_PI3K_recruit * total_p85_RTK * y[47] + kMTOR_Feedback * y[55] * y[48]
    dydt[48] = k_PI3K_recruit * total_p85_RTK * y[47] - kMTOR_Feedback * y[55] * y[48]
    dydt[49] = -k_PIP2_to_PIP3 * y[48] * y[49] + k_PTEN * y[50]
    dydt[50] = k_PIP2_to_PIP3 * y[48] * y[49] - k_PTEN * y[50]
    dydt[51] = -kAkt * y[50] * y[51] + kdegradAKT * y[52]
    dydt[52] = kAkt * y[50] * y[51] - kdegradAKT * y[52]
    dydt[53] = (1 - y[52]) * ka1 / (1 + (y[53] / 15e-5))
    dydt[54] = -ka1 * y[52] * y[54] + kdegrad * y[55]
    dydt[55] = ka1 * y[52] * y[54] - kdegrad * y[55]
    dydt[56] = -k4ebp1 * y[55] * y[56] + kb1 * y[57] + k_4EBP1_dephos * y[58]
    dydt[57] = k4ebp1 * y[55] * y[56] - kb1 * y[57] - k43b1 * y[57]
    dydt[58] = k43b1 * y[57] - k_4EBP1_dephos * y[58]
    
    # KSR
    dydt[59] = -kKSRphos * (y[15] + y[17] + y[20]) * y[59] + kKSRdephos * y[60]
    dydt[60] = kKSRphos * (y[15] + y[17] + y[20]) * y[59] - kKSRdephos * y[60]
    
    return dydt

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
def normit(v, mode):
    arr = np.array(v)
    if mode == 'first': denom = arr[0]
    elif mode == 'last': denom = arr[-1]
    elif mode == 'max': denom = np.max(arr)
    else: denom = 1.0
    return arr / max(1e-16, denom)

def objective_function_wrapped(p_resid, p_full_init, fix_idx, fix_val, time_seconds, exp_data_norm, y0):
    """
    Wrapper for objective function to handle fixed parameter
    p_resid: The parameters being optimized
    fix_idx: Index of parameter to fix
    fix_val: Value where it is fixed
    """
    # Reconstruct full parameter vector
    p_full = np.insert(p_resid, fix_idx, fix_val)
    
    # Run Solver
    # Using 'Radau' for robustness, but 'LSODA' might be faster for profile loops
    sol = solve_ivp(mapk_ode, [time_seconds[0], time_seconds[-1]], y0, 
                    args=(p_full,), t_eval=time_seconds, method='Radau', rtol=1e-4, atol=1e-6)
    
    if not sol.success or sol.y.shape[1] != len(time_seconds):
        return 1e9

    Y = sol.y.T 
    
    # Calculate Error
    m_pEGFR_n = normit(Y[:, 2], 'max')
    m_pMEK_n  = normit(Y[:, 26], 'first')
    m_pERK_n  = normit(Y[:, 28], 'first')
    m_DUSP_n  = normit(Y[:, 30], 'max')
    m_pCRAF_n = normit(Y[:, 22], 'max')
    m_pAKT_n  = normit(Y[:, 52], 'max')
    m_p4EBP1_n= normit(Y[:, 58], 'max')
    
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

def run_profile_likelihood(p_opt, param_idx, param_name, bounds, data_pack, steps=20, local_factor=0.2):
    """
    Calculates profile likelihood for a single parameter.
    local_factor: +/- X% to scan around the optimum
    """
    time_seconds, exp_data_norm, y0 = data_pack
    
    opt_val = p_opt[param_idx]
    
    # Define range (log space is usually better for biological parameters)
    lower = max(bounds[param_idx][0], opt_val * (1 - local_factor))
    upper = min(bounds[param_idx][1], opt_val * (1 + local_factor))
    
    # If the user's bounds are super tight (1e-12 to 1e-3), ensure we don't scan practically 0 range
    if upper < lower: upper = lower * 2
    
    scan_values = np.linspace(lower, upper, steps)
    chi2_values = []
    
    print(f"\nScanning {param_name} (idx {param_idx}) from {lower:.2e} to {upper:.2e}")
    
    # Prepare bounds for residual parameters (remove the fixed one)
    resid_bounds = bounds[:param_idx] + bounds[param_idx+1:]
    
    # Initial guess for residuals: use p_opt (minus the fixed one)
    p_resid_init = np.delete(p_opt, param_idx)

    for val in scan_values:
        res = minimize(objective_function_wrapped, p_resid_init, 
                      args=(p_opt, param_idx, val, time_seconds, exp_data_norm, y0),
                      method='SLSQP', bounds=resid_bounds, options={'maxiter': 50})
        chi2_values.append(res.fun)
        print(f"  {param_name}={val:.2e} -> Chi2={res.fun:.4f}")
        
    return scan_values, np.array(chi2_values)

# =============================================================================
# 3. MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    
    # --- SETUP DATA ---
    y0 = np.zeros(62)
    y0[0]=1.0; y0[3]=1.0; y0[6]=1.0; y0[9]=1.0; y0[11]=1.0; y0[18]=1.0; y0[20]=1.0
    y0[21]=0.8; y0[22]=0.2; y0[24]=1.0; y0[25]=1.0; y0[26]=1.0; y0[27]=1.0; y0[28]=0.8
    y0[29]=1.0; y0[30]=1.0; y0[31]=1.0; y0[32]=1.0; y0[33]=1.0; y0[34]=1.0; y0[35]=1.0; y0[36]=1.0
    y0[37]=1.0; y0[40]=1.0; y0[42]=1.0; y0[47]=1.0; y0[49]=1.0; y0[51]=1.0; y0[54]=1.0
    y0[56]=1.0; y0[59]=1.0; y0[61]=0.0

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
    data_pack = (time_seconds, exp_data_norm, y0)

    # --- INPUT OPTIMIZED PARAMETERS ---
    # NOTE: PLEASE UPDATE THIS IF YOU HAVE A NEWER/BETTER FIT
    # Using p_opt from previous successful run (~9.11 Error)
    p_opt = np.array([
        9.55502554, 7.39658806, 1.59626843, 0.90341034, 2.86945236, 1.71809685,
        1.65579974, 4.93931776, 8.12826981, 5.90168714, 4.68920853, 3.44885134,
        5.90111866, 4.42624268, 4.64513604, 9.05643887, 9.60998322, 9.56217679,
        2.35963287, 6.79951863, 5.87946535, 0.09250034, 1.90424759, 5.10442432,
        0.56169963, 3.18469232, 0.35158817, 0.01492954, 9.03274576, 2.23363394,
        9.31953107, 2.72173506, 5.36260683, 6.42877017, 8.46997403, 4.54857409,
        8.86568513, 9.44897184, 3.52400048, 4.43764163, 9.19995713, 0.35835363,
        6.19233967, 3.13842072, 9.42080297, 9.74279328, 5.22233274, 3.84020375,
        6.94715053, 9.64756984, 9.12178204, 5.43729891, 2.68886829, 7.16506273,
        1.70351907, 3.39321894, 5.44138224, 0.6710561,  9.86835713, 8.83050831,
        0.52693555, 0.68326328, 4.88175978
    ])
    
    bnds = [(1e-12, 10.0) for _ in range(len(p_opt))]
    
    # --- SELECT PARAMETERS TO PROFILE ---
    # We choose a few critical ones to demonstrate. 
    # Index 5: kpErk (Phosphorylation of ERK)
    # Index 57: Vemurafenib (Drug strength parameter, effectively)
    # Index 60: kParadoxCRAF (The paradoxical activation term)
    
    params_to_scan = [
        (5, 'kpErk (ERK Activation)'),
        (22, 'kDusps (Feedback)'),
        (60, 'kParadoxCRAF (Drug Paradox)')
    ]
    
    plt.figure(figsize=(15, 5))
    
    for i, (idx, name) in enumerate(params_to_scan):
        print(f"--- Profiling {name} ---")
        vals, chis = run_profile_likelihood(p_opt, idx, name, bnds, data_pack, steps=15, local_factor=0.5)
        
        # Plot
        ax = plt.subplot(1, 3, i+1)
        ax.plot(vals, chis, 'b.-', linewidth=2)
        ax.plot(p_opt[idx], np.min(chis), 'r*', markersize=15, label='Optimum')
        
        # Threshold for 95% CI (approx Chi2 with 1 dof = 3.84)
        threshold = np.min(chis) + 3.84
        ax.axhline(threshold, color='k', linestyle='--', label='95% Threshold')
        
        ax.set_title(f"Profile: {name}")
        ax.set_xlabel("Parameter Value")
        ax.set_ylabel("Chi-Squared Error")
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Profile_Likelihood_Analysis.png')
    print("\nAnalysis Complete. Saved to Profile_Likelihood_Analysis.png")
