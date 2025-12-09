
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# MODEL DEFINITION
# =============================================================================
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
# HELPER FUNCTIONS
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
# MAIN PLOTTING SCRIPT
# =============================================================================
if __name__ == "__main__":
    
    # 1. OPTIMIZED PARAMETERS (From optimization run)
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

    # 2. EXPERIMENTAL DATA
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
    
    # Normalize experimental data
    exp_data_norm = {k: (v-np.min(v))/(np.max(v)-np.min(v)+1e-16) for k,v in exp_data.items()}

    # 3. INITIAL CONDITIONS
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

    # 4. RUN SIMULATION
    print("Running simulation with optimized parameters...")
    t_fine = np.linspace(0, 48*3600, 400)
    sol = solve_ivp(mapk_ode, [t_fine[0], t_fine[-1]], y0, args=(p_opt,), t_eval=t_fine, method='Radau')
    
    if not sol.success:
        print(f"Warning: Integration did not succeed explicitly: {sol.message}")

    # 5. GENERATE PLOTS
    plot_t = sol.t / 3600 # Convert to hours
    
    # Extract states
    Y = sol.y.T
    
    # Proteins required: pEGFR, pMEK, pERK, DUSP, pCRAF, pAKT, p4EBP1
    # Indices from model:
    # m_pEGFR = Y[:, 2]; m_pCRAF = Y[:, 22]; m_pMEK = Y[:, 26]
    # m_pERK = Y[:, 28]; m_DUSP = Y[:, 30]; m_pAKT = Y[:, 52]; m_p4EBP1 = Y[:, 58]
    
    # Normalizations:
    # m_pEGFR_n = normit(m_pEGFR, 'max')
    # m_pMEK_n  = normit(m_pMEK, 'first')
    # m_pERK_n  = normit(m_pERK, 'first')
    # m_DUSP_n  = normit(m_DUSP, 'max')
    # m_pCRAF_n = normit(m_pCRAF, 'max')
    # m_pAKT_n  = normit(m_pAKT, 'max')
    # m_p4EBP1_n= normit(m_p4EBP1, 'max')

    proteins = [
        {'name': 'pEGFR',  'idx': 2,  'norm': 'max',   'color': 'orange'},
        {'name': 'pMEK',   'idx': 26, 'norm': 'first', 'color': 'green'},
        {'name': 'pERK',   'idx': 28, 'norm': 'first', 'color': 'blue'},
        {'name': 'DUSP',   'idx': 30, 'norm': 'max',   'color': 'purple'},
        {'name': 'pCRAF',  'idx': 22, 'norm': 'max',   'color': 'cyan'},
        {'name': 'pAKT',   'idx': 52, 'norm': 'max',   'color': 'red'},
        {'name': 'p4EBP1', 'idx': 58, 'norm': 'max',   'color': 'brown'}
    ]
    
    # Create a 4x2 grid (8 slots for 7 plots)
    fig, axes = plt.subplots(4, 2, figsize=(15, 20))
    axes = axes.flatten()
    
    for i, prot in enumerate(proteins):
        ax = axes[i]
        
        # Model Data
        model_vals = normit(Y[:, prot['idx']], prot['norm'])
        
        # Exp Data
        exp_vals = exp_data_norm[prot['name']]
        
        # Plotting
        ax.plot(plot_t, model_vals, '-', color=prot['color'], linewidth=2.5, label='Model Fit')
        ax.plot(time_hours, exp_vals, 'o', color='black', markersize=8, markeredgewidth=2, fillstyle='none', label='Experiment')
        
        ax.set_title(f"{prot['name']} Dynamics", fontsize=14, fontweight='bold')
        ax.set_xlabel('Time (Hours)', fontsize=12)
        ax.set_ylabel('Normalized Level', fontsize=12)
        ax.legend(loc='best')
        ax.grid(True, linestyle='--', alpha=0.5)
        
    # Remove the 8th empty subplot
    axes[7].axis('off')
    
    plt.tight_layout()
    plt.suptitle("Model Fit vs Experimental Data (All Proteins)",  y=1.02, fontsize=18)
    
    output_path = 'Full_Model_Fit_Analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_path}")
