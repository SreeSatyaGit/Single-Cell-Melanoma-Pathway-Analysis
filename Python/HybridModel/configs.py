import torch

def get_baseline_configs():
    """
    Returns the baseline parameter vector and initial conditions from MATLAB VemandTram.m
    """
    # 68 Parameters (Ordered as in params_vector)
    # Note: These are the guesses before optimization
    p = [
        10e-5, 10e-6, 10e-3,          # RTK on, off, cat
        10e-4, 7e-6, 20e-5,           # CRAF act, MEK phos, ERK phos
        10e-7, 20e-6, 10e-5, 7e-7, 36e-6, # RTK degrad, ERK inhib, Shc dephos, ptp, Grb2 bind
        10e-5, 10e-5, 17e-6,          # Sprty inhib, SOS bind, ERK phos SOS
        15.7342e-6, 1.7342e-4, 20.7342e-5, 1.7342e-8, # ERK phos CRAF, CRAF degrad, ERK phos MEK, MEK degrad
        2.7342e-4, 1.7342e-7, 10e-2, 1.7342e-04, 1.11e-8, # DUSP inhib ERK, ERK degrad, BRAF inhib, DUSP stop, DUSP max tx
        1.7e-06, 5.5000e-5, 10e-7, 8e-5, # SPRY max tx, SPRY come down, Degrad general, mRNA decay
        1.11e-8, 1.7e-06,             # DUSP max tx (duplicated in vector?), SPRY max tx (duplicated)
        1e-5, 1e-6,                   # ERK dephos, DUSP deg
        10e-5, 10e-5,                 # Her2 act, Her3 act
        10e-4, 10e-4, 10e-4, 10e-4,   # p85 bind EGFR, Her2, Her3, IGFR
        10e-5, 10e-4,                 # p85 unbind, PI3K recruit
        10e-4,                        # mTOR feedback
        10e-4, 10e-5,                 # PIP2 to PIP3, PTEN
        10e-5, 10e-7,                 # AKT act, AKT degrad
        10e-8, 10e-3, 10e-5, 10e-5,   # mTOR kb1, k43b1, 4EBP1 phos, 4EBP1 dephos
        5e-6, 5e-6,                   # KSR phos, dephos
        8e-6, 8e-6, 3e-6,             # BRAF route, CRAF route, KSR route
        1e-6, 1e-9, 1e-9, 2.0,        # Tram conc, Ki_RAF, Ki_KSR, n_tram
        1.0, 6e-6, 1e-5,              # Vem conc, dimer form, dimer dissoc
        0.5, 0.4, 1.5,                # Paradox gamma, Vem IC50, Vem Hill
        10e-5, 10e-4,                 # PDGFR act, p85 bind PDGFR
        10e-5, 10e-5,                 # S6K phos, S6K dephos
        0.05                          # K_displace
    ]

    # Initial Conditions (68 Species)
    IC_EGFR = [1.0, 0.35, 0.35]
    IC_Her2 = [1.0, 0.245, 0.245]
    IC_Her3 = [1.0, 0.203, 0.203]
    IC_SHC = [1.0, 0.0, 1.0]
    IC_Grb2_SOS = [0.0, 0.0]
    IC_HRAS = [1.0, 0.0] # Fixed from Matlab vector construction y0
    IC_NRAS = [1.0, 0.0]
    IC_KRAS = [1.0, 0.0, 0.0] # Matlab: [1.0, 0.0, 1.0] but wait, let's re-check y0 unpacking
    
    # Re-packing exactly as VemandTram.m Section 2
    y0_list = [
        1.0, 0.35, 0.35,              # EGFR [0:2]
        1.0, 0.245, 0.245,            # Her2 [3:5]
        1.0, 0.203, 0.203,            # Her3 [6:8]
        1.0, 0.0, 1.0,                # SHC [9:11]
        0.0, 0.0,                     # Grb2_SOS [12:13]
        0.0, 0.0,                     # HRAS [14:15]
        0.0, 0.0,                     # NRAS [16:17]
        1.0, 0.0, 1.0,                # KRAS [18:20] -> y[20] is active
        0.8, 0.366,                   # CRAF [21:22] -> y[22] is pCRAF
        1.0, 1.0,                     # BRAF [23:24] -> y[24] is BRAF^P
        1.0, 1.759,                   # MEK [25:26]
        1.0, 2.903,                   # ERK [27:28]
        1.0, 2.677,                   # DUSP [29:30]
        1.0, 1.0,                     # SPRY [31:32]
        1.0,                          # pERK_degrad [33]
        1.0,                          # pMEK_degrad [34]
        1.0,                          # pCRAF_degrad [35]
        1.0,                          # DUSP_stop [36]
        1.0, 0.0, 0.0,                # IGFR [37:39]
        1.0, 0.0,                     # IRS [40:41]
        1.0,                          # p85 [42]
        0.1, 0.1, 0.1, 0.1,           # p85_EGFR, Her2, Her3, IGFR [43:46]
        1.0, 0.2,                     # PI3K [47:48]
        1.0, 0.1,                     # PIP [49:50]
        1.0, 0.513,                   # AKT [51:52]
        0.0,                          # FOXO [53]
        1.0, 0.5,                     # mTORC [54:55]
        1.0, 0.5, 1.002,              # frebp1 [56:58]
        1.0, 0.0,                     # KSR [59:60]
        0.0,                          # BRAF_CRAF_dimer [61]
        1.0, 0.474, 0.474,            # PDGFR [62:64]
        1.0, 1.432,                   # S6K [65:66]
        0.1                           # p85_PDGFR [67]
    ]
    
    return torch.tensor(p, dtype=torch.float32), torch.tensor(y0_list, dtype=torch.float32)
