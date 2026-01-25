#!/usr/bin/env python3
"""
MAPK/PI3K Signaling Pathway - Vemurafenib Treatment Fitting
============================================================

This script fits the MAPK/PI3K model to experimental Western blot data
from Vemurafenib treatment. Both experimental and model outputs are
normalized to 0-1 scale for comparison.

Experimental Data:
- 6 time points (assumed: 0, 15, 30, 60, 120, 240 min)
- 12 readouts: panRAS, pMEK, pERK, DUSP, pEGFR, pCRAF, pAKT, p4EBP1, pS6K, Her2, Her3, pDGFR
- Western blot data normalized to vinculin
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn

from torchdiffeq import odeint

torch.set_default_dtype(torch.float64)

# ============================================================================
# EXPERIMENTAL DATA (Vemurafenib treatment, Western blot normalized to vinculin)
# ============================================================================
# Time points: 0, 1, 4, 8, 24, 48 hours (converted to minutes for model consistency)
EXPERIMENTAL_TIME_POINTS_HOURS = torch.tensor([0.0, 1.0, 4.0, 8.0, 24.0, 48.0])
EXPERIMENTAL_TIME_POINTS = EXPERIMENTAL_TIME_POINTS_HOURS * 60.0  # Convert to minutes

# Raw experimental data
EXPERIMENTAL_DATA_RAW = {
    'panRAS': torch.tensor([0.839133594, 0.833259289, 0.919508516, 1.240888235, 1.582734859, 1.468310571]),
    'pMEK':   torch.tensor([1.75938884, 0.170160085, 0.095112609, 0.201000276, 0.219207054, 0.502831668]),
    'pERK':   torch.tensor([2.903209735, 0.207867788, 0.303586121, 0.805254439, 1.408362153, 1.847606441]),
    'DUSP':   torch.tensor([2.677161325, 2.782754577, 1.130758062, 0.395642757, 0.828575853, 0.916618219]),
    'pEGFR':  torch.tensor([0.291928893, 0.392400458, 0.265016688, 0.394238749, 0.006158316, 0.008115099]),
    'pCRAF':  torch.tensor([0.366397596, 0.537106733, 0.465541704, 0.586732657, 1.102322681, 0.269181259]),
    'pAKT':   torch.tensor([0.513544148, 0.613178403, 1.03451863, 1.113391047, 0.535242724, 0.538273551]),
    'p4EBP1': torch.tensor([1.002468056, 1.276793699, 1.252681407, 1.707504483, 1.271216967, 0.61389625]),
    'pS6K':   torch.tensor([1.432459522, 1.520433646, 1.542177411, 1.248505245, 0.109963216, 0.013374136]),
    'Her2':   torch.tensor([0.245236744, 0.177917339, 0.239075259, 0.306884773, 1.066654783, 1.005085151]),
    'Her3':   torch.tensor([0.203233765, 0.194358998, 0.303475212, 0.674083831, 0.89702403, 0.459831389]),
    'pPDGFR': torch.tensor([0.474174188, 0.492132953, 0.743620725, 1.266460499, 2.514722273, 2.482761079]),
}

# Mapping from experimental readout names to model state indices
# Based on the 68-species state vector in the MAPK/PI3K model
READOUT_TO_STATE_INDEX = {
    'panRAS': [15, 17, 20],  # HRAS-GTP + NRAS-GTP + KRAS-GTP (sum)
    'pMEK':   [26],          # pMEK
    'pERK':   [28],          # pERK
    'DUSP':   [30],          # DUSP
    'pEGFR':  [2],           # pEGFR
    'pCRAF':  [22],          # pCRAF
    'pAKT':   [52],          # pAKT
    'p4EBP1': [58],          # p4EBP1
    'pS6K':   [66],          # pS6K
    'Her2':   [5],           # pHer2
    'Her3':   [8],           # pHer3
    'pPDGFR': [64],          # pPDGFR
}


def normalize_to_01(data: torch.Tensor) -> tuple:
    """
    Normalize data to 0-1 range using min-max normalization.
    
    Args:
        data: Input tensor
        
    Returns:
        normalized_data: Data scaled to [0, 1]
        min_val: Minimum value (for inverse transform)
        max_val: Maximum value (for inverse transform)
    """
    min_val = data.min()
    max_val = data.max()
    range_val = max_val - min_val
    
    # Handle constant data
    if range_val == 0:
        return torch.zeros_like(data), min_val, max_val
    
    normalized = (data - min_val) / range_val
    return normalized, min_val, max_val


def get_normalized_experimental_data():
    """
    Normalize all experimental data to 0-1 scale.
    
    Returns:
        dict: Normalized experimental data for each readout
        dict: Normalization parameters (min, max) for each readout
    """
    normalized_data = {}
    norm_params = {}
    
    for name, raw_data in EXPERIMENTAL_DATA_RAW.items():
        normalized, min_val, max_val = normalize_to_01(raw_data)
        normalized_data[name] = normalized
        norm_params[name] = {'min': min_val, 'max': max_val}
    
    return normalized_data, norm_params


def extract_model_readouts(y: torch.Tensor, t: torch.Tensor, 
                           exp_times: torch.Tensor) -> dict:
    """
    Extract model outputs at experimental time points.
    
    Args:
        y: Full model solution (n_points, 68)
        t: Time points from simulation
        exp_times: Experimental time points
        
    Returns:
        dict: Model outputs for each readout at experimental time points
    """
    # Interpolate to experimental time points
    readouts = {}
    
    for name, indices in READOUT_TO_STATE_INDEX.items():
        values = []
        for exp_t in exp_times:
            # Find closest time point in simulation
            idx = torch.argmin(torch.abs(t - exp_t))
            
            # Sum if multiple indices (e.g., panRAS)
            val = sum(y[idx, i].item() for i in indices)
            values.append(val)
        
        readouts[name] = torch.tensor(values)
    
    return readouts


def normalize_model_readouts(model_readouts: dict, 
                             use_experimental_scale: bool = False,
                             exp_norm_params: dict = None) -> dict:
    """
    Normalize model readouts to 0-1 scale.
    
    Args:
        model_readouts: Raw model outputs
        use_experimental_scale: If True, use experimental min/max for scaling
        exp_norm_params: Normalization parameters from experimental data
        
    Returns:
        dict: Normalized model outputs
    """
    normalized = {}
    
    for name, data in model_readouts.items():
        if use_experimental_scale and exp_norm_params and name in exp_norm_params:
            # Use experimental scale
            min_val = exp_norm_params[name]['min']
            max_val = exp_norm_params[name]['max']
            range_val = max_val - min_val
            
            if range_val == 0:
                normalized[name] = torch.zeros_like(data)
            else:
                normalized[name] = (data - min_val) / range_val
        else:
            # Self-normalize
            normalized[name], _, _ = normalize_to_01(data)
    
    return normalized


# ============================================================================
# MODEL DEFINITION (Simplified for faster simulation)
# ============================================================================
class MAPKParameters(nn.Module):
    """Container for all kinetic parameters of the MAPK/PI3K pathway."""
    
    def __init__(self):
        super().__init__()
        
        # RTK parameters
        self.ka1 = nn.Parameter(torch.tensor([0.1]))
        self.kr1 = nn.Parameter(torch.tensor([0.05]))
        self.kc1 = nn.Parameter(torch.tensor([0.2]))
        
        # RAF/MEK/ERK cascade
        self.kpCraf = nn.Parameter(torch.tensor([0.1]))
        self.kpMek = nn.Parameter(torch.tensor([0.1]))
        self.kpErk = nn.Parameter(torch.tensor([0.1]))
        
        # Degradation and feedback
        self.kDegradEgfr = nn.Parameter(torch.tensor([0.01]))
        self.kErkInbEgfr = nn.Parameter(torch.tensor([0.01]))
        self.kShcDephos = nn.Parameter(torch.tensor([0.05]))
        self.kptpDeg = nn.Parameter(torch.tensor([0.01]))
        self.kGrb2CombShc = nn.Parameter(torch.tensor([0.1]))
        self.kSprtyInbGrb2 = nn.Parameter(torch.tensor([0.01]))
        self.kSosCombGrb2 = nn.Parameter(torch.tensor([0.1]))
        self.kErkPhosSos = nn.Parameter(torch.tensor([0.01]))
        self.kErkPhosPcraf = nn.Parameter(torch.tensor([0.01]))
        self.kPcrafDegrad = nn.Parameter(torch.tensor([0.01]))
        self.kErkPhosMek = nn.Parameter(torch.tensor([0.01]))
        self.kMekDegrad = nn.Parameter(torch.tensor([0.01]))
        self.kDuspInbErk = nn.Parameter(torch.tensor([0.05]))
        self.kErkDeg = nn.Parameter(torch.tensor([0.01]))
        self.kinbBraf = nn.Parameter(torch.tensor([0.01]))
        self.kDuspStop = nn.Parameter(torch.tensor([0.01]))
        self.kDusps = nn.Parameter(torch.tensor([0.1]))
        self.kSproutyForm = nn.Parameter(torch.tensor([0.1]))
        self.kSprtyComeDown = nn.Parameter(torch.tensor([0.01]))
        self.kdegrad = nn.Parameter(torch.tensor([0.01]))
        self.km_Sprty_decay = nn.Parameter(torch.tensor([0.01]))
        self.km_Dusp = nn.Parameter(torch.tensor([0.1]))
        self.km_Sprty = nn.Parameter(torch.tensor([0.1]))
        self.kErkDephos = nn.Parameter(torch.tensor([0.05]))
        self.kDuspDeg = nn.Parameter(torch.tensor([0.01]))
        
        # Her2/Her3 activation
        self.kHer2_act = nn.Parameter(torch.tensor([0.05]))
        self.kHer3_act = nn.Parameter(torch.tensor([0.05]))
        
        # p85 binding parameters
        self.k_p85_bind_EGFR = nn.Parameter(torch.tensor([0.1]))
        self.k_p85_bind_Her2 = nn.Parameter(torch.tensor([0.1]))
        self.k_p85_bind_Her3 = nn.Parameter(torch.tensor([0.1]))
        self.k_p85_bind_IGFR = nn.Parameter(torch.tensor([0.1]))
        self.k_p85_unbind = nn.Parameter(torch.tensor([0.05]))
        self.k_PI3K_recruit = nn.Parameter(torch.tensor([0.1]))
        
        # mTOR feedback
        self.kMTOR_Feedback = nn.Parameter(torch.tensor([0.01]))
        
        # PIP conversion
        self.k_PIP2_to_PIP3 = nn.Parameter(torch.tensor([0.1]))
        self.k_PTEN = nn.Parameter(torch.tensor([0.05]))
        
        # AKT
        self.kAkt = nn.Parameter(torch.tensor([0.1]))
        self.kdegradAKT = nn.Parameter(torch.tensor([0.01]))
        
        # 4EBP1
        self.kb1 = nn.Parameter(torch.tensor([0.05]))
        self.k43b1 = nn.Parameter(torch.tensor([0.05]))
        self.k4ebp1 = nn.Parameter(torch.tensor([0.1]))
        self.k_4EBP1_dephos = nn.Parameter(torch.tensor([0.05]))
        
        # KSR
        self.kKSRphos = nn.Parameter(torch.tensor([0.1]))
        self.kKSRdephos = nn.Parameter(torch.tensor([0.05]))
        
        # MEK phosphorylation
        self.kMekByBraf = nn.Parameter(torch.tensor([0.1]))
        self.kMekByCraf = nn.Parameter(torch.tensor([0.1]))
        self.kMekByKSR = nn.Parameter(torch.tensor([0.05]))
        
        # Trametinib parameters
        self.Tram = nn.Parameter(torch.tensor([0.0]))
        self.K_tram_RAF = nn.Parameter(torch.tensor([1.0]))
        self.K_tram_KSR = nn.Parameter(torch.tensor([1.0]))
        self.n_tram = nn.Parameter(torch.tensor([1.0]))
        
        # Vemurafenib parameters
        self.Vemurafenib = nn.Parameter(torch.tensor([1.0]))  # Default: drug applied
        self.kDimerForm = nn.Parameter(torch.tensor([0.1]))
        self.kDimerDissoc = nn.Parameter(torch.tensor([0.05]))
        self.kParadoxCRAF = nn.Parameter(torch.tensor([0.1]))
        self.IC50_vem = nn.Parameter(torch.tensor([1.0]))
        self.Hill_n_vem = nn.Parameter(torch.tensor([1.0]))
        
        # PDGFR and S6K
        self.kPDGFR_act = nn.Parameter(torch.tensor([0.05]))
        self.k_p85_bind_PDGFR = nn.Parameter(torch.tensor([0.1]))
        self.kS6K_phos = nn.Parameter(torch.tensor([0.1]))
        self.kS6K_dephos = nn.Parameter(torch.tensor([0.05]))
        
        # Crosstalk parameters
        self.kRAS_PI3K = nn.Parameter(torch.tensor([0.05]))
        self.kERK_IRS_inhibit = nn.Parameter(torch.tensor([0.02]))
        self.kERK_PTEN_activate = nn.Parameter(torch.tensor([0.01]))
        self.kAKT_CRAF_inhibit = nn.Parameter(torch.tensor([0.02]))
        self.kS6K_IRS_inhibit = nn.Parameter(torch.tensor([0.02]))
        self.kERK_GAB1_inhibit = nn.Parameter(torch.tensor([0.01]))
        self.kAKT_TSC2_phos = nn.Parameter(torch.tensor([0.05]))
        self.kERK_RSK_activate = nn.Parameter(torch.tensor([0.03]))


class MAPKPI3KODE(nn.Module):
    """ODE system for the MAPK/PI3K signaling pathway."""
    
    def __init__(self, params=None):
        super().__init__()
        self.params = params if params is not None else MAPKParameters()
        self.eps = 1e-10
    
    def forward(self, t, y):
        p = self.params
        y = torch.clamp(y, min=0)
        dydt = torch.zeros_like(y)
        
        # EGFR [0-2]
        dydt[0] = -p.ka1 * y[0] + p.kr1 * y[1]
        dydt[1] = p.ka1 * y[0] - p.kr1 * y[1] - p.kc1 * y[1]
        dydt[2] = p.kc1 * y[1] - p.kDegradEgfr * y[2] - p.kErkInbEgfr * y[28] * y[2]
        
        # Her2 [3-5]
        dydt[3] = -p.kHer2_act * y[3] + p.kr1 * y[4]
        dydt[4] = p.kHer2_act * y[3] - p.kr1 * y[4] - p.kc1 * y[4]
        dydt[5] = p.kc1 * y[4] - p.kDegradEgfr * y[5] - p.kErkInbEgfr * y[28] * y[5]
        
        # Her3 [6-8]
        dydt[6] = -p.kHer3_act * y[6] + p.kr1 * y[7]
        dydt[7] = p.kHer3_act * y[6] - p.kr1 * y[7] - p.kc1 * y[7]
        dydt[8] = p.kc1 * y[7] - p.kDegradEgfr * y[8] - p.kErkInbEgfr * y[28] * y[8]
        
        # Shc/Grb2/SOS [9-13]
        dydt[9] = -p.ka1 * y[2] * y[9]
        dydt[10] = p.ka1 * y[2] * y[9] - p.kShcDephos * y[11] * y[10]
        dydt[11] = -p.kptpDeg * y[10] * y[11]
        dydt[12] = p.kGrb2CombShc * y[10] * y[2] - p.kSprtyInbGrb2 * y[26] * y[12]
        dydt[13] = p.kSosCombGrb2 * y[12] * y[10] - p.kErkPhosSos * y[24] * y[13]
        
        # RAS [14-20]
        dydt[14] = -p.ka1 * y[13] * y[14]
        dydt[15] = p.ka1 * y[13] * y[14]
        dydt[16] = -p.ka1 * y[13] * y[16]
        dydt[17] = p.ka1 * y[13] * y[16]
        dydt[18] = -p.ka1 * y[13] * y[18]
        dydt[19] = p.ka1 * y[13] * y[18] - p.ka1 * y[19] * y[20]
        dydt[20] = -p.ka1 * y[19] * y[20]
        
        # RAF with Vemurafenib [21-24, 61]
        IC50_n = p.IC50_vem ** p.Hill_n_vem
        Vem_n = p.Vemurafenib ** p.Hill_n_vem
        kBRAF_eff = p.ka1 * IC50_n / (IC50_n + Vem_n + self.eps)
        paradox_activation = p.kParadoxCRAF * p.Vemurafenib * y[61]
        AKT_inhibition_CRAF = p.kAKT_CRAF_inhibit * y[52] * y[21]
        
        dydt[21] = (-p.kpCraf * y[19] * y[21] 
                    + p.kErkPhosPcraf * y[28] * y[22] 
                    + p.kPcrafDegrad * y[22] * y[35]
                    - p.kDimerForm * y[24] * y[21] * p.Vemurafenib 
                    + p.kDimerDissoc * y[61]
                    - AKT_inhibition_CRAF)
        dydt[22] = (p.kpCraf * y[19] * y[21] 
                    - p.kErkPhosPcraf * y[28] * y[22] 
                    - p.kPcrafDegrad * y[22] * y[35]
                    + paradox_activation)
        dydt[23] = (-kBRAF_eff * y[23] * y[19] 
                    - p.kDimerForm * y[24] * y[21] * p.Vemurafenib 
                    + p.kDimerDissoc * y[61])
        dydt[24] = (kBRAF_eff * y[23] * y[19] 
                    - p.kinbBraf * y[24] 
                    - p.kDimerForm * y[24] * y[21] * p.Vemurafenib 
                    + p.kDimerDissoc * y[61])
        dydt[61] = (p.kDimerForm * y[24] * y[21] * p.Vemurafenib 
                    - p.kDimerDissoc * y[61] 
                    - p.kPcrafDegrad * y[61] * y[35])
        
        # MEK [25-26]
        raf_to_mek = p.kpMek * y[22] + p.kMekByBraf * y[24] + p.kMekByCraf * y[22]
        ksr_to_mek = p.kMekByKSR * y[60]
        dydt[25] = (-(raf_to_mek + ksr_to_mek) * y[25] 
                    + p.kErkPhosMek * y[28] * y[26] 
                    + p.kMekDegrad * y[26] * y[34])
        dydt[26] = ((raf_to_mek + ksr_to_mek) * y[25] 
                    - p.kErkPhosMek * y[28] * y[26] 
                    - p.kMekDegrad * y[26] * y[34])
        
        # ERK [27-28]
        dydt[27] = (-p.kpErk * y[26] * y[27] 
                    + p.kDuspInbErk * y[30] * y[28] 
                    + p.kErkDeg * y[28] * y[33] 
                    + p.kErkDephos * y[30] * y[28])
        dydt[28] = (p.kpErk * y[26] * y[27] 
                    - p.kDuspInbErk * y[30] * y[28] 
                    - p.kErkDeg * y[28] * y[33] 
                    - p.kErkDephos * y[30] * y[28])
        
        # DUSP [29-30]
        denom_dusp = 1 + (p.km_Dusp / (p.kDusps + self.eps)) * y[28]
        dydt[29] = (p.km_Dusp * y[28] / (denom_dusp + self.eps) 
                    - p.kDuspStop * y[29] * y[36] 
                    - p.kDuspDeg * y[29] * y[28])
        dydt[30] = -p.kDuspStop * y[29] * y[30]
        
        # SPRY [31-32]
        denom_spry = 1 + (p.km_Sprty / (p.kSproutyForm + self.eps)) * y[28]
        dydt[31] = (p.km_Sprty * y[28] / (denom_spry + self.eps) 
                    - p.kSprtyComeDown * y[31] * y[32])
        dydt[32] = -p.kSprtyComeDown * y[31] * y[32]
        
        # Degradation trackers [33-36]
        dydt[33] = -p.kErkDeg * y[28] * y[33]
        dydt[34] = -p.kMekDegrad * y[26] * y[34]
        dydt[35] = -p.kPcrafDegrad * y[22] * y[35]
        dydt[36] = -p.kDuspStop * y[29] * y[36]
        
        # IGFR [37-39]
        dydt[37] = -p.ka1 * y[37] + p.kr1 * y[38]
        dydt[38] = p.ka1 * y[37] - p.kr1 * y[38] - p.kc1 * y[38]
        dydt[39] = p.kc1 * y[38] - p.kErkInbEgfr * y[28] * y[39]
        
        # IRS [40-41]
        ERK_inhibit_IRS = p.kERK_IRS_inhibit * y[28] * y[41]
        S6K_inhibit_IRS = p.kS6K_IRS_inhibit * y[66] * y[41]
        dydt[40] = (-p.ka1 * y[2] * y[40] + ERK_inhibit_IRS + S6K_inhibit_IRS)
        dydt[41] = (p.ka1 * y[2] * y[40] - ERK_inhibit_IRS - S6K_inhibit_IRS)
        
        # p85 [42-46, 67]
        ERK_GAB1_inhibition_factor = 1.0 / (1.0 + p.kERK_GAB1_inhibit * y[28])
        dydt[42] = (-p.k_p85_bind_EGFR * y[2] * y[42] * ERK_GAB1_inhibition_factor
                    - p.k_p85_bind_Her2 * y[5] * y[42] * ERK_GAB1_inhibition_factor
                    - p.k_p85_bind_Her3 * y[8] * y[42] * ERK_GAB1_inhibition_factor
                    - p.k_p85_bind_IGFR * y[39] * y[42]
                    - p.k_p85_bind_PDGFR * y[64] * y[42]
                    + p.k_p85_unbind * (y[43] + y[44] + y[45] + y[46] + y[67]))
        dydt[43] = p.k_p85_bind_EGFR * y[2] * y[42] * ERK_GAB1_inhibition_factor - p.k_p85_unbind * y[43]
        dydt[44] = p.k_p85_bind_Her2 * y[5] * y[42] * ERK_GAB1_inhibition_factor - p.k_p85_unbind * y[44]
        dydt[45] = p.k_p85_bind_Her3 * y[8] * y[42] * ERK_GAB1_inhibition_factor - p.k_p85_unbind * y[45]
        dydt[46] = p.k_p85_bind_IGFR * y[39] * y[42] - p.k_p85_unbind * y[46]
        dydt[67] = p.k_p85_bind_PDGFR * y[64] * y[42] - p.k_p85_unbind * y[67]
        
        # PI3K [47-48]
        total_RAS_GTP = y[15] + y[17] + y[20]
        RAS_PI3K_activation = p.kRAS_PI3K * total_RAS_GTP * y[47]
        total_p85_RTK = y[43] + y[44] + y[45] + y[46] + y[67]
        dydt[47] = (-p.k_PI3K_recruit * total_p85_RTK * y[47]
                    - RAS_PI3K_activation
                    + p.kMTOR_Feedback * y[55] * y[48])
        dydt[48] = (p.k_PI3K_recruit * total_p85_RTK * y[47]
                    + RAS_PI3K_activation
                    - p.kMTOR_Feedback * y[55] * y[48])
        
        # PIP [49-50]
        PTEN_activity = p.k_PTEN * (1.0 + p.kERK_PTEN_activate * y[28])
        dydt[49] = -p.k_PIP2_to_PIP3 * y[48] * y[49] + PTEN_activity * y[50]
        dydt[50] = p.k_PIP2_to_PIP3 * y[48] * y[49] - PTEN_activity * y[50]
        
        # AKT [51-52]
        dydt[51] = -p.kAkt * y[50] * y[51] + p.kdegradAKT * y[52]
        dydt[52] = p.kAkt * y[50] * y[51] - p.kdegradAKT * y[52]
        
        # FOXO [53]
        denom_foxo = 1 + (y[53] / 15e-5)
        dydt[53] = (1 - y[52]) * p.ka1 / (denom_foxo + self.eps)
        
        # mTORC [54-55]
        AKT_mTORC_activation = p.ka1 * y[52] * y[54]
        RSK_mTORC_activation = p.kERK_RSK_activate * y[28] * y[54]
        TSC2_inactivation = p.kAKT_TSC2_phos * y[52]
        dydt[54] = (-AKT_mTORC_activation - RSK_mTORC_activation + p.kdegrad * y[55])
        dydt[55] = (AKT_mTORC_activation + RSK_mTORC_activation 
                    + TSC2_inactivation * y[55] - p.kdegrad * y[55])
        
        # 4EBP1 [56-58]
        dydt[56] = (-p.k4ebp1 * y[55] * y[56] + p.kb1 * y[57] + p.k_4EBP1_dephos * y[58])
        dydt[57] = (p.k4ebp1 * y[55] * y[56] - p.kb1 * y[57] - p.k43b1 * y[57])
        dydt[58] = p.k43b1 * y[57] - p.k_4EBP1_dephos * y[58]
        
        # KSR [59-60]
        dydt[59] = -p.kKSRphos * total_RAS_GTP * y[59] + p.kKSRdephos * y[60]
        dydt[60] = p.kKSRphos * total_RAS_GTP * y[59] - p.kKSRdephos * y[60]
        
        # PDGFR [62-64]
        dydt[62] = -p.kPDGFR_act * y[62] + p.kr1 * y[63]
        dydt[63] = p.kPDGFR_act * y[62] - p.kr1 * y[63] - p.kc1 * y[63]
        dydt[64] = p.kc1 * y[63] - p.kDegradEgfr * y[64] - p.kErkInbEgfr * y[28] * y[64]
        
        # S6K [65-66]
        dydt[65] = -p.kS6K_phos * y[55] * y[65] + p.kS6K_dephos * y[66]
        dydt[66] = p.kS6K_phos * y[55] * y[65] - p.kS6K_dephos * y[66]
        
        return dydt


def get_default_initial_conditions():
    """Get default initial conditions for all 68 species."""
    y0 = torch.zeros(68)
    
    # RTKs
    y0[0] = 1.0      # EGFR
    y0[3] = 0.5      # Her2
    y0[6] = 0.5      # Her3
    y0[37] = 0.5     # IGFR
    y0[62] = 0.3     # PDGFR
    
    # Signaling intermediates
    y0[9] = 1.0      # Shc
    y0[12] = 1.0     # Grb2
    
    # RAS
    y0[14] = 1.0     # HRAS-GDP
    y0[16] = 1.0     # NRAS-GDP
    y0[18] = 1.0     # KRAS-GDP
    y0[20] = 0.1     # KRAS intermediate
    
    # RAF
    y0[21] = 1.0     # CRAF
    y0[23] = 1.0     # BRAF
    
    # MEK/ERK
    y0[25] = 1.0     # MEK
    y0[27] = 1.0     # ERK
    
    # Feedback
    y0[30] = 0.1     # DUSP
    y0[32] = 0.1     # SPRY
    
    # Degradation trackers
    y0[33] = 1.0
    y0[34] = 1.0
    y0[35] = 1.0
    y0[36] = 1.0
    
    # PI3K pathway
    y0[40] = 0.5     # IRS
    y0[42] = 1.0     # p85
    y0[47] = 1.0     # PI3K
    y0[49] = 1.0     # PIP2
    y0[51] = 1.0     # AKT
    y0[53] = 0.1     # FOXO
    y0[54] = 1.0     # mTORC
    y0[56] = 1.0     # 4EBP1
    y0[59] = 0.5     # KSR
    y0[65] = 1.0     # S6K
    
    return y0


def simulate_pathway(model, y0, t_span, n_points=100, method='dopri5'):
    """Simulate the MAPK/PI3K pathway."""
    t = torch.linspace(t_span[0], t_span[1], n_points)
    
    with torch.no_grad():
        y = odeint(model, y0, t, method=method)
    
    return t, y


def plot_comparison(exp_data_norm, model_data_norm, exp_times, save_dir=None):
    """
    Plot comparison between experimental and model data.
    Both should already be normalized to 0-1 scale.
    """
    n_readouts = len(exp_data_norm)
    n_cols = 4
    n_rows = (n_readouts + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten()
    
    exp_times_np = exp_times.numpy()
    
    for i, (name, exp_vals) in enumerate(exp_data_norm.items()):
        ax = axes[i]
        
        exp_vals_np = exp_vals.numpy()
        model_vals_np = model_data_norm[name].numpy()
        
        # Plot experimental data
        ax.scatter(exp_times_np, exp_vals_np, s=80, c='red', marker='o', 
                   label='Experimental', zorder=5, edgecolors='darkred', linewidths=1)
        
        # Plot model prediction
        ax.plot(exp_times_np, model_vals_np, 'b-', linewidth=2, label='Model')
        ax.scatter(exp_times_np, model_vals_np, s=40, c='blue', marker='s', 
                   zorder=4, edgecolors='darkblue', linewidths=1)
        
        ax.set_xlabel('Time (min)', fontsize=11)
        ax.set_ylabel('Normalized (0-1)', fontsize=11)
        ax.set_title(name, fontsize=13, fontweight='bold')
        ax.legend(fontsize=9)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Hide empty subplots
    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)
    
    plt.suptitle('Vemurafenib Treatment: Model vs Experimental Data\n(Both normalized to 0-1 scale)', 
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, 'vemurafenib_comparison.png')
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, axes


def compute_loss(exp_data_norm, model_data_norm):
    """Compute MSE loss between experimental and model data."""
    total_loss = 0.0
    losses = {}
    
    for name in exp_data_norm:
        mse = torch.mean((exp_data_norm[name] - model_data_norm[name])**2).item()
        losses[name] = mse
        total_loss += mse
    
    return total_loss / len(exp_data_norm), losses


def main():
    parser = argparse.ArgumentParser(description='MAPK/PI3K Model - Vemurafenib Fitting')
    parser.add_argument('--t_end', type=float, default=2880.0,  # 48 hours in minutes
                        help='Simulation end time (min)')
    parser.add_argument('--n_points', type=int, default=100, 
                        help='Number of time points')
    parser.add_argument('--method', type=str, default='dopri5',
                        choices=['dopri5', 'rk4', 'euler', 'midpoint'],
                        help='ODE solver method')
    parser.add_argument('--vemurafenib', type=float, default=1.0,
                        help='Vemurafenib concentration (nM)')
    parser.add_argument('--save_dir', type=str, default=None,
                        help='Directory to save figures')
    args = parser.parse_args()
    
    print("=" * 70)
    print("MAPK/PI3K Model - Vemurafenib Treatment Fitting")
    print("=" * 70)
    print(f"\nSimulation: 0 to {args.t_end:.1f} min ({args.n_points} points)")
    print(f"Solver: {args.method}")
    print(f"Vemurafenib: {args.vemurafenib} nM")
    print("=" * 70)
    
    # Get normalized experimental data
    print("\nNormalizing experimental data to 0-1 scale...")
    exp_data_norm, exp_norm_params = get_normalized_experimental_data()
    
    print("\nExperimental data (normalized):")
    for name, data in exp_data_norm.items():
        print(f"  {name}: min={data.min().item():.3f}, max={data.max().item():.3f}")
    
    # Create model
    params = MAPKParameters()
    params.Vemurafenib = nn.Parameter(torch.tensor([args.vemurafenib]))
    model = MAPKPI3KODE(params=params)
    
    # Get initial conditions
    y0 = get_default_initial_conditions()
    
    # Simulate
    print("\nRunning simulation...")
    t, y = simulate_pathway(
        model, y0,
        t_span=(0.0, args.t_end),
        n_points=args.n_points,
        method=args.method
    )
    print("Simulation complete!")
    
    # Extract model readouts at experimental time points
    print("\nExtracting model outputs at experimental time points...")
    model_readouts_raw = extract_model_readouts(y, t, EXPERIMENTAL_TIME_POINTS)
    
    # Normalize model outputs (self-normalize each readout to 0-1)
    print("Normalizing model outputs to 0-1 scale...")
    model_readouts_norm = normalize_model_readouts(model_readouts_raw)
    
    print("\nModel predictions (normalized):")
    for name, data in model_readouts_norm.items():
        print(f"  {name}: min={data.min().item():.3f}, max={data.max().item():.3f}")
    
    # Compute loss
    avg_loss, losses = compute_loss(exp_data_norm, model_readouts_norm)
    print(f"\nAverage MSE Loss: {avg_loss:.4f}")
    print("\nPer-readout MSE losses:")
    for name, loss in sorted(losses.items(), key=lambda x: x[1], reverse=True):
        print(f"  {name}: {loss:.4f}")
    
    # Plot comparison
    print("\nGenerating comparison plots...")
    plot_comparison(exp_data_norm, model_readouts_norm, 
                   EXPERIMENTAL_TIME_POINTS, save_dir=args.save_dir)
    
    plt.show()
    
    print("\nDone!")


if __name__ == "__main__":
    main()
