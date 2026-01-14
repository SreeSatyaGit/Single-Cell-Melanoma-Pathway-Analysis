import torch
from torch import nn

class MAPKPhysics(nn.Module):
    def __init__(self, params_vector, trainable_indices=None):
        super().__init__()
        if not isinstance(params_vector, torch.Tensor):
            params_vector = torch.tensor(params_vector, dtype=torch.float32)
        
        # Original baseline parameters (frozen)
        self.register_buffer('p_base', params_vector)
        
        # Selective refinement: only these indices become trainable
        self.trainable_indices = trainable_indices if trainable_indices is not None else []
        if self.trainable_indices:
            # Initialize deltas at 0.0 (log-space or additive? Let's do additive to start)
            # Actually, multi-scale parameters are better handled in log-space or as multipliers
            self.p_multipliers = nn.Parameter(torch.ones(len(self.trainable_indices)))
        else:
            self.p_multipliers = None

    @property
    def p(self):
        """Effective parameter vector"""
        if self.p_multipliers is None:
            return self.p_base
        
        # Start with baseline
        p_eff = self.p_base.clone()
        # Apply multipliers to selected indices
        # Clamp multipliers to [0.1, 10.0] to keep physics realistic
        multipliers = torch.clamp(self.p_multipliers, 0.1, 10.0)
        for i, idx in enumerate(self.trainable_indices):
            p_eff[idx] = self.p_base[idx] * multipliers[i]
        return p_eff

    def forward(self, t, y, drug_conc):
        """
        Calculates dy/dt for a batch of drug conditions.
        y: [Batch, 68]
        drug_conc: [Batch, 4]
        """
        # Clamp y to prevent NaNs
        y_pos = torch.clamp(y, 1e-12, 10.0)
        
        # Unpack species
        s = [y_pos[:, i:i+1] for i in range(68)]
        
        # Parameters (now immutable)
        ka1, kr1, kc1 = self.p[0], self.p[1], self.p[2]
        kpCraf, kpMek, kpErk = self.p[3], self.p[4], self.p[5]
        kDegradEgfr, kErkInbEgfr = self.p[6], self.p[7]
        kShcDephos, kptpDeg, kGrb2CombShc = self.p[8], self.p[9], self.p[10]
        kSprtyInbGrb2, kSosCombGrb2, kErkPhosSos = self.p[11], self.p[12], self.p[13]
        kErkPhosPcraf, kPcrafDegrad, kErkPhosMek, kMekDegrad = self.p[14], self.p[15], self.p[16], self.p[17]
        kDuspInbErk, kErkDeg, kinbBraf = self.p[18], self.p[19], self.p[20]
        kDuspStop, kDusps, kSproutyForm = self.p[21], self.p[22], self.p[23]
        kSprtyComeDown, kdegrad_ras_gap = self.p[24], self.p[25]
        km_Dusp, km_Sprty, kErkDephos, kDuspDeg = self.p[27], self.p[28], self.p[29], self.p[30]
        kHer2_act, kHer3_act = self.p[31], self.p[32]
        k_p85_bind_EGFR, k_p85_bind_Her2, k_p85_bind_Her3, k_p85_bind_IGFR = self.p[33], self.p[34], self.p[35], self.p[36]
        k_p85_unbind, k_PI3K_recruit, kMTOR_Feedback = self.p[37], self.p[38], self.p[39]
        k_PIP2_to_PIP3, k_PTEN, kAkt, kdegradAKT = self.p[40], self.p[41], self.p[42], self.p[43]
        kb1, k43b1, k4ebp1, k_4EBP1_dephos = self.p[44], self.p[45], self.p[46], self.p[47]
        kKSRphos, kKSRdephos = self.p[48], self.p[49]
        kMekByBraf, kMekByCraf, kMekByKSR = self.p[50], self.p[51], self.p[52]
        
        K_tram_RAF, K_tram_KSR, n_tram = self.p[54], self.p[55], self.p[56]
        kDimerForm, kDimerDissoc, kParadoxCRAF = self.p[58], self.p[59], self.p[60]
        IC50_vem, Hill_n_vem = self.p[61], self.p[62]
        kPDGFR_act, k_p85_bind_PDGFR, kS6K_phos, kS6K_dephos = self.p[63], self.p[64], self.p[65], self.p[66]
        K_displace = self.p[67]

        # Drugs
        Tram_ext, Vem_ext, RASi_ext, PI3Ki_ext = drug_conc[:, 0:1], drug_conc[:, 1:2], drug_conc[:, 2:3], drug_conc[:, 3:4]
        
        # ODE Logic (scaled to hours)
        d = [None] * 68
        scale = 3600.0
        
        # 1. EGFR
        d[0] = -ka1*s[0] + kr1*s[1]
        d[1] = ka1*s[0] - kr1*s[1] - kc1*s[1]
        d[2] = kc1*s[1] - kDegradEgfr*s[2] - kErkInbEgfr*s[28]*s[2]
        
        # 2. Her2
        d[3] = -kHer2_act*s[3] + kr1*s[4]
        d[4] = kHer2_act*s[3] - kr1*s[4] - kc1*s[4]
        d[5] = kc1*s[4] - kDegradEgfr*s[5] - kErkInbEgfr*s[28]*s[5]
        
        # 3. Her3
        d[6] = -kHer3_act*s[6] + kr1*s[7]
        d[7] = kHer3_act*s[6] - kr1*s[7] - kc1*s[7]
        d[8] = kc1*s[7] - kDegradEgfr*s[8] - kErkInbEgfr*s[28]*s[8]
        
        # 4. SOS/Grb2
        d[9] = -ka1*s[2]*s[9]
        d[10] = ka1*s[2]*s[9] - kShcDephos*s[11]*s[10]
        d[11] = -kptpDeg*s[10]*s[11] 
        d[12] = kGrb2CombShc*s[10]*s[2] - kSprtyInbGrb2*s[32]*s[12]
        d[13] = kSosCombGrb2*s[12]*s[10] - kErkPhosSos*s[28]*s[13]
        
        # 5. RAS
        panRAS_act = s[15] + s[17] + s[20] 
        d[14] = -ka1*s[13]*s[14] + kdegrad_ras_gap*s[15]
        d[15] = ka1*s[13]*s[14] - kdegrad_ras_gap*s[15]
        d[16] = -ka1*s[13]*s[16] + kdegrad_ras_gap*s[17]
        d[17] = ka1*s[13]*s[16] - kdegrad_ras_gap*s[17]
        d[18] = -ka1*s[13]*s[18] + kdegrad_ras_gap*s[20]
        d[19] = torch.zeros_like(s[0])
        d[20] = ka1*s[13]*s[18] - kdegrad_ras_gap*s[20] - (RASi_ext * 2.0 * s[20])
        
        # 6. RAF/Vem
        IC50_n = IC50_vem**Hill_n_vem
        Vem_n = Vem_ext**Hill_n_vem
        kBRAF_eff = ka1 * IC50_n / (IC50_n + Vem_n + 1e-9)
        dimerization = kDimerForm * s[23] * s[21] * Vem_ext
        d[21] = -kpCraf*panRAS_act*s[21] + kErkPhosPcraf*s[28]*s[22] + kPcrafDegrad*s[22]*s[35] - dimerization + kDimerDissoc*s[61]
        d[22] = kpCraf*panRAS_act*s[21] - kErkPhosPcraf*s[28]*s[22] - kPcrafDegrad*s[22]*s[35] + kParadoxCRAF*Vem_ext*s[61]
        d[23] = -kBRAF_eff*s[23]*panRAS_act - dimerization + kDimerDissoc*s[61]
        d[24] = kBRAF_eff*s[23]*panRAS_act - kinbBraf*s[24] - dimerization + kDimerDissoc*s[61]
        d[61] = dimerization - kDimerDissoc*s[61] - kPcrafDegrad*s[61]*s[35]
        
        # 7. MEK/Tram
        Load = panRAS_act + (s[21]+s[22]) + (s[23]+s[24])
        Ki = K_tram_RAF * (1.0 + (Load / (K_displace + 1e-9))**2)
        fMEK = 1.0 / (torch.clamp(1.0 + (Tram_ext / (Ki + 1e-9))**n_tram, 1.0, 100.0))
        raf_to_mek = (kpMek*s[22] + kMekByBraf*s[24] + kMekByCraf*s[22] + kpMek*s[61] + kMekByKSR*s[60])
        d[25] = -raf_to_mek*s[25] + kErkPhosMek*s[28]*s[26] + kMekDegrad*s[26]*s[34]
        d[26] = raf_to_mek*s[25] - kErkPhosMek*s[28]*s[26] - kMekDegrad*s[26]*s[34]
        
        # 8. ERK
        erk_act = kpErk * s[26] * s[27] * fMEK
        d[27] = -erk_act + kErkDephos*s[30]*s[28] + kErkDeg*s[28]*s[33]
        d[28] = erk_act - kErkDephos*s[30]*s[28] - kErkDeg*s[28]*s[33]
        
        # 9. Feedbacks
        d[29] = km_Dusp*s[28]/(1.0 + (km_Dusp/kDusps)*s[28] + 1e-9) - kDuspStop*s[29]*s[36] - kDuspDeg*s[29]*s[28]
        d[30] = -kDuspStop*s[29]*s[30]
        d[31] = km_Sprty*s[28]/(1.0 + (km_Sprty/kSproutyForm)*s[28] + 1e-9) - kSprtyComeDown*s[31]*s[32]
        d[32] = -kSprtyComeDown*s[31]*s[32]
        
        # Intermediates
        d[33] = -kErkDeg*s[28]*s[33]
        d[34] = -kMekDegrad*s[26]*s[34]
        d[35] = -kPcrafDegrad*s[22]*s[35]
        d[36] = -kDuspStop*s[29]*s[36]
        
        # 11. PI3K
        d[37] = -ka1*s[37] + kr1*s[38]
        d[38] = ka1*s[37] - kr1*s[38] - kc1*s[38]
        d[39] = kc1*s[38] - kErkInbEgfr*s[28]*s[39]
        d[40] = -ka1*s[2]*s[40]
        d[41] = ka1*s[2]*s[40]
        
        # 12. p85
        tp85 = s[43] + s[44] + s[45] + s[46] + s[67]
        d[43] = k_p85_bind_EGFR*s[2]*s[42] - k_p85_unbind*s[43]
        d[44] = k_p85_bind_Her2*s[5]*s[42] - k_p85_unbind*s[44]
        d[45] = k_p85_bind_Her3*s[8]*s[42] - k_p85_unbind*s[45]
        d[46] = k_p85_bind_IGFR*s[39]*s[42] - k_p85_unbind*s[46]
        d[67] = k_p85_bind_PDGFR*s[64]*s[42] - k_p85_unbind*s[67]
        d[42] = -(d[43] + d[44] + d[45] + d[46] + d[67])
        
        # 13. AKT
        pi3k_r = (k_PI3K_recruit * tp85 * s[47]) / (1.0 + PI3Ki_ext / 0.1)
        d[47] = -pi3k_r + kMTOR_Feedback*s[55]*s[48]
        d[48] = pi3k_r - kMTOR_Feedback*s[55]*s[48]
        d[49] = -k_PIP2_to_PIP3*s[48]*s[49] + k_PTEN*s[50]
        d[50] = k_PIP2_to_PIP3*s[48]*s[49] - k_PTEN*s[50]
        d[51] = -kAkt*s[50]*s[51] + kdegradAKT*s[52]
        d[52] = kAkt*s[50]*s[51] - kdegradAKT*s[52]
        
        # 14. mTOR
        d[53] = (torch.max(torch.tensor(0.0), 1.0 - s[52])) * kAkt / (1.0 + s[53]/(15e-5+1e-9))
        d[54] = -kAkt*s[52]*s[54] + kPcrafDegrad*s[55]
        d[55] = kAkt*s[52]*s[54] - kPcrafDegrad*s[55]
        d[56] = -k4ebp1*s[55]*s[56] + kb1*s[57] + k_4EBP1_dephos*s[58]
        d[57] = k4ebp1*s[55]*s[56] - kb1*s[57] - k43b1*s[57]
        d[58] = k43b1*s[57] - k_4EBP1_dephos*s[58]
        
        # 15. KSR
        d[59] = -kKSRphos*panRAS_act*s[59] + kKSRdephos*s[60]
        d[60] = kKSRphos*panRAS_act*s[59] - kKSRdephos*s[60]
        
        # 16. PDGFR
        d[62] = -kPDGFR_act*s[62] + kr1*s[63]
        d[63] = kPDGFR_act*s[62] - kr1*s[63] - kc1*s[63]
        d[64] = kc1*s[63] - kDegradEgfr*s[64] - kErkInbEgfr*s[28]*s[64]
        
        # 17. S6K
        d[65] = -kS6K_phos*s[55]*s[65] + kS6K_dephos*s[66]
        d[66] = kS6K_phos*s[55]*s[65] - kS6K_dephos*s[66]
        
        return torch.cat(d, dim=-1) * scale

    def get_sensors(self, Y):
        sensors = {}
        def norm(v): return (v - v.min(dim=0, keepdim=True)[0]) / (v.max(dim=0, keepdim=True)[0] - v.min(dim=0, keepdim=True)[0] + 1e-6)
        sensors['pEGFR'] = norm(Y[:, :, 2])
        sensors['panRAS'] = norm(Y[:, :, 15] + Y[:, :, 17] + Y[:, :, 20])
        sensors['pCRAF'] = norm(Y[:, :, 22])
        sensors['pMEK'] = norm(Y[:, :, 26])
        sensors['pERK'] = norm(Y[:, :, 28])
        sensors['DUSP'] = norm(Y[:, :, 30])
        sensors['pAKT'] = norm(Y[:, :, 52])
        sensors['p4ebp1'] = norm(Y[:, :, 58])
        sensors['her2'] = norm(Y[:, :, 5])
        sensors['her3'] = norm(Y[:, :, 8])
        sensors['pDGFR'] = norm(Y[:, :, 64])
        sensors['pS6k'] = norm(Y[:, :, 66])
        return sensors
