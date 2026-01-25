#!/usr/bin/env python3
"""
EGF-EGFR Signaling Pathway ODE Model
=====================================

This example models the EGFR (Epidermal Growth Factor Receptor) signaling pathway:

    EGF + EGFR <-> EGF*EGFR -> pEGFR

Reactions:
    1. EGF + EGFR -> EGF*EGFR  (binding, rate k1)
    2. EGF*EGFR -> EGF + EGFR  (unbinding, rate k2)  
    3. EGF*EGFR -> pEGFR       (phosphorylation, rate k3)
    4. pEGFR -> EGFR           (dephosphorylation/recycling, rate k4)

Species:
    - EGF: Epidermal Growth Factor (ligand)
    - EGFR: EGF Receptor (unbound)
    - EGF_EGFR: Bound ligand-receptor complex
    - pEGFR: Phosphorylated (activated) EGFR

"""

import argparse
import matplotlib.pyplot as plt
import torch
import torch.nn as nn

from torchdiffeq import odeint

torch.set_default_dtype(torch.float64)


class EGFRSignalingODE(nn.Module):
    """
    ODE system for EGF-EGFR signaling pathway.
    
    State vector: [EGF, EGFR, EGF_EGFR, pEGFR]
    """
    
    def __init__(self, k1=0.1, k2=0.05, k3=0.2, k4=0.01):
        """
        Initialize the EGFR signaling model with kinetic parameters.
        
        Args:
            k1: Binding rate constant (EGF + EGFR -> EGF*EGFR)
            k2: Unbinding rate constant (EGF*EGFR -> EGF + EGFR)
            k3: Phosphorylation rate constant (EGF*EGFR -> pEGFR)
            k4: Dephosphorylation rate constant (pEGFR -> EGFR)
        """
        super().__init__()
        # Kinetic rate constants as learnable parameters
        self.k1 = nn.Parameter(torch.tensor([k1]))  # Binding
        self.k2 = nn.Parameter(torch.tensor([k2]))  # Unbinding
        self.k3 = nn.Parameter(torch.tensor([k3]))  # Phosphorylation
        self.k4 = nn.Parameter(torch.tensor([k4]))  # Dephosphorylation

    def forward(self, t, y):
        """
        Compute derivatives for the ODE system.
        
        Args:
            t: Current time
            y: State vector [EGF, EGFR, EGF_EGFR, pEGFR]
            
        Returns:
            Derivatives [dEGF/dt, dEGFR/dt, dEGF_EGFR/dt, dpEGFR/dt]
        """
        # Unpack state variables
        EGF = y[..., 0]
        EGFR = y[..., 1]
        EGF_EGFR = y[..., 2]
        pEGFR = y[..., 3]
        
        # Ensure non-negative concentrations for stable kinetics
        EGF = torch.clamp(EGF, min=0)
        EGFR = torch.clamp(EGFR, min=0)
        EGF_EGFR = torch.clamp(EGF_EGFR, min=0)
        pEGFR = torch.clamp(pEGFR, min=0)
        
        # Reaction rates (mass action kinetics)
        # Ensure positive rate constants
        k1 = torch.abs(self.k1)
        k2 = torch.abs(self.k2)
        k3 = torch.abs(self.k3)
        k4 = torch.abs(self.k4)
        
        v1 = k1 * EGF * EGFR           # EGF + EGFR -> EGF*EGFR
        v2 = k2 * EGF_EGFR             # EGF*EGFR -> EGF + EGFR
        v3 = k3 * EGF_EGFR             # EGF*EGFR -> pEGFR
        v4 = k4 * pEGFR                # pEGFR -> EGFR

        # ODEs (mass balance)
        dEGF_dt = -v1 + v2             # EGF consumed by binding, released by unbinding
        dEGFR_dt = -v1 + v2 + v4       # EGFR consumed by binding, released by unbinding and recycling
        dEGF_EGFR_dt = v1 - v2 - v3    # Complex formed by binding, consumed by unbinding and phosphorylation
        dpEGFR_dt = v3 - v4            # pEGFR formed by phosphorylation, consumed by dephosphorylation

        return torch.stack([dEGF_dt, dEGFR_dt, dEGF_EGFR_dt, dpEGFR_dt], dim=-1)


def simulate_egfr_pathway(model, y0, t_span, n_points=100, method='dopri5'):
    """
    Simulate the EGFR signaling pathway.
    
    Args:
        model: EGFRSignalingODE instance
        y0: Initial conditions [EGF0, EGFR0, EGF_EGFR0, pEGFR0]
        t_span: Tuple of (t_start, t_end)
        n_points: Number of time points
        method: ODE solver method
        
    Returns:
        t: Time points
        y: Solution array (n_points, 4)
    """
    t = torch.linspace(t_span[0], t_span[1], n_points)
    
    with torch.no_grad():
        y = odeint(model, y0, t, method=method)
    
    return t, y


def plot_signaling_dynamics(t, y, title="EGF-EGFR Signaling Dynamics", save_path=None):
    """
    Plot the dynamics of all species in the signaling pathway.
    """
    t_np = t.numpy()
    y_np = y.numpy()
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    species = ['EGF (Ligand)', 'EGFR (Receptor)', 'EGF·EGFR (Complex)', 'pEGFR (Activated)']
    colors = ['#2ecc71', '#3498db', '#9b59b6', '#e74c3c']
    
    for i, (ax, name, color) in enumerate(zip(axes.flat, species, colors)):
        ax.plot(t_np, y_np[:, i], color=color, linewidth=2.5, label=name)
        ax.fill_between(t_np, 0, y_np[:, i], color=color, alpha=0.2)
        ax.set_xlabel('Time (min)', fontsize=12)
        ax.set_ylabel('Concentration (nM)', fontsize=12)
        ax.set_title(name, fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([t_np[0], t_np[-1]])
        ax.set_ylim(bottom=0)
        
        # Style the plot
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, axes


def plot_combined_dynamics(t, y, save_path=None):
    """
    Plot all species on a single figure for comparison.
    """
    t_np = t.numpy()
    y_np = y.numpy()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    species = ['EGF', 'EGFR', 'EGF·EGFR', 'pEGFR']
    colors = ['#2ecc71', '#3498db', '#9b59b6', '#e74c3c']
    linestyles = ['-', '--', '-.', ':']
    
    for i, (name, color, ls) in enumerate(zip(species, colors, linestyles)):
        ax.plot(t_np, y_np[:, i], color=color, linewidth=2.5, 
                label=name, linestyle=ls)
    
    ax.set_xlabel('Time (min)', fontsize=14)
    ax.set_ylabel('Concentration (nM)', fontsize=14)
    ax.set_title('EGF-EGFR Signaling Pathway Dynamics', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([t_np[0], t_np[-1]])
    ax.set_ylim(bottom=0)
    
    # Style the plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, ax


def main():
    parser = argparse.ArgumentParser(description='EGF-EGFR Signaling Pathway ODE Model')
    parser.add_argument('--egf0', type=float, default=10.0, 
                        help='Initial EGF concentration (nM)')
    parser.add_argument('--egfr0', type=float, default=100.0, 
                        help='Initial EGFR concentration (nM)')
    parser.add_argument('--t_end', type=float, default=60.0, 
                        help='Simulation end time (min)')
    parser.add_argument('--n_points', type=int, default=200, 
                        help='Number of time points')
    parser.add_argument('--method', type=str, default='dopri5',
                        choices=['dopri5', 'rk4', 'euler', 'midpoint'],
                        help='ODE solver method')
    parser.add_argument('--save', action='store_true', 
                        help='Save figures')
    args = parser.parse_args()
    
    print("=" * 60)
    print("EGF-EGFR Signaling Pathway Model")
    print("=" * 60)
    print("\nReaction Scheme:")
    print("  EGF + EGFR <-> EGF·EGFR -> pEGFR")
    print("\nInitial Conditions:")
    print(f"  EGF:      {args.egf0:.1f} nM")
    print(f"  EGFR:     {args.egfr0:.1f} nM")
    print(f"  EGF·EGFR: 0.0 nM")
    print(f"  pEGFR:    0.0 nM")
    print(f"\nSimulation: 0 to {args.t_end:.1f} min ({args.n_points} points)")
    print(f"Solver: {args.method}")
    print("=" * 60)
    
    # Create model with default kinetic parameters
    model = EGFRSignalingODE(
        k1=0.1,    # Binding rate
        k2=0.05,   # Unbinding rate
        k3=0.2,    # Phosphorylation rate
        k4=0.01    # Dephosphorylation rate
    )
    
    print("\nKinetic Parameters:")
    print(f"  k1 (binding):         {model.k1.item():.4f} nM⁻¹·min⁻¹")
    print(f"  k2 (unbinding):       {model.k2.item():.4f} min⁻¹")
    print(f"  k3 (phosphorylation): {model.k3.item():.4f} min⁻¹")
    print(f"  k4 (dephosphorylation): {model.k4.item():.4f} min⁻¹")
    
    # Initial conditions: [EGF, EGFR, EGF_EGFR, pEGFR]
    y0 = torch.tensor([args.egf0, args.egfr0, 0.0, 0.0])
    
    # Simulate
    print("\nRunning simulation...")
    t, y = simulate_egfr_pathway(
        model, y0, 
        t_span=(0.0, args.t_end), 
        n_points=args.n_points,
        method=args.method
    )
    
    # Report final concentrations
    print("\nFinal Concentrations:")
    print(f"  EGF:      {y[-1, 0].item():.4f} nM")
    print(f"  EGFR:     {y[-1, 1].item():.4f} nM")
    print(f"  EGF·EGFR: {y[-1, 2].item():.4f} nM")
    print(f"  pEGFR:    {y[-1, 3].item():.4f} nM")
    
    # Conservation check
    total_egf = y[-1, 0] + y[-1, 2]  # Free EGF + bound EGF
    total_receptor = y[-1, 1] + y[-1, 2] + y[-1, 3]  # Free + bound + phosphorylated
    print(f"\nConservation Check:")
    print(f"  Total EGF (should be ~{args.egf0:.1f}):   {total_egf.item():.4f} nM")
    print(f"  Total EGFR (should be ~{args.egfr0:.1f}): {total_receptor.item():.4f} nM")
    
    # Plot results
    save_path1 = 'egfr_signaling_subplots.png' if args.save else None
    save_path2 = 'egfr_signaling_combined.png' if args.save else None
    
    plot_signaling_dynamics(t, y, save_path=save_path1)
    plot_combined_dynamics(t, y, save_path=save_path2)
    
    plt.show()
    
    print("\nSimulation complete!")
    

if __name__ == "__main__":
    main()
