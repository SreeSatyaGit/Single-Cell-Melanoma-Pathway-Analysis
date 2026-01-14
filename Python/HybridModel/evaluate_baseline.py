import torch
from data_loader import get_experimental_data
from configs import get_baseline_configs
from model import MAPKPhysics
from torchdiffeq import odeint
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path

def evaluate_baseline():
    """
    Evaluate pure MATLAB-calibrated physics on the 3 conditions.
    No training, just assessment.
    """
    print("="*70)
    print("BASELINE PHYSICS EVALUATION")
    print("="*70)
    
    # Setup
    p_baseline, y0_baseline = get_baseline_configs()
    physics = MAPKPhysics(p_baseline)
    
    # Data
    data = get_experimental_data()
    t_span = data['times'] / 3600.0
    
    conditions = {
        'Trametinib': torch.tensor([1.0, 0.0, 0.0, 0.0]),
        'Vemurafenib': torch.tensor([0.0, 1.0, 0.0, 0.0]),
        'Combination': torch.tensor([1.0, 1.0, 0.0, 0.0])
    }
    
    # Simulate all 3 conditions
    batch_drugs = torch.stack(list(conditions.values()))
    batch_y0 = y0_baseline.repeat(3, 1)
    
    print("\nRunning 48-hour simulations...")
    yhat = odeint(lambda t, y: physics(t, y, batch_drugs), 
                   batch_y0, t_span, method='midpoint', options={'step_size': 0.05})
    print("✓ Simulations complete!")
    
    # Extract sensors
    sensors = physics.get_sensors(yhat)  # [6, 3] for each protein
    
    # Calculate condition-specific losses
    species_keys = ['pEGFR', 'pCRAF', 'pMEK', 'pERK', 'DUSP', 'pAKT', 'p4ebp1', 'her2', 'her3', 'pDGFR', 'pS6k', 'panRAS']
    condition_names = ['tram', 'vem', 'combo']
    
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    
    overall_loss = 0
    for cond_idx, cond_name in enumerate(condition_names):
        cond_loss = 0
        for species in species_keys:
            if species in sensors and species in data[cond_name]:
                pred = sensors[species][:, cond_idx]  # [6]
                target = data[cond_name][species]
                mse = torch.mean((pred - target)**2).item()
                cond_loss += mse
        
        avg_loss = cond_loss / len(species_keys)
        overall_loss += avg_loss
        print(f"{list(conditions.keys())[cond_idx]:15s}: MSE = {avg_loss:.6f}")
    
    print(f"\nOverall Average MSE: {overall_loss / 3:.6f}")
    
    # Visualization
    print("\n" + "="*70)
    print("GENERATING PLOTS")
    print("="*70)
    
    fig, axes = plt.subplots(3, 4, figsize=(16, 10))
    axes = axes.flatten()
    
    key_sensors = ['pERK', 'pMEK', 'pCRAF', 'pAKT', 'DUSP', 'panRAS', 'p4ebp1', 'pEGFR', 'her2', 'her3', 'pDGFR', 'pS6k']
    colors = ['#e74c3c', '#3498db', '#2ecc71']
    condition_labels = list(conditions.keys())
    
    for idx, sensor in enumerate(key_sensors[:12]):
        ax = axes[idx]
        
        if sensor in sensors:
            time_hrs = t_span.numpy()
            
            for cond_idx, cond_name in enumerate(condition_names):
                # Predicted
                pred = sensors[sensor][:, cond_idx].detach().numpy()
                ax.plot(time_hrs, pred, '-', color=colors[cond_idx], linewidth=2, 
                       label=f'{condition_labels[cond_idx]} (Model)')
                
                # Experimental
                if sensor in data[cond_name]:
                    target = data[cond_name][sensor].numpy()
                    ax.scatter(time_hrs, target, color=colors[cond_idx], s=60, 
                             marker='o', edgecolor='white', linewidth=1.5, zorder=10)
            
            ax.set_title(sensor, fontweight='bold', fontsize=10)
            ax.set_xlabel('Time (hours)', fontsize=8)
            ax.set_ylabel('Normalized Activity', fontsize=8)
            ax.grid(True, alpha=0.2)
            ax.set_xlim(-2, 50)
            ax.set_ylim(-0.1, 1.1)
            
            if idx == 0:
                ax.legend(fontsize=7, loc='upper right')
    
    plt.tight_layout()
    
    # Create output directory
    output_dir = Path(__file__).parent / 'outputs'
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / 'baseline_physics_evaluation.png'
    
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"✓ Plots saved to: {output_path}")
    
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    print("The plots show:")
    print("  • Lines = Model predictions (pure MATLAB physics)")
    print("  • Dots  = Your experimental Western Blot data")
    print("\nA good fit means the MATLAB model already captures the biology well.")
    print("Discrepancies indicate where a Neural Network could add value later.")
    print("="*70)

if __name__ == "__main__":
    evaluate_baseline()
