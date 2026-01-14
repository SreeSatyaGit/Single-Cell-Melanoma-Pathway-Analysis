# Notebook Visualization Code

Copy and paste these cells into your Jupyter notebook to visualize the trained model predictions.

---

## Cell 1: Setup and Import

```python
import os
os.chdir('/Users/bharadwajanandivada/SCMPA/Python/HybridModel')

import torch
import matplotlib.pyplot as plt
import numpy as np
from visualize_trained_model import visualize_trained_model

# Set matplotlib style for better plots
plt.style.use('seaborn-v0_8-darkgrid')
%matplotlib inline
```

---

## Cell 2: Run Full Visualization (Easiest - One Line!)

```python
# This runs everything: loads model, makes predictions, generates all plots
model, predictions, metrics = visualize_trained_model()
```

---

## Cell 3: Manual Step-by-Step (For More Control)

```python
from visualize_trained_model import load_trained_model, run_predictions, calculate_metrics, create_publication_plots

# Load model
model, y0_baseline = load_trained_model()

# Run predictions
predictions, data, t_span_hours = run_predictions(model, y0_baseline)

# Calculate metrics
metrics, overall_mse = calculate_metrics(predictions, data)

# Create plots
fig = create_publication_plots(predictions, data, t_span_hours)
plt.show()

print(f"\n📊 Overall MSE: {overall_mse:.6f}")
```

---

## Cell 4: Create Custom Plots for Individual Conditions

```python
import torch
import matplotlib.pyplot as plt
from windowed_hybrid import WindowedHybridModel
from data_loader import get_experimental_data
from configs import get_baseline_configs

# Load model
p_baseline, y0_baseline = get_baseline_configs()
model = WindowedHybridModel(p_baseline)
checkpoint = torch.load('windowed_hybrid_trained.pth')
model.load_state_dict(checkpoint['model_state_dict'])
model.eval()

# Get data
data = get_experimental_data()
t_span_hours = data['times'] / 3600.0

# Pick a condition (0=Tram, 1=Vem, 2=Combo)
condition_name = "Combination"
drug_conc = torch.tensor([[1.0, 1.0, 0.0, 0.0]])  # Tram + Vem
condition_idx = torch.tensor([2])

# Run prediction
with torch.no_grad():
    y0 = y0_baseline.unsqueeze(0)
    output = model(y0, t_span_hours, drug_conc, condition_idx)
    sensors = model.physics.get_sensors(output.squeeze(0))

# Plot specific proteins
proteins_to_plot = ['pERK', 'pMEK', 'pAKT', 'DUSP']
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, protein in enumerate(proteins_to_plot):
    ax = axes[idx]
    
    if protein in sensors:
        # Model prediction
        pred = sensors[protein].detach().numpy()
        time = t_span_hours.numpy()
        
        ax.plot(time, pred, '-', linewidth=3, label='Hybrid Model', color='#2ecc71')
        
        # Experimental data
        if protein in data['combo']:
            exp_data = data['combo'][protein].numpy()
            ax.scatter(time, exp_data, s=100, label='Experimental', 
                      color='#e74c3c', edgecolor='white', linewidth=2, zorder=10)
        
        ax.set_title(f'{protein} - {condition_name}', fontweight='bold', fontsize=13)
        ax.set_xlabel('Time (hours)', fontsize=11)
        ax.set_ylabel('Normalized Activity', fontsize=11)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-2, 50)
        ax.set_ylim(-0.1, 1.1)

plt.tight_layout()
plt.savefig('outputs/custom_plot.png', dpi=300, bbox_inches='tight')
plt.show()
```

---

## Cell 5: Compare All Three Conditions Side-by-Side

```python
from visualize_trained_model import load_trained_model, run_predictions

# Load and run predictions
model, y0_baseline = load_trained_model()
predictions, data, t_span_hours = run_predictions(model, y0_baseline)

# Focus on one protein across all conditions
protein = 'pERK'
time = t_span_hours.numpy()

fig, ax = plt.subplots(figsize=(12, 6))

colors = {'Trametinib': '#e74c3c', 'Vemurafenib': '#3498db', 'Combination': '#2ecc71'}
data_keys = {'Trametinib': 'tram', 'Vemurafenib': 'vem', 'Combination': 'combo'}

for cond_name, pred_info in predictions.items():
    sensors = pred_info['sensors']
    data_key = pred_info['data_key']
    color = colors[cond_name]
    
    if protein in sensors:
        # Predictions
        pred = sensors[protein].detach().numpy()
        ax.plot(time, pred, '-', linewidth=3, label=f'{cond_name} (Model)', 
               color=color, alpha=0.8)
        
        # Experimental
        if protein in data[data_key]:
            exp = data[data_key][protein].numpy()
            ax.scatter(time, exp, s=120, color=color, edgecolor='white', 
                      linewidth=2, zorder=10, label=f'{cond_name} (Exp)')

ax.set_title(f'{protein} Response Across Treatment Conditions', 
            fontweight='bold', fontsize=15)
ax.set_xlabel('Time (hours)', fontsize=13)
ax.set_ylabel('Normalized Activity', fontsize=13)
ax.legend(fontsize=11, loc='best', frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(-2, 50)
ax.set_ylim(-0.1, 1.1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(f'outputs/{protein}_all_conditions.png', dpi=300, bbox_inches='tight')
plt.show()
```

---

## Cell 6: Calculate and Display Improvement Over Baseline

```python
from visualize_trained_model import load_trained_model, run_predictions
from configs import get_baseline_configs
from model import MAPKPhysics
from torchdiffeq import odeint

# Load trained model
model, y0_baseline = load_trained_model()
predictions_hybrid, data, t_span_hours = run_predictions(model, y0_baseline)

# Run baseline physics
p_baseline, _ = get_baseline_configs()
physics_baseline = MAPKPhysics(p_baseline)

conditions_map = {
    'Trametinib': (torch.tensor([[1.0, 0.0, 0.0, 0.0]]), 'tram'),
    'Vemurafenib': (torch.tensor([[0.0, 1.0, 0.0, 0.0]]), 'vem'),
    'Combination': (torch.tensor([[1.0, 1.0, 0.0, 0.0]]), 'combo')
}

print("=" * 70)
print("BASELINE vs HYBRID MODEL COMPARISON")
print("=" * 70)

for cond_name, (drug, data_key) in conditions_map.items():
    # Baseline prediction
    with torch.no_grad():
        traj_baseline = odeint(
            lambda t, y: physics_baseline(t, y, drug),
            y0_baseline.unsqueeze(0), t_span_hours,
            method='midpoint', options={'step_size': 0.05}
        )
        sensors_baseline = physics_baseline.get_sensors(traj_baseline.squeeze(0))
    
    # Calculate MSE for baseline
    baseline_mse = 0
    hybrid_mse = 0
    count = 0
    
    for protein in ['pERK', 'pMEK', 'pAKT', 'DUSP', 'pCRAF']:
        if protein in sensors_baseline and protein in data[data_key]:
            # Baseline
            pred_base = sensors_baseline[protein]
            target = data[data_key][protein]
            baseline_mse += torch.mean((pred_base - target)**2).item()
            
            # Hybrid
            pred_hybrid = predictions_hybrid[cond_name]['sensors'][protein]
            hybrid_mse += torch.mean((pred_hybrid - target)**2).item()
            
            count += 1
    
    if count > 0:
        baseline_mse /= count
        hybrid_mse /= count
        improvement = ((baseline_mse - hybrid_mse) / baseline_mse) * 100
        
        print(f"\n{cond_name}:")
        print(f"  Baseline MSE:  {baseline_mse:.6f}")
        print(f"  Hybrid MSE:    {hybrid_mse:.6f}")
        print(f"  Improvement:   {improvement:+.2f}%")

print("\n" + "=" * 70)
```

---

## Quick Start (Copy This!)

For the fastest results, just run this single cell:

```python
import os
os.chdir('/Users/bharadwajanandivada/SCMPA/Python/HybridModel')

from visualize_trained_model import visualize_trained_model
model, predictions, metrics = visualize_trained_model()
```

All plots will be saved to the `outputs/` directory!
