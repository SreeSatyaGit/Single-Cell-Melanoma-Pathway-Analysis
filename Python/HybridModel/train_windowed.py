import torch
from torch import optim
from data_loader import get_experimental_data
from configs import get_baseline_configs
from windowed_hybrid import WindowedHybridModel
import matplotlib.pyplot as plt
import time
from pathlib import Path

def train_windowed_model(epochs=200, lr=1e-3):
    """
    Train the windowed hybrid model with stable gradients.
    
    Key innovation: Gradients only flow through 8h windows + NN corrections,
    NOT through full 48h ODE integration.
    """
    print("="*70)
    print("WINDOWED HYBRID MODEL TRAINING")
    print("="*70)
    
    # Setup
    p_baseline, y0_baseline = get_baseline_configs()
    model = WindowedHybridModel(p_baseline)
    
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Trainable parameters: {trainable_params}")
    print(f"  - Condition embeddings: {model.condition_encoder.embeddings.numel()}")
    print(f"  - Correction network: {sum(p.numel() for p in model.correction_net.parameters())}")
    
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    # Data
    data = get_experimental_data()
    t_span_hours = data['times'] / 3600.0
    
    # Prepare batched data
    condition_map = {'tram': 0, 'vem': 1, 'combo': 2}
    batch_y0 = []
    batch_drugs = []
    batch_conditions = []
    batch_targets = []
    
    species_keys = ['pEGFR', 'pCRAF', 'pMEK', 'pERK', 'DUSP', 'pAKT', 'p4ebp1', 'her2', 'her3', 'pDGFR', 'pS6k', 'panRAS']
    
    for cond_name, cond_idx in condition_map.items():
        batch_y0.append(y0_baseline)
        
        # Drug concentrations
        if cond_name == 'tram':
            batch_drugs.append(torch.tensor([1.0, 0.0, 0.0, 0.0]))
        elif cond_name == 'vem':
            batch_drugs.append(torch.tensor([0.0, 1.0, 0.0, 0.0]))
        else:  # combo
            batch_drugs.append(torch.tensor([1.0, 1.0, 0.0, 0.0]))
        
        batch_conditions.append(cond_idx)
        
        # Target: [num_times, num_sensors]
        target = torch.stack([data[cond_name][s] for s in species_keys], dim=-1)  # [6, 12]
        batch_targets.append(target)
    
    batch_y0 = torch.stack(batch_y0)  # [3, 68]
    batch_drugs = torch.stack(batch_drugs)  # [3, 4]
    batch_conditions = torch.tensor(batch_conditions, dtype=torch.long)  # [3]
    batch_targets = torch.stack(batch_targets)  # [3, 6, 12]
    
    # Training loop
    print("\n" + "="*70)
    print("TRAINING")
    print("="*70)
    
    loss_history = []
    start_time = time.time()
    
    for epoch in range(epochs):
        optimizer.zero_grad()
        
        try:
            # Forward pass: [3, 6, 68]
            predictions = model(batch_y0, t_span_hours, batch_drugs, batch_conditions)
            
            # Extract sensors
            sensors_pred = model.physics.get_sensors(predictions.permute(1, 0, 2))  # Needs [T, B, 68]
            
            # Stack predictions: [3, 6, 12]
            preds_stacked = torch.stack([sensors_pred[s].T for s in species_keys], dim=-1)
            
            # MSE loss
            loss = torch.mean((preds_stacked - batch_targets)**2)
            
            if torch.isnan(loss) or torch.isinf(loss):
                print(f"\n✗ NaN/Inf at epoch {epoch}")
                print(f"  Correction scale: {model.correction_net.correction_scale.item():.4f}")
                break
            
            loss.backward()
            
            # Moderate gradient clipping
            grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            
            optimizer.step()
            loss_history.append(loss.item())
            
            if epoch % 20 == 0:
                corr_scale = model.correction_net.correction_scale.item()
                print(f"Epoch {epoch:03d} | Loss: {loss.item():.6f} | Grad: {grad_norm:.2f} | CorrScale: {corr_scale:.4f}")
        
        except Exception as e:
            print(f"\n✗ Exception at epoch {epoch}: {e}")
            import traceback
            traceback.print_exc()
            break
    
    print("\n" + "="*70)
    print("TRAINING COMPLETE")
    print("="*70)
    print(f"Final Loss: {loss.item():.6f}")
    print(f"Baseline Loss: 0.294")
    print(f"Improvement: {(0.294 - loss.item())/0.294 * 100:.1f}%")
    print(f"Total Time: {time.time() - start_time:.1f}s")
    
    # Create output directory
    output_dir = Path(__file__).parent / 'outputs'
    output_dir.mkdir(exist_ok=True)
    
    # Save model
    model_path = output_dir / 'windowed_hybrid_trained.pth'
    torch.save({
        'model_state_dict': model.state_dict(),
        'loss_history': loss_history,
        'final_loss': loss.item(),
        'baseline_loss': 0.294,
        'improvement_pct': (0.294 - loss.item())/0.294 * 100,
        'training_time': time.time() - start_time
    }, model_path)
    print(f"Model saved to: {model_path}")
    
    # Plot
    if len(loss_history) > 1:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        
        axes[0].plot(loss_history, linewidth=2)
        axes[0].axhline(y=0.294, color='r', linestyle='--', label='Baseline')
        axes[0].set_xlabel('Epoch')
        axes[0].set_ylabel('MSE Loss')
        axes[0].set_title('Training Progress')
        axes[0].grid(True, alpha=0.3)
        axes[0].legend()
        
        axes[1].plot(loss_history, linewidth=2)
        axes[1].axhline(y=0.294, color='r', linestyle='--', label='Baseline')
        axes[1].set_yscale('log')
        axes[1].set_xlabel('Epoch')
        axes[1].set_ylabel('MSE Loss (log)')
        axes[1].set_title('Training Progress (Log Scale)')
        axes[1].grid(True, alpha=0.3)
        axes[1].legend()
        
        plt.tight_layout()
        plot_path = output_dir / 'windowed_training.png'
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Training curve saved to: {plot_path}")

if __name__ == "__main__":
    train_windowed_model(epochs=200, lr=1e-3)
