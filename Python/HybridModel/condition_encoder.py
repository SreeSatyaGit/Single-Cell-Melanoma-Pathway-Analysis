import torch
from torch import nn

class ConditionEncoder(nn.Module):
    """
    Learns embeddings for each drug condition to capture:
    - Drug-specific effects (Tram vs Vem)
    - Synergistic/antagonistic interactions
    - Temporal resistance patterns
    """
    def __init__(self, num_conditions=3, embedding_dim=16):
        super().__init__()
        self.embeddings = nn.Parameter(torch.randn(num_conditions, embedding_dim) * 0.1)
        
    def forward(self, condition_idx):
        """
        Args:
            condition_idx: Integer tensor [Batch] with values 0,1,2 for Tram/Vem/Combo
        Returns:
            Condition embedding [Batch, embedding_dim]
        """
        return self.embeddings[condition_idx]
    
    def get_all_embeddings(self):
        """Returns all condition embeddings for analysis"""
        return self.embeddings
