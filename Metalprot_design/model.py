"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This file contains neural net classes from Metalprot_learning
"""

#imports
import torch

class DoubleLayerNet(torch.nn.Module):
    def __init__(self, input_dim: int, l1: int, l2: int, l3: int, output_dim: int, input_dropout: float, hidden_dropout: float):
        super(DoubleLayerNet, self).__init__()

        self.block1 = torch.nn.Sequential(
            torch.nn.Linear(input_dim, l1),
            torch.nn.Dropout(input_dropout),
            torch.nn.ReLU(),
            torch.nn.Linear(l1,l2),
            torch.nn.Dropout(hidden_dropout),
            torch.nn.ReLU(),
            torch.nn.Linear(l2,l3),
            torch.nn.Dropout(hidden_dropout),
            torch.nn.ReLU(),
            torch.nn.Linear(l3,output_dim),
            torch.nn.ReLU()
        )

    def forward(self, x):
        y = self.block1(x.float())
        return y 