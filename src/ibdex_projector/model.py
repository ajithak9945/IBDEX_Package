from __future__ import annotations

import json

import torch
from torch import nn

from .artifacts import artifact_path


class IbdexCvae(nn.Module):
    def __init__(self, input_dim: int, encoder_condition_dim: int, decoder_condition_dim: int,
                 latent_dim: int = 16, hidden_dim: int = 512, depth: int = 2, dropout: float = 0.1):
        super().__init__()
        enc = []
        in_dim = input_dim + encoder_condition_dim
        for _ in range(depth):
            enc.extend([nn.Linear(in_dim, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(), nn.Dropout(dropout)])
            in_dim = hidden_dim
        self.encoder = nn.Sequential(*enc)
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)

        self.dec_blocks = nn.ModuleList()
        in_dim = latent_dim + decoder_condition_dim
        for _ in range(depth):
            self.dec_blocks.append(nn.Sequential(
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.GELU(),
                nn.Dropout(dropout),
            ))
            in_dim = hidden_dim
        self.dec_out = nn.Linear(hidden_dim, input_dim)

        self.dataset_adv_head = nn.Sequential(nn.Linear(latent_dim, 32), nn.GELU(), nn.Dropout(dropout), nn.Linear(32, 4))
        self.dataset_adv_linear_head = nn.Linear(latent_dim, 4)
        self.tissue_adv_head = nn.Sequential(nn.Linear(latent_dim, 32), nn.GELU(), nn.Dropout(dropout), nn.Linear(32, 7))

    def encode(self, x: torch.Tensor, c_encoder: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.encoder(torch.cat([x, c_encoder], dim=1))
        return self.fc_mu(h), torch.clamp(self.fc_logvar(h), min=-10.0, max=10.0)

    def decode(self, z: torch.Tensor, c_decoder: torch.Tensor) -> torch.Tensor:
        h = torch.cat([z, c_decoder], dim=1)
        for block in self.dec_blocks:
            h = block(h)
        return self.dec_out(h)


def load_model(device: str = "cpu") -> tuple[IbdexCvae, dict]:
    config = json.loads(artifact_path("model_config.json").read_text())
    params = config["best_params"]
    model = IbdexCvae(
        input_dim=config["input_dim"],
        encoder_condition_dim=config["encoder_condition_dim"],
        decoder_condition_dim=config["decoder_condition_dim"],
        latent_dim=params["latent_dim"],
        hidden_dim=params["hidden_dim"],
        depth=params["depth"],
        dropout=params["dropout"],
    )
    state = torch.load(artifact_path("model_final_train_val.pt"), map_location=device)
    model.load_state_dict(state)
    model.to(device)
    model.eval()
    return model, config
