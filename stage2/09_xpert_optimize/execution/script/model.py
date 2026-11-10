"""XPert with a frozen continuous-expression identity and shared Morgan input."""

import torch
from torch import nn
from torch.nn import functional as F


CHANNELS = 256
HEADS = 8
DROPOUT = 0.1
LOSS_WEIGHTS = (0.2, 0.003, 0.2, 1.0)
MLP_WIDTHS = {f"XPert_MLP{width}": width for width in (512, 1024, 2048, 4096)}


class Attention(nn.Module):
    def __init__(self):
        super().__init__()
        self.query = nn.Linear(CHANNELS, CHANNELS)
        self.key = nn.Linear(CHANNELS, CHANNELS)
        self.value = nn.Linear(CHANNELS, CHANNELS)

    def forward(self, query, context):
        def heads(value):
            return value.unflatten(-1, (HEADS, CHANNELS // HEADS)).transpose(1, 2)

        value = F.scaled_dot_product_attention(
            heads(self.query(query)), heads(self.key(context)), heads(self.value(context)),
            dropout_p=DROPOUT if self.training else 0.0,
        )
        return value.transpose(1, 2).flatten(2)


class Residual(nn.Module):
    def __init__(self):
        super().__init__()
        self.norm = nn.LayerNorm(CHANNELS, eps=1e-12)
        self.dropout = nn.Dropout(DROPOUT)

    def forward(self, update, tokens):
        return tokens + self.dropout(self.norm(update))


class SelfBlock(nn.Module):
    def __init__(self):
        super().__init__()
        self.attention = Attention()
        self.attention_residual = Residual()
        self.feedforward = nn.Sequential(nn.Linear(CHANNELS, CHANNELS * 2), nn.ReLU(),
                                         nn.Dropout(DROPOUT), nn.Linear(CHANNELS * 2, CHANNELS))
        self.output_residual = Residual()

    def attend(self, tokens):
        return self.attention_residual(self.attention(tokens, tokens), tokens)

    def output(self, tokens):
        return self.output_residual(self.feedforward(tokens), tokens)

    def forward(self, tokens):
        return self.output(self.attend(tokens))


class CrossBlock(SelfBlock):
    """Author CA layer: drug SA, gene SA, then gene-to-drug attention."""

    def __init__(self):
        super().__init__()
        self.drug_attention = SelfBlock()
        self.cross_attention = Attention()

    def forward(self, genes, drug):
        drug = self.drug_attention(drug)
        genes = self.attend(genes)
        genes = self.attention_residual(self.cross_attention(genes, drug), genes)
        return self.output(genes), drug


class TokenMap(nn.Sequential):
    """Apply the same two-layer map along the token axis of every channel."""

    def __init__(self, input_tokens, hidden_tokens, output_tokens):
        super().__init__(nn.Linear(input_tokens, hidden_tokens), nn.GELU(),
                         nn.Linear(hidden_tokens, output_tokens))

    def forward(self, matrix):
        mapped = super().forward(matrix.transpose(1, 2).contiguous())
        return mapped.transpose(1, 2).contiguous()


class XPert(nn.Module):
    def __init__(self, identity, latent_tokens=None):
        super().__init__()
        if identity.ndim != 2 or identity.shape[1] != CHANNELS:
            raise ValueError("identity must have one 256-channel vector per gene")
        self.register_buffer("identity", identity, persistent=False)
        self.expression = nn.Linear(1, CHANNELS)
        self.morgan = nn.Linear(2048, CHANNELS)
        self.dose = nn.Embedding(10, CHANNELS)
        self.time = nn.Embedding(3, CHANNELS)
        self.drug_positions = nn.Parameter(torch.empty(3, CHANNELS))
        nn.init.normal_(self.drug_positions)
        self.cell_input = nn.Sequential(nn.LayerNorm(CHANNELS, eps=1e-12), nn.Dropout(DROPOUT))
        self.drug_input = nn.Sequential(nn.LayerNorm(CHANNELS, eps=1e-12), nn.Dropout(DROPOUT))

        self.treatment = nn.ModuleList([CrossBlock(), SelfBlock(), SelfBlock(), CrossBlock()])
        self.baseline = nn.Sequential(*(SelfBlock() for _ in range(4)))
        self.treated_head = self.make_head()
        self.baseline_head = self.make_head()
        self.delta_head = self.make_head()

        # Create arm-specific parameters after the shared XPert initialization.
        self.latent_tokens = latent_tokens
        genes = len(identity)
        self.encoder = TokenMap(genes, latent_tokens, latent_tokens) if latent_tokens else nn.Identity()
        self.decoder = TokenMap(latent_tokens, latent_tokens, genes) if latent_tokens else nn.Identity()

    @staticmethod
    def make_head():
        return nn.Sequential(nn.Linear(CHANNELS, 64), nn.ReLU(), nn.Dropout(DROPOUT), nn.Linear(64, 1))

    def forward(self, control, morgan, dose):
        cell = self.expression(control.unsqueeze(-1))
        genes = self.encoder(self.cell_input(cell + self.identity.to(cell.dtype)))
        drug = torch.stack((self.dose(dose), self.time(torch.full_like(dose, 2)), self.morgan(morgan)), dim=1)
        drug = self.drug_input(drug + self.drug_positions)
        treated, drug = self.treatment[0](genes, drug)
        treated = self.treatment[1](treated)
        treated = self.treatment[2](treated)
        treated, _ = self.treatment[3](treated, drug)
        treated = self.decoder(treated)
        baseline = self.decoder(self.baseline(genes))
        return (self.treated_head(treated).squeeze(-1), self.baseline_head(baseline).squeeze(-1),
                self.delta_head(treated - baseline).squeeze(-1))

    def architecture(self):
        return {
            "gene_count": len(self.identity),
            "channels": CHANNELS,
            "attention_heads": HEADS,
            "attention_tokens": self.latent_tokens or len(self.identity),
            "latent_tokens": self.latent_tokens,
            "token_mlp_hidden": self.latent_tokens,
            "token_mlp_activation": "GELU" if self.latent_tokens else None,
            "parameter_count": sum(p.numel() for p in self.parameters()),
            "compression_parameter_count": sum(p.numel() for m in (self.encoder, self.decoder)
                                               for p in m.parameters()),
        }


def loss_terms(prediction, treated, control):
    """Per-profile gene means and PCC, averaged equally across the batch."""
    expression, reconstruction, delta = (value.float() for value in prediction)
    target = treated - control
    predicted_centered = delta - delta.mean(dim=1, keepdim=True)
    target_centered = target - target.mean(dim=1, keepdim=True)
    denominator = (predicted_centered.square().sum(dim=1) * target_centered.square().sum(dim=1) + 1e-6).sqrt()
    pearson_loss = 1 - (predicted_centered * target_centered).sum(dim=1) / denominator
    return torch.stack((F.mse_loss(expression, treated), F.mse_loss(reconstruction, control),
                        F.mse_loss(delta, target), pearson_loss.mean()))
