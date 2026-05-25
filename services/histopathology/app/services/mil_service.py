"""Attention-based Multiple Instance Learning (ABMIL) aggregator.

A WSI is a "bag" of many patch features with no per-patch labels -- only a
slide-level label (mutation present / absent). ABMIL (Ilse et al., 2018) learns
attention weights over patches, pools them into a single slide-level embedding,
and classifies per gene.

IMPORTANT: there are NO trained weights here. The module runs with
randomly-initialized (but DETERMINISTICALLY seeded) weights, so predictions are
illustrative scaffolding only. main.py always appends a warning to that effect.

SCAFFOLD markers throughout indicate where trained ABMIL weights / a real
classification head would replace the random init.
"""
import logging
from typing import Dict, List, Any

import numpy as np

logger = logging.getLogger(__name__)

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    _TORCH_AVAILABLE = True
except Exception as e:  # pragma: no cover
    torch = None
    nn = None
    F = None
    _TORCH_AVAILABLE = False
    logger.warning("torch not available (%s); MIL aggregation disabled.", e)

from .cnn_service import FEATURE_DIM

_SEED = 42

# Genes the (scaffold) classification head is built for. Requested genes outside
# this set still get a deterministic prediction via a per-gene hash offset.
SUPPORTED_GENES = ["TP53", "IDH1", "KRAS"]


if _TORCH_AVAILABLE:

    class ABMIL(nn.Module):
        """Attention-based MIL: gated attention pooling + per-gene heads.

        Architecture (Ilse et al. gated attention):
          patch feats (N, D) -> embed (N, H)
          attention: a = softmax( w^T ( tanh(V h) * sigmoid(U h) ) )  over N
          slide embedding z = sum_i a_i * h_i    (1, H)
          per-gene head: sigmoid(linear(z)) -> P(mutated) per gene
        """

        def __init__(self, in_dim: int = FEATURE_DIM, hidden: int = 256, n_genes: int = 3):
            super().__init__()
            self.embed = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU())
            # Gated attention branches.
            self.att_V = nn.Linear(hidden, 128)
            self.att_U = nn.Linear(hidden, 128)
            self.att_w = nn.Linear(128, 1)
            # SCAFFOLD: one logit per gene. Trained weights would replace this.
            self.classifier = nn.Linear(hidden, n_genes)

        def forward(self, feats):  # feats: (N, in_dim)
            h = self.embed(feats)                       # (N, H)
            a = self.att_w(torch.tanh(self.att_V(h)) * torch.sigmoid(self.att_U(h)))  # (N,1)
            a = torch.softmax(a, dim=0)                 # attention over patches
            z = torch.sum(a * h, dim=0, keepdim=True)   # (1, H) slide embedding
            logits = self.classifier(z)                 # (1, n_genes)
            probs = torch.sigmoid(logits)               # (1, n_genes)
            return probs.squeeze(0), a.squeeze(1)       # (n_genes,), (N,)


class MILService:
    """Wraps the ABMIL module and produces per-gene predictions."""

    def __init__(self):
        self.available = _TORCH_AVAILABLE
        if not self.available:
            self.device = None
            self.model = None
            return

        torch.manual_seed(_SEED)  # deterministic random init for stable demos
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        # SCAFFOLD: random-init ABMIL. Load trained weights here in production:
        #   self.model.load_state_dict(torch.load(MIL_WEIGHTS, map_location=...))
        self.model = ABMIL(in_dim=FEATURE_DIM, hidden=256, n_genes=len(SUPPORTED_GENES))
        self.model.to(self.device)
        self.model.eval()
        logger.info("MILService initialized (UNTRAINED scaffold) on %s", self.device)

    def _gene_prob(self, base_probs: "np.ndarray", gene: str) -> float:
        """Map a requested gene to a probability.

        Supported genes read directly from the head's output. Unsupported genes
        get a deterministic value derived from a stable hash of the gene name so
        the demo is reproducible without inventing a real signal.
        """
        if gene in SUPPORTED_GENES:
            return float(base_probs[SUPPORTED_GENES.index(gene)])
        # Deterministic pseudo-probability for an out-of-set gene.
        seed = sum(ord(c) for c in gene)
        return float((seed % 1000) / 1000.0)

    def predict(self, features: "np.ndarray", target_mutations: List[str]) -> Dict[str, Any]:
        """Run ABMIL over patch features.

        Args:
            features: (N, FEATURE_DIM) per-patch feature array.
            target_mutations: genes to predict.

        Returns dict:
            {
              "predictions": [ {gene, mutated, probability, confidence}, ... ],
              "attention": np.ndarray (N,)  # per-patch attention weights
            }

        Raises RuntimeError if torch unavailable so main.py can degrade gracefully.
        """
        if not self.available:
            raise RuntimeError("torch unavailable; cannot run MIL aggregation.")
        if features is None or len(features) == 0:
            raise RuntimeError("No patch features provided to MIL aggregator.")

        feats_t = torch.from_numpy(np.asarray(features, dtype=np.float32)).to(self.device)
        with torch.no_grad():
            base_probs, attention = self.model(feats_t)
        base_probs = base_probs.cpu().numpy()
        attention = attention.cpu().numpy()

        predictions = []
        for gene in target_mutations:
            p = self._gene_prob(base_probs, gene)
            predictions.append(
                {
                    "gene": gene,
                    "mutated": bool(p >= 0.5),
                    "probability": round(p, 4),
                    "confidence": confidence_from_prob(p),
                }
            )

        return {"predictions": predictions, "attention": attention}


def confidence_from_prob(prob: float) -> str:
    """Confidence bucketing per the service spec.

    "high" if |prob - 0.5| > 0.35, "medium" if > 0.15, else "low".
    """
    margin = abs(prob - 0.5)
    if margin > 0.35:
        return "high"
    if margin > 0.15:
        return "medium"
    return "low"


# Module-level singleton (matches the gnn_service convention).
mil_service = MILService()
