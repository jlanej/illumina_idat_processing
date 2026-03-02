#!/usr/bin/env python3
"""
ml_cnv_calling.py

Experimental Deep Learning-based CNV caller using a 1D CNN + Bi-LSTM model
(CNVSegmenter) to predict copy number states from BAF and LRR signals in the
pipeline's output BCF files.

Architecture:
    - 1D CNN layers: local smoothing and feature extraction of LRR/BAF noise
    - Bi-LSTM: captures sequence context and identifies transition boundaries
    - Fully connected output: 5 copy number classes (CN=0, 1, 2, 3, 4)

Input channels (3):
    - LRR  : Log R Ratio (normalized total intensity deviation)
    - BAF  : B Allele Frequency
    - Distance: log10-scaled genomic distance to the next probe
                (makes the model robust to variable Illumina array probe densities)

Usage:
    # Train
    python3 scripts/ml_cnv_calling.py train \\
        --bcf output/stage2/vcf/stage2_reclustered.bcf \\
        --truth truth_cnvs.bed \\
        --output model.pt

    # Predict
    python3 scripts/ml_cnv_calling.py predict \\
        --bcf output/stage2/vcf/stage2_reclustered.bcf \\
        --model model.pt \\
        --output predicted_cnvs.bed

Requires: torch, pysam, numpy, pandas, scikit-learn
"""

import argparse
import sys

import numpy as np

try:
    import pysam
except ImportError:
    pysam = None  # Checked at runtime

try:
    import pandas as pd
except ImportError:
    pd = None  # Checked at runtime

try:
    from sklearn.utils.class_weight import compute_class_weight
except ImportError:
    compute_class_weight = None  # Checked at runtime

try:
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset
    _TORCH_AVAILABLE = True
except ImportError:
    torch = None  # Checked at runtime
    _TORCH_AVAILABLE = False


# =============================================================================
# Model Definition
# =============================================================================

def _get_module_base():
    """Return nn.Module when torch is available, otherwise a plain object."""
    if _TORCH_AVAILABLE:
        import torch.nn as _nn
        return _nn.Module
    return object


class CNVSegmenter(_get_module_base()):
    """
    Sequence-to-sequence CNV segmentation model.

    Combines a 1D CNN (for local smoothing and feature extraction of LRR/BAF
    noise) and a Bi-LSTM (to capture long-range sequence context and identify
    copy number transition boundaries).

    Input:  (batch_size, 3, seq_len)  -- channels: LRR, BAF, Distance
    Output: (batch_size, seq_len, 5) -- logits for CN states 0, 1, 2, 3, 4
    """

    def __init__(self, cnn_channels=64, lstm_hidden=128, lstm_layers=2,
                 num_classes=5, dropout=0.3):
        super().__init__()

        # 1D CNN for local feature extraction and noise smoothing
        self.cnn = nn.Sequential(
            nn.Conv1d(in_channels=3, out_channels=cnn_channels,
                      kernel_size=7, padding=3),
            nn.BatchNorm1d(cnn_channels),
            nn.ReLU(),
            nn.Conv1d(in_channels=cnn_channels, out_channels=cnn_channels,
                      kernel_size=5, padding=2),
            nn.BatchNorm1d(cnn_channels),
            nn.ReLU(),
            nn.Dropout(p=dropout),
        )

        # Bi-LSTM for sequence-level context
        self.lstm = nn.LSTM(
            input_size=cnn_channels,
            hidden_size=lstm_hidden,
            num_layers=lstm_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if lstm_layers > 1 else 0.0,
        )

        # Fully connected output head
        self.fc = nn.Linear(lstm_hidden * 2, num_classes)

    def forward(self, x):
        """
        Args:
            x: Tensor of shape (batch_size, 3, seq_len)
        Returns:
            logits: Tensor of shape (batch_size, seq_len, num_classes)
        """
        # CNN: (batch, 3, seq_len) -> (batch, cnn_channels, seq_len)
        features = self.cnn(x)

        # Transpose for LSTM: (batch, seq_len, cnn_channels)
        features = features.transpose(1, 2)

        # Bi-LSTM: (batch, seq_len, lstm_hidden * 2)
        lstm_out, _ = self.lstm(features)

        # FC: (batch, seq_len, num_classes)
        logits = self.fc(lstm_out)
        return logits


# =============================================================================
# Data Loading from BCF
# =============================================================================

def load_bcf_data(bcf_path, sample=None):
    """
    Load BAF, LRR, and genomic positions from a BCF file.

    Reads FORMAT/BAF and FORMAT/LRR arrays from the pipeline's BCF output.
    Per-probe positions are used to compute the log10-scaled inter-probe
    distance (3rd input channel), which encodes array probe density.

    Args:
        bcf_path: Path to the BCF/VCF file.
        sample: Sample name to extract. If None, uses the first sample.

    Returns:
        dict with keys:
            'chroms'   : list of chromosome strings per probe
            'positions': np.ndarray of int, genomic position per probe
            'lrr'      : np.ndarray of float32, LRR values
            'baf'      : np.ndarray of float32, BAF values
            'distance' : np.ndarray of float32, log10(distance to next probe)
    """
    if pysam is None:
        raise ImportError("pysam is required: pip install pysam")

    chroms = []
    positions = []
    lrr_vals = []
    baf_vals = []

    with pysam.VariantFile(bcf_path) as vcf:
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError(f"No samples found in {bcf_path}")
        if sample is None:
            sample = samples[0]
        elif sample not in samples:
            raise ValueError(f"Sample '{sample}' not found in {bcf_path}. "
                             f"Available: {samples}")

        for rec in vcf.fetch():
            fmt = rec.samples[sample]
            lrr = fmt.get("LRR")
            baf = fmt.get("BAF")

            # Skip probes where both values are missing
            if lrr is None and baf is None:
                continue

            lrr_v = float(lrr) if lrr is not None and lrr != "." else float("nan")
            baf_v = float(baf) if baf is not None and baf != "." else float("nan")

            chroms.append(rec.chrom)
            positions.append(rec.pos)
            lrr_vals.append(lrr_v)
            baf_vals.append(baf_v)

    if not positions:
        raise ValueError(f"No LRR/BAF data found in {bcf_path} for sample {sample}. "
                         "Ensure FORMAT/LRR and FORMAT/BAF fields are present.")

    positions = np.array(positions, dtype=np.int64)
    lrr_arr = np.array(lrr_vals, dtype=np.float32)
    baf_arr = np.array(baf_vals, dtype=np.float32)

    # Compute log10-scaled distance to next probe (inter-probe distance channel)
    # For the last probe on each chromosome, use the median distance as fallback
    raw_dist = np.diff(positions, append=positions[-1] + 1)

    # Reset distance at chromosome boundaries
    for i in range(len(chroms) - 1):
        if chroms[i] != chroms[i + 1]:
            raw_dist[i] = 1  # Boundary: set to 1 (log10(1)=0)

    # Clamp to at least 1 bp and compute log10 distance
    raw_dist = np.maximum(raw_dist, 1).astype(np.float32)
    distance = np.log10(raw_dist)

    # Replace NaN LRR/BAF with neutral values (missing data)
    # LRR=0 represents diploid baseline (no copy number deviation)
    # BAF=0.5 represents the expected heterozygous state for diploid probes
    lrr_arr = np.nan_to_num(lrr_arr, nan=0.0)
    baf_arr = np.nan_to_num(baf_arr, nan=0.5)

    return {
        "chroms": chroms,
        "positions": positions,
        "lrr": lrr_arr,
        "baf": baf_arr,
        "distance": distance,
    }


def load_truth_labels(bed_path, chroms, positions):
    """
    Assign copy number labels to probes from a truth-set BED file.

    The BED file must have at least 4 columns:
        chrom, start (0-based), end, copy_number

    Probes not overlapping any BED interval are assigned CN=2 (diploid).

    Args:
        bed_path: Path to the truth BED file.
        chroms: List of chromosome strings per probe.
        positions: np.ndarray of 1-based genomic positions per probe.

    Returns:
        np.ndarray of int64, copy number labels (0-4) per probe.
    """
    if pd is None:
        raise ImportError("pandas is required for label loading: pip install pandas")

    bed = pd.read_csv(bed_path, sep="\t", header=None,
                      names=["chrom", "start", "end", "cn"],
                      usecols=[0, 1, 2, 3], comment="#",
                      dtype={"chrom": str, "start": int, "end": int, "cn": int})

    # Build interval lookup per chromosome
    intervals = {}
    for _, row in bed.iterrows():
        chrom = str(row["chrom"])
        cn = int(np.clip(row["cn"], 0, 4))
        intervals.setdefault(chrom, []).append((row["start"], row["end"], cn))

    labels = np.full(len(positions), fill_value=2, dtype=np.int64)

    for i, (chrom, pos) in enumerate(zip(chroms, positions)):
        # pos is 1-based (VCF/BCF), BED is 0-based half-open: [start, end)
        pos0 = pos - 1
        for start, end, cn in intervals.get(str(chrom), []):
            if start <= pos0 < end:
                labels[i] = cn
                break

    return labels


# =============================================================================
# Sliding Window Utilities
# =============================================================================

def make_windows(data, window_size, step):
    """
    Segment a 1D data array into overlapping windows.

    Args:
        data: np.ndarray of shape (n_probes, n_channels)
        window_size: Number of probes per window.
        step: Step size between consecutive windows.

    Returns:
        windows: np.ndarray of shape (n_windows, n_channels, window_size)
        starts:  np.ndarray of int, start index of each window in data.
    """
    n = len(data)
    starts = list(range(0, max(n - window_size + 1, 1), step))
    windows = []
    for s in starts:
        end = min(s + window_size, n)
        chunk = data[s:end]
        # Pad with zeros if shorter than window_size
        if len(chunk) < window_size:
            pad = np.zeros((window_size - len(chunk), data.shape[1]),
                           dtype=np.float32)
            chunk = np.concatenate([chunk, pad], axis=0)
        windows.append(chunk.T)  # (n_channels, window_size)
    return np.stack(windows, axis=0), np.array(starts, dtype=np.int64)


def aggregate_window_predictions(all_preds, starts, n_probes, window_size):
    """
    Aggregate per-probe class predictions from overlapping windows by majority vote.

    Args:
        all_preds: np.ndarray of shape (n_windows, window_size) -- predicted classes
        starts: np.ndarray of window start indices
        n_probes: Total number of probes
        window_size: Window size used during inference

    Returns:
        np.ndarray of shape (n_probes,) -- final predicted copy number per probe
    """
    # Accumulate votes
    votes = np.zeros((n_probes, 5), dtype=np.int32)
    for win_preds, start in zip(all_preds, starts):
        for offset, cls in enumerate(win_preds):
            idx = start + offset
            if idx < n_probes:
                votes[idx, int(cls)] += 1

    # Majority vote; default to CN=2 for probes with no coverage
    final = np.where(votes.sum(axis=1) > 0, votes.argmax(axis=1), 2)
    return final.astype(np.int64)


# =============================================================================
# BED Output
# =============================================================================

def predictions_to_bed(chroms, positions, predicted_cn, output_path):
    """
    Collapse adjacent probes with the same non-diploid predicted state into
    BED intervals and write to a BED file.

    Only segments with CN != 2 (i.e., gains CN>2 or losses CN<2) are output.

    Args:
        chroms: List of chromosome strings per probe.
        positions: np.ndarray of 1-based positions per probe.
        predicted_cn: np.ndarray of predicted copy numbers per probe.
        output_path: Path to write the output BED file.
    """
    n = len(chroms)
    records = []

    i = 0
    while i < n:
        cn = predicted_cn[i]
        if cn == 2:
            i += 1
            continue

        # Find the end of this run (same chrom and same CN)
        j = i + 1
        while j < n and predicted_cn[j] == cn and chroms[j] == chroms[i]:
            j += 1

        # BED: 0-based start, end is exclusive (use pos - 1 for start)
        start = int(positions[i]) - 1
        end = int(positions[j - 1])  # last probe position (1-based) is end
        records.append((chroms[i], start, end, int(cn)))
        i = j

    with open(output_path, "w") as f:
        f.write("# chrom\tstart\tend\tcopy_number\n")
        for chrom, start, end, cn in records:
            f.write(f"{chrom}\t{start}\t{end}\t{cn}\n")

    print(f"Wrote {len(records)} CNV segments to {output_path}", file=sys.stderr)


# =============================================================================
# Train subcommand
# =============================================================================

def cmd_train(args):
    """Train the CNVSegmenter on labelled BCF + truth BED data."""
    if torch is None:
        raise ImportError("torch is required: pip install torch")

    print(f"Loading BCF data from: {args.bcf}", file=sys.stderr)
    data = load_bcf_data(args.bcf, sample=args.sample)

    n_probes = len(data["positions"])
    print(f"  Loaded {n_probes} probes", file=sys.stderr)

    print(f"Loading truth labels from: {args.truth}", file=sys.stderr)
    labels = load_truth_labels(args.truth, data["chroms"], data["positions"])

    # Stack channels: (n_probes, 3)
    features = np.stack([data["lrr"], data["baf"], data["distance"]], axis=1)

    # Compute class weights to handle CN=2 imbalance
    if compute_class_weight is None:
        raise ImportError("scikit-learn is required: pip install scikit-learn")
    unique_classes = np.unique(labels)
    weights = compute_class_weight("balanced", classes=unique_classes, y=labels)
    class_weights = np.ones(5, dtype=np.float32)
    for cls, w in zip(unique_classes, weights):
        class_weights[cls] = float(w)

    print(f"  Class weights: {class_weights}", file=sys.stderr)

    # Build windows
    print("Building training windows...", file=sys.stderr)
    windows, starts = make_windows(features, args.window_size, args.step)
    win_labels = []
    for start in starts:
        end = min(start + args.window_size, n_probes)
        lbl = labels[start:end]
        if len(lbl) < args.window_size:
            pad = np.full(args.window_size - len(lbl), 2, dtype=np.int64)
            lbl = np.concatenate([lbl, pad])
        win_labels.append(lbl)
    win_labels = np.stack(win_labels, axis=0)

    print(f"  {len(windows)} windows of size {args.window_size}", file=sys.stderr)

    # Tensors
    X = torch.tensor(windows, dtype=torch.float32)
    Y = torch.tensor(win_labels, dtype=torch.long)

    dataset = TensorDataset(X, Y)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True,
                        drop_last=False)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  Training on device: {device}", file=sys.stderr)

    model = CNVSegmenter(
        cnn_channels=args.cnn_channels,
        lstm_hidden=args.lstm_hidden,
        lstm_layers=args.lstm_layers,
        dropout=args.dropout,
    ).to(device)

    weight_tensor = torch.tensor(class_weights, dtype=torch.float32).to(device)
    criterion = nn.CrossEntropyLoss(weight=weight_tensor)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        n_batches = 0
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            logits = model(xb)  # (batch, seq_len, 5)
            # CrossEntropyLoss expects (batch, classes, seq_len)
            loss = criterion(logits.permute(0, 2, 1), yb)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            total_loss += loss.item()
            n_batches += 1
        avg_loss = total_loss / max(n_batches, 1)
        print(f"  Epoch {epoch}/{args.epochs}  loss={avg_loss:.4f}", file=sys.stderr)

    torch.save(model.state_dict(), args.output)
    print(f"Model saved to: {args.output}", file=sys.stderr)


# =============================================================================
# Predict subcommand
# =============================================================================

def cmd_predict(args):
    """Run sliding window inference with a trained CNVSegmenter."""
    if torch is None:
        raise ImportError("torch is required: pip install torch")

    print(f"Loading BCF data from: {args.bcf}", file=sys.stderr)
    data = load_bcf_data(args.bcf, sample=args.sample)

    n_probes = len(data["positions"])
    print(f"  Loaded {n_probes} probes", file=sys.stderr)

    features = np.stack([data["lrr"], data["baf"], data["distance"]], axis=1)

    # Build sliding windows
    windows, starts = make_windows(features, args.window_size, args.step)
    print(f"  {len(windows)} inference windows of size {args.window_size}",
          file=sys.stderr)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = CNVSegmenter().to(device)
    state = torch.load(args.model, map_location=device)
    model.load_state_dict(state)
    model.eval()

    X = torch.tensor(windows, dtype=torch.float32)
    dataset = TensorDataset(X)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False)

    all_preds = []
    with torch.no_grad():
        for (xb,) in loader:
            xb = xb.to(device)
            logits = model(xb)  # (batch, seq_len, 5)
            preds = logits.argmax(dim=-1).cpu().numpy()  # (batch, seq_len)
            all_preds.append(preds)

    all_preds = np.concatenate(all_preds, axis=0)  # (n_windows, window_size)

    final_cn = aggregate_window_predictions(all_preds, starts, n_probes,
                                            args.window_size)

    predictions_to_bed(data["chroms"], data["positions"], final_cn, args.output)


# =============================================================================
# CLI
# =============================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        prog="ml_cnv_calling.py",
        description=(
            "Experimental Deep Learning CNV caller (1D CNN + Bi-LSTM). "
            "Predicts copy number states (0-4) from BAF and LRR signals "
            "in the pipeline's BCF output."
        ),
    )
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # ------------------------------------------------------------------
    # train subcommand
    # ------------------------------------------------------------------
    train_p = subparsers.add_parser(
        "train",
        help="Train CNVSegmenter on a labelled BCF + truth BED file.",
        description=(
            "Train the CNVSegmenter model. Reads LRR/BAF from a BCF file, "
            "assigns copy-number labels from a truth BED file, and trains "
            "using Weighted Cross-Entropy Loss to handle CN=2 class imbalance."
        ),
    )
    train_p.add_argument("--bcf", required=True,
                         help="Input BCF/VCF file with FORMAT/LRR and FORMAT/BAF fields")
    train_p.add_argument("--truth", required=True,
                         help="Truth-set BED file: chrom, start, end, copy_number")
    train_p.add_argument("--output", required=True,
                         help="Output path for trained model weights (.pt)")
    train_p.add_argument("--sample", default=None,
                         help="Sample name to use (default: first sample in BCF)")
    train_p.add_argument("--window-size", type=int, default=256,
                         help="Probes per training window (default: 256)")
    train_p.add_argument("--step", type=int, default=128,
                         help="Step size between windows (default: 128)")
    train_p.add_argument("--batch-size", type=int, default=32,
                         help="Training batch size (default: 32)")
    train_p.add_argument("--epochs", type=int, default=20,
                         help="Number of training epochs (default: 20)")
    train_p.add_argument("--lr", type=float, default=1e-3,
                         help="Learning rate (default: 0.001)")
    train_p.add_argument("--cnn-channels", type=int, default=64,
                         help="CNN output channels (default: 64)")
    train_p.add_argument("--lstm-hidden", type=int, default=128,
                         help="Bi-LSTM hidden size per direction (default: 128)")
    train_p.add_argument("--lstm-layers", type=int, default=2,
                         help="Number of Bi-LSTM layers (default: 2)")
    train_p.add_argument("--dropout", type=float, default=0.3,
                         help="Dropout rate (default: 0.3)")
    train_p.set_defaults(func=cmd_train)

    # ------------------------------------------------------------------
    # predict subcommand
    # ------------------------------------------------------------------
    pred_p = subparsers.add_parser(
        "predict",
        help="Run sliding window inference on a BCF file.",
        description=(
            "Predict CNV regions using a trained CNVSegmenter model. "
            "Applies sliding window inference over the BCF probes and "
            "outputs a BED file of predicted CNV regions (CN != 2)."
        ),
    )
    pred_p.add_argument("--bcf", required=True,
                        help="Input BCF/VCF file with FORMAT/LRR and FORMAT/BAF fields")
    pred_p.add_argument("--model", required=True,
                        help="Trained model weights file (.pt)")
    pred_p.add_argument("--output", required=True,
                        help="Output BED file for predicted CNV segments")
    pred_p.add_argument("--sample", default=None,
                        help="Sample name to use (default: first sample in BCF)")
    pred_p.add_argument("--window-size", type=int, default=256,
                        help="Probes per inference window (default: 256)")
    pred_p.add_argument("--step", type=int, default=64,
                        help="Step size between windows for sliding inference (default: 64)")
    pred_p.add_argument("--batch-size", type=int, default=64,
                        help="Inference batch size (default: 64)")
    pred_p.set_defaults(func=cmd_predict)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
