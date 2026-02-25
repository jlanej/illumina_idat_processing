#!/usr/bin/env python3
"""
ml_cnv_calling.py

Experimental deep learning-based CNV caller using a 1D-CNN + Bi-LSTM
sequence-to-sequence model (CNVSegmenter) to predict copy number states
from BAF and LRR intensities produced by the Illumina IDAT processing
pipeline.

This approach is designed as a SOTA replacement for traditional HMM-based
CNV callers (PennCNV, QuantiSNP). It uses:
  - 1D CNN layers for local smoothing and feature extraction of noisy
    LRR/BAF signals
  - Bidirectional LSTM to capture long-range sequence context and identify
    copy number transition boundaries
  - An explicit inter-probe distance channel (log-scaled) to handle
    variable probe densities across different Illumina array designs

Input channels (per probe):
  1. LRR  — Log R Ratio
  2. BAF  — B Allele Frequency
  3. Distance — log10(distance to next probe + 1)

Output: 5-class per-probe prediction (CN states 0, 1, 2, 3, 4)

Training labels come from truth-set BED files such as the 1000 Genomes
Project structural variant call sets.

Usage:
    # Train
    python3 ml_cnv_calling.py train \\
        --bcf stage2_reclustered.bcf \\
        --truth-bed 1000g_cnv_truthset.bed \\
        --sample NA12878 \\
        --output-model cnv_model.pt

    # Predict
    python3 ml_cnv_calling.py predict \\
        --bcf stage2_reclustered.bcf \\
        --model cnv_model.pt \\
        --sample NA12878 \\
        --output-bed predicted_cnvs.bed

Requires: torch, pysam, numpy, pandas, scikit-learn
Install:  pip install -r requirements_ml.txt
"""

import argparse
import sys

import numpy as np
import pandas as pd
import pysam
import torch
import torch.nn as nn
from sklearn.utils.class_weight import compute_class_weight


# =============================================================================
# Model Definition
# =============================================================================

class CNVSegmenter(nn.Module):
    """
    Sequence-to-sequence model for per-probe copy number prediction.

    Architecture:
        Input (3 channels: LRR, BAF, Distance)
          -> 1D CNN stack (local feature extraction / denoising)
          -> Bidirectional LSTM (sequence context)
          -> Fully connected output (5-class CN prediction per probe)

    The 1D CNN smooths noisy array signals and extracts local patterns,
    while the BiLSTM captures long-range dependencies and transition
    boundaries that are critical for accurate segmentation.
    """

    def __init__(self, input_channels=3, cnn_channels=64, lstm_hidden=128,
                 lstm_layers=2, num_classes=5, dropout=0.3):
        super().__init__()

        # 1D CNN for local feature extraction
        self.cnn = nn.Sequential(
            nn.Conv1d(input_channels, cnn_channels, kernel_size=7, padding=3),
            nn.BatchNorm1d(cnn_channels),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Conv1d(cnn_channels, cnn_channels, kernel_size=5, padding=2),
            nn.BatchNorm1d(cnn_channels),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Conv1d(cnn_channels, cnn_channels, kernel_size=3, padding=1),
            nn.BatchNorm1d(cnn_channels),
            nn.ReLU(),
        )

        # Bidirectional LSTM for sequence context
        self.lstm = nn.LSTM(
            input_size=cnn_channels,
            hidden_size=lstm_hidden,
            num_layers=lstm_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if lstm_layers > 1 else 0.0,
        )

        # Fully connected output layer (BiLSTM outputs 2 * hidden_size)
        self.fc = nn.Sequential(
            nn.Linear(lstm_hidden * 2, lstm_hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(lstm_hidden, num_classes),
        )

    def forward(self, x):
        """
        Args:
            x: Tensor of shape (batch, channels=3, seq_len)

        Returns:
            Tensor of shape (batch, seq_len, num_classes) — logits
        """
        # CNN: (batch, channels, seq_len) -> (batch, cnn_channels, seq_len)
        cnn_out = self.cnn(x)

        # Transpose for LSTM: (batch, seq_len, cnn_channels)
        lstm_in = cnn_out.permute(0, 2, 1)

        # LSTM: (batch, seq_len, 2 * lstm_hidden)
        lstm_out, _ = self.lstm(lstm_in)

        # FC per time step: (batch, seq_len, num_classes)
        logits = self.fc(lstm_out)

        return logits


# =============================================================================
# Data Loading
# =============================================================================

def load_bcf_data(bcf_path, sample_id):
    """
    Read LRR, BAF, chromosome, and position from a BCF file for one sample.

    Uses pysam to read FORMAT/LRR and FORMAT/BAF arrays from the pipeline's
    output .bcf files.

    Args:
        bcf_path: Path to .bcf file produced by the pipeline.
        sample_id: Sample identifier in the BCF.

    Returns:
        DataFrame with columns: chrom, pos, lrr, baf
    """
    vcf = pysam.VariantFile(bcf_path)

    if sample_id not in vcf.header.samples:
        available = list(vcf.header.samples)
        raise ValueError(
            f"Sample '{sample_id}' not found in BCF. "
            f"Available samples: {available[:10]}{'...' if len(available) > 10 else ''}"
        )

    records = []
    for rec in vcf.fetch():
        sample_data = rec.samples[sample_id]
        lrr = sample_data.get("LRR", None)
        baf = sample_data.get("BAF", None)
        if lrr is None or baf is None:
            continue
        # Handle tuple values from pysam
        if isinstance(lrr, tuple):
            lrr = lrr[0]
        if isinstance(baf, tuple):
            baf = baf[0]
        if lrr is None or baf is None:
            continue
        records.append({
            "chrom": rec.chrom,
            "pos": rec.pos,
            "lrr": float(lrr),
            "baf": float(baf),
        })

    vcf.close()

    if not records:
        raise ValueError(f"No records with LRR/BAF found for sample '{sample_id}'")

    df = pd.DataFrame(records)
    df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)
    return df


def compute_distance_channel(df):
    """
    Compute log-scaled inter-probe distance channel.

    For probes on the same chromosome, distance = log10(pos_next - pos_current + 1).
    At chromosome boundaries, distance is set to a large sentinel value (6.0,
    corresponding to ~1 Mb) to signal a break in genomic continuity.

    This channel makes the model robust to different Illumina array probe
    densities (e.g., GSA ~700K probes vs Omni2.5 ~2.5M probes).

    Args:
        df: DataFrame with chrom and pos columns.

    Returns:
        numpy array of log-scaled distances, length = len(df).
    """
    positions = df["pos"].values
    chroms = df["chrom"].values
    n = len(df)
    distances = np.full(n, 6.0, dtype=np.float32)  # sentinel for last probe

    for i in range(n - 1):
        if chroms[i] == chroms[i + 1]:
            dist = positions[i + 1] - positions[i]
            distances[i] = np.log10(max(dist, 1) + 1)
        # else: keep sentinel value (chromosome boundary)

    return distances


def load_truth_bed(bed_path):
    """
    Load a truth-set BED file defining known CNV regions.

    Expected BED format (tab-separated):
        chrom  start  end  cn_state

    where cn_state is the integer copy number (0, 1, 3, or 4).
    Regions not covered by the BED file are assumed to be CN=2 (diploid).

    Compatible with 1000 Genomes Project structural variant truth sets
    and other standard CNV benchmark resources.

    Args:
        bed_path: Path to the truth-set .bed file.

    Returns:
        DataFrame with columns: chrom, start, end, cn
    """
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        names=["chrom", "start", "end", "cn"],
        comment="#",
        dtype={"chrom": str, "start": int, "end": int, "cn": int},
    )
    return df


def assign_labels(probe_df, truth_df):
    """
    Assign CN labels to each probe based on truth-set BED regions.

    Probes falling within a truth-set region get that region's CN state.
    All other probes default to CN=2 (normal diploid).

    Args:
        probe_df: DataFrame with chrom, pos columns.
        truth_df: DataFrame with chrom, start, end, cn columns.

    Returns:
        numpy array of integer labels (0–4), length = len(probe_df).
    """
    labels = np.full(len(probe_df), 2, dtype=np.int64)  # default CN=2

    for _, region in truth_df.iterrows():
        mask = (
            (probe_df["chrom"] == region["chrom"])
            & (probe_df["pos"] >= region["start"])
            & (probe_df["pos"] <= region["end"])
        )
        labels[mask.values] = region["cn"]

    return labels


def create_windows(lrr, baf, distances, labels, window_size, step_size):
    """
    Create sliding windows of input features and labels for training.

    Args:
        lrr: 1D array of LRR values.
        baf: 1D array of BAF values.
        distances: 1D array of log-scaled inter-probe distances.
        labels: 1D array of CN labels (or None for inference).
        window_size: Number of probes per window.
        step_size: Stride between consecutive windows.

    Returns:
        Tuple of (features, labels) numpy arrays.
        features shape: (num_windows, 3, window_size)
        labels shape: (num_windows, window_size) or None
    """
    n = len(lrr)
    if n < window_size:
        # Pad to window_size
        pad_len = window_size - n
        lrr = np.concatenate([lrr, np.zeros(pad_len)])
        baf = np.concatenate([baf, np.full(pad_len, 0.5)])
        distances = np.concatenate([distances, np.full(pad_len, 6.0)])
        if labels is not None:
            labels = np.concatenate([labels, np.full(pad_len, 2, dtype=np.int64)])
        n = window_size

    windows_x = []
    windows_y = []

    for start in range(0, n - window_size + 1, step_size):
        end = start + window_size
        window = np.stack([
            lrr[start:end],
            baf[start:end],
            distances[start:end],
        ], axis=0)  # shape: (3, window_size)
        windows_x.append(window)
        if labels is not None:
            windows_y.append(labels[start:end])

    # Include the last window if it wasn't covered
    if n > window_size and (n - window_size) % step_size != 0:
        start = n - window_size
        window = np.stack([
            lrr[start:],
            baf[start:],
            distances[start:],
        ], axis=0)
        windows_x.append(window)
        if labels is not None:
            windows_y.append(labels[start:])

    features = np.array(windows_x, dtype=np.float32)
    if labels is not None:
        label_arr = np.array(windows_y, dtype=np.int64)
        return features, label_arr
    return features, None


# =============================================================================
# Training
# =============================================================================

def train_model(args):
    """Train the CNVSegmenter model."""
    print(f"Loading BCF data from: {args.bcf}", file=sys.stderr)
    probe_df = load_bcf_data(args.bcf, args.sample)

    print(f"  Loaded {len(probe_df)} probes for sample {args.sample}", file=sys.stderr)

    print(f"Loading truth-set BED: {args.truth_bed}", file=sys.stderr)
    truth_df = load_truth_bed(args.truth_bed)
    print(f"  Loaded {len(truth_df)} truth regions", file=sys.stderr)

    # Compute features
    lrr = probe_df["lrr"].values.astype(np.float32)
    baf = probe_df["baf"].values.astype(np.float32)
    distances = compute_distance_channel(probe_df)
    labels = assign_labels(probe_df, truth_df)

    print(f"  Label distribution: {dict(zip(*np.unique(labels, return_counts=True)))}", file=sys.stderr)

    # Create sliding windows
    window_size = args.window_size
    step_size = args.step_size
    features, label_windows = create_windows(
        lrr, baf, distances, labels, window_size, step_size
    )
    print(f"  Created {len(features)} training windows "
          f"(window={window_size}, step={step_size})", file=sys.stderr)

    # Compute class weights for weighted cross-entropy (handle CN=2 imbalance)
    all_labels = label_windows.flatten()
    classes_present = np.unique(all_labels)
    weights = compute_class_weight("balanced", classes=classes_present, y=all_labels)
    class_weights = np.ones(5, dtype=np.float32)
    for cls, w in zip(classes_present, weights):
        class_weights[cls] = w

    print(f"  Class weights: {dict(enumerate(class_weights))}", file=sys.stderr)

    # Set up device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  Using device: {device}", file=sys.stderr)

    # Build model
    model = CNVSegmenter(
        input_channels=3,
        cnn_channels=args.cnn_channels,
        lstm_hidden=args.lstm_hidden,
        lstm_layers=args.lstm_layers,
        num_classes=5,
        dropout=args.dropout,
    ).to(device)

    weight_tensor = torch.tensor(class_weights, dtype=torch.float32).to(device)
    criterion = nn.CrossEntropyLoss(weight=weight_tensor)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    # Convert to tensors
    X = torch.tensor(features, dtype=torch.float32)
    Y = torch.tensor(label_windows, dtype=torch.long)
    dataset = torch.utils.data.TensorDataset(X, Y)
    loader = torch.utils.data.DataLoader(
        dataset, batch_size=args.batch_size, shuffle=True
    )

    # Training loop
    print(f"\nTraining for {args.epochs} epochs...", file=sys.stderr)
    model.train()
    for epoch in range(args.epochs):
        epoch_loss = 0.0
        num_batches = 0
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)

            logits = model(batch_x)  # (batch, seq_len, 5)
            # Reshape for cross-entropy: (batch * seq_len, 5) vs (batch * seq_len,)
            loss = criterion(
                logits.reshape(-1, 5),
                batch_y.reshape(-1),
            )

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()
            num_batches += 1

        avg_loss = epoch_loss / max(num_batches, 1)
        if (epoch + 1) % max(1, args.epochs // 10) == 0 or epoch == 0:
            print(f"  Epoch {epoch + 1}/{args.epochs}  Loss: {avg_loss:.4f}", file=sys.stderr)

    # Save model
    torch.save(model.state_dict(), args.output_model)
    print(f"\nModel saved to: {args.output_model}", file=sys.stderr)


# =============================================================================
# Prediction
# =============================================================================

def predict_cnv(args):
    """Run CNV prediction using a trained model."""
    print(f"Loading BCF data from: {args.bcf}", file=sys.stderr)
    probe_df = load_bcf_data(args.bcf, args.sample)
    print(f"  Loaded {len(probe_df)} probes for sample {args.sample}", file=sys.stderr)

    # Compute features
    lrr = probe_df["lrr"].values.astype(np.float32)
    baf = probe_df["baf"].values.astype(np.float32)
    distances = compute_distance_channel(probe_df)

    # Set up device and model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  Using device: {device}", file=sys.stderr)

    model = CNVSegmenter(
        input_channels=3,
        cnn_channels=args.cnn_channels,
        lstm_hidden=args.lstm_hidden,
        lstm_layers=args.lstm_layers,
        num_classes=5,
        dropout=0.0,  # No dropout at inference
    ).to(device)

    model.load_state_dict(torch.load(args.model, map_location=device, weights_only=True))
    model.eval()

    # Sliding window inference
    window_size = args.window_size
    step_size = args.step_size
    features, _ = create_windows(lrr, baf, distances, None, window_size, step_size)

    n_probes = len(probe_df)
    vote_counts = np.zeros((n_probes, 5), dtype=np.float32)

    print(f"  Running inference on {len(features)} windows...", file=sys.stderr)

    with torch.no_grad():
        for i in range(0, len(features), args.batch_size):
            batch = torch.tensor(
                features[i:i + args.batch_size], dtype=torch.float32
            ).to(device)
            logits = model(batch)  # (batch, window_size, 5)
            probs = torch.softmax(logits, dim=-1).cpu().numpy()

            for j in range(len(probs)):
                win_idx = i + j
                start = win_idx * step_size
                end = min(start + window_size, n_probes)
                actual_len = end - start
                vote_counts[start:end] += probs[j, :actual_len, :]

    # Pick the class with the highest accumulated probability
    predictions = np.argmax(vote_counts, axis=1)

    # Collapse adjacent probes with the same non-diploid CN into BED regions
    cnv_regions = collapse_predictions(probe_df, predictions)

    # Write output BED
    with open(args.output_bed, "w") as f:
        f.write("#chrom\tstart\tend\tcn_state\tnum_probes\n")
        for region in cnv_regions:
            f.write(f"{region['chrom']}\t{region['start']}\t{region['end']}\t"
                    f"{region['cn']}\t{region['n_probes']}\n")

    print(f"\nPredicted {len(cnv_regions)} CNV regions -> {args.output_bed}", file=sys.stderr)
    print(f"  CN distribution: {dict(zip(*np.unique(predictions, return_counts=True)))}", file=sys.stderr)


def collapse_predictions(probe_df, predictions):
    """
    Collapse adjacent probes with the same predicted non-diploid CN state
    into contiguous BED regions.

    Only reports regions where CN != 2 (i.e., deletions CN<2 or gains CN>2).

    Args:
        probe_df: DataFrame with chrom and pos columns.
        predictions: Array of predicted CN states per probe.

    Returns:
        List of dicts with chrom, start, end, cn, n_probes.
    """
    regions = []
    n = len(predictions)
    i = 0

    while i < n:
        cn = int(predictions[i])
        if cn == 2:
            i += 1
            continue

        # Start a new CNV region
        chrom = probe_df.iloc[i]["chrom"]
        start_pos = int(probe_df.iloc[i]["pos"])
        end_pos = start_pos
        n_probes = 1
        j = i + 1

        while j < n:
            if (int(predictions[j]) == cn
                    and probe_df.iloc[j]["chrom"] == chrom):
                end_pos = int(probe_df.iloc[j]["pos"])
                n_probes += 1
                j += 1
            else:
                break

        regions.append({
            "chrom": chrom,
            "start": start_pos,
            "end": end_pos,
            "cn": cn,
            "n_probes": n_probes,
        })
        i = j

    return regions


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Experimental deep learning CNV caller (1D-CNN + BiLSTM)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Train on 1000 Genomes truth set
  python3 ml_cnv_calling.py train \\
      --bcf stage2_reclustered.bcf \\
      --truth-bed 1000g_cnv_truthset.bed \\
      --sample NA12878 \\
      --output-model cnv_model.pt

  # Predict CNVs
  python3 ml_cnv_calling.py predict \\
      --bcf stage2_reclustered.bcf \\
      --model cnv_model.pt \\
      --sample NA12878 \\
      --output-bed predicted_cnvs.bed

Truth-set BED format (tab-separated, no header):
  chrom  start  end  cn_state

  The highly curated 1000 Genomes structural variant truth sets
  (e.g., from the 1000 Genomes Phase 3 integrated SV call set or
  the Chaisson et al. multi-platform validated set) can be converted
  to this format for training labels.
""",
    )

    # Common model architecture arguments
    arch_args = argparse.ArgumentParser(add_help=False)
    arch_args.add_argument("--cnn-channels", type=int, default=64,
                           help="CNN feature channels (default: 64)")
    arch_args.add_argument("--lstm-hidden", type=int, default=128,
                           help="LSTM hidden size (default: 128)")
    arch_args.add_argument("--lstm-layers", type=int, default=2,
                           help="Number of LSTM layers (default: 2)")
    arch_args.add_argument("--window-size", type=int, default=512,
                           help="Sliding window size in probes (default: 512)")
    arch_args.add_argument("--step-size", type=int, default=256,
                           help="Sliding window step size (default: 256)")
    arch_args.add_argument("--batch-size", type=int, default=32,
                           help="Batch size (default: 32)")

    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    # --- train subcommand ---
    train_parser = subparsers.add_parser(
        "train",
        parents=[arch_args],
        help="Train the CNV model on labeled data",
        description="Train the CNVSegmenter model using a BCF file and truth-set BED.",
    )
    train_parser.add_argument("--bcf", required=True,
                              help="Input BCF file with LRR/BAF FORMAT fields")
    train_parser.add_argument("--truth-bed", required=True,
                              help="Truth-set BED file (chrom, start, end, cn_state)")
    train_parser.add_argument("--sample", required=True,
                              help="Sample ID to extract from BCF")
    train_parser.add_argument("--output-model", required=True,
                              help="Output path for trained model weights (.pt)")
    train_parser.add_argument("--epochs", type=int, default=50,
                              help="Number of training epochs (default: 50)")
    train_parser.add_argument("--lr", type=float, default=1e-3,
                              help="Learning rate (default: 1e-3)")
    train_parser.add_argument("--dropout", type=float, default=0.3,
                              help="Dropout rate (default: 0.3)")

    # --- predict subcommand ---
    predict_parser = subparsers.add_parser(
        "predict",
        parents=[arch_args],
        help="Predict CNVs using a trained model",
        description="Run sliding window CNV inference on a BCF file.",
    )
    predict_parser.add_argument("--bcf", required=True,
                                help="Input BCF file with LRR/BAF FORMAT fields")
    predict_parser.add_argument("--model", required=True,
                                help="Trained model weights file (.pt)")
    predict_parser.add_argument("--sample", required=True,
                                help="Sample ID to extract from BCF")
    predict_parser.add_argument("--output-bed", required=True,
                                help="Output BED file for predicted CNV regions")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "train":
        train_model(args)
    elif args.command == "predict":
        predict_cnv(args)


if __name__ == "__main__":
    main()
