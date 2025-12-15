#!/usr/bin/env python3
import sys
import argparse
import math
import json
from collections import defaultdict, Counter
import numpy as np
import matplotlib.pyplot as plt
import pathlib

def parse_hiseq_line(line):
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    parts = line.split("\t")
    if len(parts) < 5:
        return None
    try:
        strand = parts[1]
        virus  = parts[2]
        start  = int(parts[3])
        seq    = parts[4]
        if strand not in {"+","-"}:
            return None
        return (virus, strand, start, seq)
    except Exception:
        return None

def load_reads(paths, min_len=18, max_len=30):
    plus_reads = defaultdict(list)   # (virus,L) -> list[(start, seq)]
    minus_reads_5p = defaultdict(lambda: defaultdict(list))  # (virus,L) -> start -> seqs
    minus_reads_LM = defaultdict(lambda: defaultdict(list))  # (virus,L) -> start -> seqs
    total, kept = 0, 0
    for path in paths:
        fh = sys.stdin if path == "-" else open(path, "r")
        with fh:
            for line in fh:
                total += 1
                rec = parse_hiseq_line(line)
                if rec is None:
                    continue
                virus, strand, start, seq = rec
                L = len(seq)
                if L < min_len or L > max_len:
                    continue
                kept += 1
                key = (virus, L)
                if strand == "+":
                    plus_reads[key].append((start, seq))
                else:
                    minus_reads_5p[key][start].append(seq)  # assume recorded start == minus 5' end (rightmost coord)
                    minus_reads_LM[key][start].append(seq)  # assume recorded start == leftmost coord
    return plus_reads, minus_reads_5p, minus_reads_LM, total, kept

def match_pairs(plus_reads, minus_5p, minus_LM):
    """
    Return matched pairs and basic tallies.
    Each match has standardized offset = (minus_left_edge - plus_left_edge).
    Under 5' convention: minus_left_edge = s_minus_recorded - (L-1)
    Under leftmost convention: minus_left_edge = s_minus_recorded
    Expected Dicer-2 signature: standardized offset == +2.
    """
    matches = []
    len_counts = Counter()
    virus_counts = Counter()
    total_plus = sum(len(v) for v in plus_reads.values())
    total_minus = 0
    for d in minus_5p.values():
        total_minus += sum(len(v) for v in d.values())

    for (virus, L), plist in plus_reads.items():
        mp5 = minus_5p.get((virus, L), {})
        mlm = minus_LM.get((virus, L), {})
        for s_plus, _seqp in plist:
            # convention A: minus start recorded is its 5' end (rightmost coord)
            e_plus = s_plus + L - 1
            s_minus_5p = e_plus + 2  # recorded start expected if perfect 2-nt 3' overhangs
            for i, seqm in enumerate(mp5.get(s_minus_5p, [])):
                minus_left = s_minus_5p - (L - 1)
                offset = minus_left - s_plus  # standardize to left-edge difference
                matches.append((virus, L, "minus_start_is_5prime", s_plus, s_minus_5p, offset))
                len_counts[L] += 1
                virus_counts[virus] += 1
            # convention B: minus start recorded is leftmost coord
            s_minus_lm = s_plus + 2
            for i, seqm in enumerate(mlm.get(s_minus_lm, [])):
                minus_left = s_minus_lm
                offset = minus_left - s_plus
                matches.append((virus, L, "minus_start_is_leftmost", s_plus, s_minus_lm, offset))
                len_counts[L] += 1
                virus_counts[virus] += 1

    return matches, {"total_plus": total_plus, "total_minus": total_minus,
                     "pairs_by_length": dict(len_counts), "pairs_by_virus": dict(virus_counts)}

def compute_histograms(matches, offset_min=-10, offset_max=10):
    """
    Build histograms:
      - overall offset histogram
      - per-length offset histogram (heatmap data)
    """
    offsets = [off for (_v, L, _c, _sp, _sm, off) in matches if offset_min <= off <= offset_max]
    overall = Counter(offsets)
    per_len = defaultdict(Counter)
    for v, L, c, sp, sm, off in matches:
        if offset_min <= off <= offset_max:
            per_len[L][off] += 1
    return overall, per_len

def permutation_enrichment(plus_reads, minus_5p, minus_LM, Ls, num_perm=200, rng_seed=1,
                           offset_min=-10, offset_max=10):
    """
    Null model: shuffle minus recorded starts within each (virus, L) bucket, preserving multiplicities.
    Recompute standardized offsets relative to plus starts, tally histograms.
    """
    rng = np.random.default_rng(rng_seed)
    hist_perms = np.zeros((num_perm, offset_max - offset_min + 1), dtype=np.float64)
    # Pre-extract structures for performance
    buckets = {}
    for (virus, L) in Ls:
        plist = plus_reads.get((virus, L), [])
        right_dict = minus_5p.get((virus, L), {})
        left_dict  = minus_LM.get((virus, L), {})
        # Expand minus positions with multiplicity (for shuffling)
        m5_list = []
        for pos, seqs in right_dict.items():
            m5_list.extend([pos] * len(seqs))
        ml_list = []
        for pos, seqs in left_dict.items():
            ml_list.extend([pos] * len(seqs))
        buckets[(virus, L)] = (plist, m5_list, ml_list)

    offsets_range = list(range(offset_min, offset_max+1))

    for p in range(num_perm):
        counts = Counter()
        for (virus, L), (plist, m5_list, ml_list) in buckets.items():
            if not plist or (not m5_list and not ml_list):
                continue
            # shuffle copies
            m5 = list(m5_list)
            ml = list(ml_list)
            if m5:
                rng.shuffle(m5)
            if ml:
                rng.shuffle(ml)
            i5 = 0
            il = 0
            # Pair each plus read with one minus "candidate" under each convention if available,
            # but *randomized* minus positions.
            for s_plus, _ in plist:
                # 5' convention random pick
                if i5 < len(m5):
                    s_minus_rec = m5[i5]
                    i5 += 1
                    minus_left = s_minus_rec - (L - 1)
                    off = minus_left - s_plus
                    if offset_min <= off <= offset_max:
                        counts[off] += 1
                # leftmost convention random pick
                if il < len(ml):
                    s_minus_lm = ml[il]
                    il += 1
                    minus_left = s_minus_lm
                    off = minus_left - s_plus
                    if offset_min <= off <= offset_max:
                        counts[off] += 1

        hist_perms[p, :] = [counts[o] for o in offsets_range]

    perm_mean = np.mean(hist_perms, axis=0)
    perm_std  = np.std(hist_perms, axis=0, ddof=1)  # sample std
    return offsets_range, perm_mean, perm_std

def binned_density(matches, bin_size=1000):
    """ Return per-virus binned counts of duplexes based on plus start coordinate. """
    bins = defaultdict(Counter)  # virus -> bin_index -> count
    for virus, L, c, s_plus, s_minus_rec, off in matches:
        bi = s_plus // bin_size
        bins[virus][bi] += 1
    return bins

def make_figure(out_prefix, overall_hist, per_len_hist, offsets_range, perm_mean, perm_std,
                len_counts, density_bins, split_panels=False):
    # Sort lengths for consistent heatmap order
    lengths = sorted(per_len_hist.keys())
    offset_list = offsets_range

    # Prepare heatmap matrix (len x offsets)
    if lengths:
        heat = np.zeros((len(lengths), len(offset_list)), dtype=float)
        for i, L in enumerate(lengths):
            for j, off in enumerate(offset_list):
                heat[i, j] = per_len_hist[L].get(off, 0.0)

    # Panel A: length distribution
    def plot_panel_A():
        plt.figure()
        xs = sorted(len_counts.keys())
        ys = [len_counts[x] for x in xs]
        plt.bar(xs, ys)
        plt.xlabel("siRNA length (nt)")
        plt.ylabel("Duplex count")
        plt.title("Matched duplex length distribution")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}.panelA_length_dist.svg")
        plt.savefig(f"{out_prefix}.panelA_length_dist.png", dpi=300)
        plt.close()

    # Panel B: offset histogram with 2-nt Z-score in title
    def plot_panel_B():
        plt.figure()
        ys_obs = [overall_hist.get(o, 0) for o in offset_list]
        plt.bar(offset_list, ys_obs)
        # Compute Z for offset==2 if available
        try:
            idx2 = offset_list.index(2)
            mu2 = perm_mean[idx2] if perm_mean is not None else float('nan')
            sd2 = perm_std[idx2] if perm_std is not None else float('nan')
            z2 = (ys_obs[idx2] - mu2) / sd2 if (sd2 and not math.isnan(sd2)) else float('nan')
            title = f"Offset histogram (peak at +2); Z@+2 = {z2:.2f}"
        except ValueError:
            title = "Offset histogram (peak at +2)"
        plt.xlabel("Minus left edge − Plus left edge (nt)")
        plt.ylabel("Duplex count")
        plt.title(title)
        plt.tight_layout()
        plt.savefig(f"{out_prefix}.panelB_offset_hist.svg")
        plt.savefig(f"{out_prefix}.panelB_offset_hist.png", dpi=300)
        plt.close()

    # Panel C: heatmap length vs offset
    def plot_panel_C():
        if not lengths:
            return
        plt.figure()
        plt.imshow(heat, aspect='auto', origin='lower',
                   extent=[offset_list[0]-0.5, offset_list[-1]+0.5, lengths[0]-0.5, lengths[-1]+0.5])
        plt.colorbar(label="Duplex count")
        plt.xlabel("Offset (minus left − plus left)")
        plt.ylabel("Length (nt)")
        plt.title("Length vs offset")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}.panelC_heatmap.svg")
        plt.savefig(f"{out_prefix}.panelC_heatmap.png", dpi=300)
        plt.close()

    # Panel D: genome density by 1kb bins, one plot per virus (stacked vertically)
    def plot_panel_D():
        if not density_bins:
            return
        viruses = sorted(density_bins.keys())
        # If lots of viruses, just plot the first 6
        max_to_plot = min(6, len(viruses))
        plt.figure(figsize=(8, 2*max_to_plot))
        for i, virus in enumerate(viruses[:max_to_plot], start=1):
            ax = plt.subplot(max_to_plot, 1, i)
            bkeys = sorted(density_bins[virus].keys())
            xs = [k for k in bkeys]
            ys = [density_bins[virus][k] for k in bkeys]
            ax.plot(xs, ys)
            ax.set_ylabel("Pairs")
            ax.set_title(virus)
            if i == max_to_plot:
                ax.set_xlabel("Genome bin (1kb)")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}.panelD_density.svg")
        plt.savefig(f"{out_prefix}.panelD_density.png", dpi=300)
        plt.close()

    # Generate panels
    plot_panel_A()
    plot_panel_B()
    plot_panel_C()
    plot_panel_D()

def write_tables(out_prefix, overall_hist, per_len_hist, offsets_range, perm_mean, perm_std):
    # Overall histogram
    with open(f"{out_prefix}.offset_hist.tsv","w") as f:
        f.write("offset\tcount\tperm_mean\tperm_std\n")
        for i, off in enumerate(offsets_range):
            f.write(f"{off}\t{overall_hist.get(off,0)}\t")
            if perm_mean is None:
                f.write("NA\tNA\n")
            else:
                f.write(f"{perm_mean[i]:.6f}\t{perm_std[i]:.6f}\n")

    # Per-length histogram
    with open(f"{out_prefix}.length_offset_heatmap.tsv","w") as f:
        f.write("length\toffset\tcount\n")
        for L in sorted(per_len_hist.keys()):
            for off, c in sorted(per_len_hist[L].items()):
                f.write(f"{L}\t{off}\t{c}\n")

def main():
    ap = argparse.ArgumentParser(description="Make publication figure showing Dicer-2 2-nt overhang enrichment.")
    ap.add_argument("inputs", nargs="*", help="HISEQ-style mapped read files (use '-' for stdin)")
    ap.add_argument("--pairs-tsv", help="Optional: precomputed pairs TSV (virus length conv s_plus e_plus s_minus_recorded pair_index)")
    ap.add_argument("--min-len", type=int, default=18)
    ap.add_argument("--max-len", type=int, default=30)
    ap.add_argument("--offset-min", type=int, default=-10)
    ap.add_argument("--offset-max", type=int, default=10)
    ap.add_argument("--permutations", type=int, default=200, help="Number of permutations for null (increase for final figure).")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", default="dicer2_results/dicer2")
    args = ap.parse_args()

    pathlib.Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)

    if args.pairs_tsv:
        # Load matches directly
        matches = []
        with open(args.pairs_tsv, "r") as fh:
            header = fh.readline().strip().split("\t")
            idx = {h:i for i,h in enumerate(header)}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                virus = parts[idx["virus"]]
                L = int(parts[idx["length"]])
                conv = parts[idx["conv"]]
                s_plus = int(parts[idx["s_plus"]])
                s_minus_rec = int(parts[idx["s_minus_recorded"]])
                # Standardize offset = minus_left - plus_left
                if conv == "minus_start_is_5prime":
                    minus_left = s_minus_rec - (L - 1)
                else:
                    minus_left = s_minus_rec
                offset = minus_left - s_plus
                matches.append((virus, L, conv, s_plus, s_minus_rec, offset))
        plus_reads = minus_5p = minus_LM = None
        summary = {"pairs_by_length": dict(Counter([L for (_v,L,_c,_sp,_sm,_o) in matches]))}
    else:
        if not args.inputs:
            print("Error: provide HISEQ inputs or --pairs-tsv", file=sys.stderr)
            sys.exit(2)
        plus_reads, minus_5p, minus_LM, total, kept = load_reads(args.inputs, args.min_len, args.max_len)
        matches, summary = match_pairs(plus_reads, minus_5p, minus_LM)

    # Histograms
    overall_hist, per_len_hist = compute_histograms(matches, args.offset_min, args.offset_max)

    # Permutation null (only possible if we have raw reads)
    if args.pairs_tsv:
        offsets_range = list(range(args.offset_min, args.offset_max+1))
        perm_mean = perm_std = None
    else:
        Ls = list(set((v,L) for (v,L,_c,_sp,_sm,_o) in matches))
        offsets_range, perm_mean, perm_std = permutation_enrichment(
            plus_reads, minus_5p, minus_LM, Ls, num_perm=args.permutations, rng_seed=args.seed,
            offset_min=args.offset_min, offset_max=args.offset_max
        )

    # Density bins (1kb) from plus starts
    dens = binned_density(matches, bin_size=1000)

    # Write tables
    write_tables(args.out_prefix, overall_hist, per_len_hist, offsets_range, perm_mean, perm_std)

    # Make figures
    make_figure(args.out_prefix, overall_hist, per_len_hist, offsets_range, perm_mean, perm_std,
                summary.get("pairs_by_length", {}), dens, split_panels=False)

    # Write a compact JSON summary
    out_json = {
        "pairs_total": len(matches),
        "pairs_by_length": summary.get("pairs_by_length", {}),
        "z_at_plus2": None
    }
    if perm_mean is not None:
        try:
            i2 = offsets_range.index(2)
            obs2 = overall_hist.get(2, 0)
            mu = perm_mean[i2]
            sd = perm_std[i2]
            z = (obs2 - mu) / sd if sd else None
            out_json["z_at_plus2"] = z
        except Exception:
            out_json["z_at_plus2"] = None

    with open(f"{args.out_prefix}.summary.json","w") as f:
        json.dump(out_json, f, indent=2)

    # Friendly stdout note
    print("Wrote:")
    print(f"  {args.out_prefix}.panelA_length_dist.svg/png")
    print(f"  {args.out_prefix}.panelB_offset_hist.svg/png")
    print(f"  {args.out_prefix}.panelC_heatmap.svg/png")
    print(f"  {args.out_prefix}.panelD_density.svg/png")
    print(f"  {args.out_prefix}.offset_hist.tsv")
    print(f"  {args.out_prefix}.length_offset_heatmap.tsv")
    print(f"  {args.out_prefix}.summary.json")

if __name__ == "__main__":
    main()
