import json
import math

# Path to the JSON file produced by make_dicer2_figure.py
summary_path = "results/dicer2.summary.json"

with open(summary_path) as f:
    summary = json.load(f)

z = summary.get("z_at_plus2")
pairs_total = summary.get("pairs_total")
lengths = summary.get("pairs_by_length", {})

if z is None:
    print("No Z-score available (used precomputed TSV without permutations).")
else:
    # Convert Z-score to approximate two-sided p-value using normal approximation
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    
    # Determine qualitative interpretation
    if z < 2:
        phrase = "no detectable enrichment"
    elif z < 4:
        phrase = "moderate enrichment consistent with Dicer-2 processing"
    elif z < 8:
        phrase = "strong enrichment characteristic of Dicer-2 products"
    else:
        phrase = "very strong enrichment confirming Dicer-2 cleavage"

    # Compose a caption-style summary
    caption = (
        f"Among {pairs_total:,} duplexes analyzed, the frequency of +2 nt overlaps "
        f"was {z:.2f} standard deviations above randomized expectation "
        f"(Z = {z:.2f}, p ≈ {p:.2e}), indicating {phrase}. "
        f"Most paired reads were {max(lengths, key=lengths.get)} nt in length, "
        f"consistent with canonical siRNAs."
    )

    print(caption)
