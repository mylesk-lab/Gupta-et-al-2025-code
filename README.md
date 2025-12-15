# Gupta-et-al-2025-code
Python Scripts for 2-nt Overhang Analysis and figure generation

Usage: python3 make_dicer2_figure.py [Name].txt (txt file with siRNAs in HiSeq format) --out-prefix results/dicer2
Writes:
  results/dicer2.panelA_length_dist.svg/png
  results/dicer2.panelB_offset_hist.svg/png
  results/dicer2.panelC_heatmap.svg/png
  results/dicer2.panelD_density.svg/png
  results/dicer2.offset_hist.tsv
  results/dicer2.length_offset_heatmap.tsv
  results/dicer2.summary.json

Usage: python3 Summary_json.py
Outputs:
Among [X number] of duplexes analyzed, the frequency of +2 nt overlaps was [X number] standard deviations..... (Z = [X], p ≈ [X]), indicating [X]. Most paired reads were [X} nt in length.......
