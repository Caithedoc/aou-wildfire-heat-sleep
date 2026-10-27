"""Palette learned from Lu et al. 2025 Nat Health (s44360-025-00011-y), sampled from the PDF on 2026-09-28.
Fig 1 bubble map: viridis.  Fig 2 choropleth: RdYlGn.  Fig 3 correlation: RdBu.
Fig 4 categorical regions: matplotlib named colours purple / gold / blue / red / orange / green.
Fig 5 fits: pure red line with pink (red, alpha .2) 95% band over grey points."""
import matplotlib.pyplot as plt
LU_PURPLE, LU_GOLD, LU_BLUE, LU_RED, LU_ORANGE, LU_GREEN = "#800080", "#FFD700", "#0000FF", "#FF0000", "#FFA500", "#008000"
LU_CAT6 = [LU_PURPLE, LU_GOLD, LU_BLUE, LU_RED, LU_ORANGE, LU_GREEN]
SEQ = "viridis"; DIV = "RdBu"; CHORO = "RdYlGn"
def seq(n, lo=0.08, hi=0.88):
    """n ordered colours from viridis, avoiding the pale-yellow end for lines on white."""
    return [plt.cm.viridis(lo + (hi - lo) * i / max(n - 1, 1)) for i in range(n)]
CAT = {"hot_only": LU_ORANGE, "smoke_only": LU_BLUE, "both": LU_RED}   # exposure categories used across figures
ACC = LU_RED; BAND_ALPHA = 0.18
