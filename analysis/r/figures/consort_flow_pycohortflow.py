#!/usr/bin/env python3
"""
V28 CONSORT Participant Flow Diagram (Fig 1A) — 大字版
=====================================================
所有字体增大至15pt左右，确保PDF中清晰可读。
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from pathlib import Path

OUTPUT_DIR = Path("/data/小米云盘/文档/论文/血磷/TTE/v28/output/figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# 1. 样式
# ============================================================================
WHITE  = "#FFFFFF"
BLACK  = "#000000"
GRAY50 = "#808080"
FONT   = "Liberation Sans"

# 方框参数 (大字版)
BOX_W = 0.32
BOX_H = 0.090
EXCL_W = 0.28
EXCL_H = 0.48
ARM_W = 0.22
ARM_H = 0.075
OUT_W = 0.22
OUT_H = 0.070

# 垂直位置
Y_ICU     = 0.90
Y_ELIG    = 0.72
Y_EXPAND  = 0.54
Y_ARMS    = 0.36
Y_OUTCOME = 0.20

# 水平位置
X_CENTER = 0.30
X_EXCL   = 0.72
X_LEFT   = 0.13
X_RIGHT  = 0.47

CORNER = 0.008

# 排除详情 (全部不加粗，紧凑排版)
EXCL_DETAIL_LINES = [
    (False, "Excluded (n = 44,210)", None),
    (False, "Demographic exclusion:", "23,000"),
    (False, "  Non-first ICU admission", "12,524"),
    (False, "  Age < 18 years", "5,468"),
    (False, "  ICU stay < 24 hours", "5,008"),
    (False, "Phosphate-related exclusion:", "6,611"),
    (False, "  (No phosphate / > 0.80 mmol/L)", None),
    (False, "Washout < 48 h:", "11,637"),
    (False, "Medical exclusion criteria:", "2,962"),
    (False, "  (8 criteria, study_cohort)", None),
    (False, "  eGFR < 15", "1,374"),
    (False, "  Short stay (elective surgery)", "796"),
    (False, "  Maintenance dialysis", "269"),
    (False, "  Malignancy hypercalcemia", "215"),
    (False, "  Hyperparathyroidism", "212"),
    (False, "  Pregnancy", "92"),
    (False, "  Phosphate binders > 30 d", "4"),
]

# ============================================================================
# 2. 绘图函数
# ============================================================================
def draw_box(ax, x, y, w, h, text, fs=15, lw=1.5, bold=False):
    """白底黑字方框"""
    box = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={CORNER}",
        facecolor=WHITE, edgecolor=BLACK, linewidth=lw, zorder=2
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=fs, color=BLACK, fontfamily=FONT,
            fontweight="bold" if bold else "normal", zorder=3)

def draw_arrow(ax, x0, y0, x1, y1, lw=1.2, color=GRAY50):
    """箭头"""
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", lw=lw, color=color),
                zorder=1)

def draw_excl_box(ax, x, y, w, h, lines):
    """排除方框 — 浅灰底 + 实线黑框"""
    box = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={CORNER}",
        facecolor="#F5F5F5", edgecolor=BLACK, linewidth=1.2, zorder=2
    )
    ax.add_patch(box)
    line_h = 0.026
    start_y = y + h/2 - 0.032
    left_x = x - w/2 + 0.02
    right_x = x + w/2 - 0.02
    for i, (is_bold, label, num) in enumerate(lines):
        if not label:
            continue
        ax.text(left_x, start_y - i * line_h, label, ha="left", va="top",
                fontsize=13, color=BLACK, fontfamily=FONT,
                fontweight="normal", zorder=3)
        if num is not None:
            ax.text(right_x, start_y - i * line_h, num, ha="right", va="top",
                    fontsize=13, color=BLACK, fontfamily=FONT, zorder=3)

def draw_title(ax):
    """标题"""
    ax.text(0.02, 0.98, "Figure 1A. Participant Flow Diagram",
            ha="left", va="top", fontsize=22, fontweight="bold",
            color=BLACK, fontfamily=FONT, transform=ax.transAxes)

# ============================================================================
# 3. 主绘图
# ============================================================================
def draw_flowchart(ax):
    """完整流程图"""
    # ---- Step 1: ICU ----
    draw_box(ax, X_CENTER, Y_ICU, BOX_W, BOX_H,
             "ICU Admissions\nAssessed for Eligibility\nn = 76,540",
             fs=16, lw=1.8, bold=True)

    # ---- Arrow: ICU -> Eligible ----
    draw_arrow(ax, X_CENTER, Y_ICU - BOX_H/2, X_CENTER, Y_ELIG + BOX_H/2)

    # ---- Step 2: Eligible ----
    draw_box(ax, X_CENTER, Y_ELIG, BOX_W, BOX_H,
             "Inclusion: Phosphate ≤ 0.80 mmol/L\nn = 46,929 eligible",
             fs=16, lw=1.8, bold=True)

    # ---- Exclusion Box ----
    excl_y = (Y_ICU - BOX_H/2 + Y_ELIG + BOX_H/2) / 2
    draw_excl_box(ax, X_EXCL, excl_y, EXCL_W, EXCL_H, EXCL_DETAIL_LINES)

    # ---- Exclusion Arrow ----
    ax.annotate("", xy=(X_EXCL - EXCL_W/2, excl_y), xytext=(X_CENTER, excl_y),
                arrowprops=dict(arrowstyle="-|>", lw=1.0, color=GRAY50), zorder=1)
    # Junction dot
    ax.add_patch(Circle((X_CENTER, excl_y), radius=0.006,
                         facecolor=GRAY50, edgecolor=WHITE, linewidth=0.5, zorder=4))

    # ---- Arrow: Eligible -> Expanded ----
    draw_arrow(ax, X_CENTER, Y_ELIG - BOX_H/2, X_CENTER, Y_EXPAND + BOX_H/2)

    # ---- Step 3: Expanded ----
    draw_box(ax, X_CENTER, Y_EXPAND, BOX_W, BOX_H,
             "Sequential Trial Expansion\n"
             "24 nested trials\n"
             "32,330 trial-periods\n"
             "(22,685 unique patients)",
             fs=15, lw=1.8, bold=True)

    # ---- Arrow: Expanded -> Arms ----
    draw_arrow(ax, X_CENTER, Y_EXPAND - BOX_H/2, X_CENTER, Y_ARMS + ARM_H/2)

    # ---- Split junction ----
    ax.add_patch(Circle((X_CENTER, Y_ARMS + ARM_H/2), radius=0.007,
                         facecolor=BLACK, edgecolor=WHITE, linewidth=0.5, zorder=4))

    # ---- Horizontal split lines ----
    ax.annotate("", xy=(X_LEFT, Y_ARMS + ARM_H/2),
                xytext=(X_CENTER, Y_ARMS + ARM_H/2),
                arrowprops=dict(arrowstyle="-", lw=1.2, color=BLACK), zorder=1)
    ax.annotate("", xy=(X_RIGHT, Y_ARMS + ARM_H/2),
                xytext=(X_CENTER, Y_ARMS + ARM_H/2),
                arrowprops=dict(arrowstyle="-", lw=1.2, color=BLACK), zorder=1)

    # ---- Step 4: Treatment ----
    draw_box(ax, X_LEFT, Y_ARMS, ARM_W, ARM_H,
             "IV Phosphate\nSupplementation\nn = 5,147 (15.9%)",
             fs=15, lw=1.5)

    # ---- Step 5: Control ----
    draw_box(ax, X_RIGHT, Y_ARMS, ARM_W, ARM_H,
             "Conservative\nManagement\nn = 27,183 (84.1%)",
             fs=15, lw=1.5)

    # ---- Arrows: Arms -> Outcomes ----
    draw_arrow(ax, X_LEFT, Y_ARMS - ARM_H/2, X_LEFT, Y_OUTCOME + OUT_H/2)
    draw_arrow(ax, X_RIGHT, Y_ARMS - ARM_H/2, X_RIGHT, Y_OUTCOME + OUT_H/2)

    # ---- Step 6: Death ----
    draw_box(ax, X_LEFT, Y_OUTCOME, OUT_W, OUT_H,
             "28-Day Death\nn = 5,610 (17.4%)",
             fs=15, lw=1.5)

    # ---- Step 7: Alive ----
    draw_box(ax, X_RIGHT, Y_OUTCOME, OUT_W, OUT_H,
             "Alive / Censored\nat Day 28\nn = 26,720 (82.6%)",
             fs=15, lw=1.5)

# ============================================================================
# 4. 导出
# ============================================================================
def export_all():
    print("===== 绘制流程图 (大字版) =====")
    print(f"  ICU: 76,540")
    print(f"  Eligible: 46,929")
    print(f"  Expanded: 32,330")
    print(f"  TX: 5,147 + CTRL: 27,183 = 32,330")
    print(f"  Death: 5,610 + Alive: 26,720 = 32,330")
    print(f"  Excluded: 44,210 = 23,000 + 6,611 + 11,637 + 2,962")

    for fmt, ext, dpi in [("png", "png", 600),
                          ("pdf", "pdf", 300),
                          ("tiff", "tiff", 300)]:
        fig, ax = plt.subplots(figsize=(22, 17))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.1)
        ax.axis("off")
        fig.patch.set_facecolor(WHITE)
        ax.set_facecolor(WHITE)
        draw_title(ax)
        draw_flowchart(ax)

        fpath = OUTPUT_DIR / f"Fig1A_Flow_Diagram.{ext}"
        fig.savefig(fpath, dpi=dpi, bbox_inches="tight",
                    format=ext, facecolor=WHITE, edgecolor="none")
        plt.close(fig)
        print(f"  {fpath.name} ({fpath.stat().st_size / 1024:.0f} KB)")

    print("===== 完成 =====")

if __name__ == "__main__":
    export_all()