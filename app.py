import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import numpy as np
import pandas as pd

# ==========================================
# 1. DATABASE & CONFIG
# ==========================================
st.set_page_config(page_title="RC Beam Designer Pro", layout="wide")

# Rebar Database (Based on Source)
BAR_INFO = {
    'RB6': {'A_cm2': 0.283, 'd_mm': 6},
    'RB9': {'A_cm2': 0.636, 'd_mm': 9},
    'DB10': {'A_cm2': 0.785, 'd_mm': 10},
    'DB12': {'A_cm2': 1.131, 'd_mm': 12},
    'DB16': {'A_cm2': 2.011, 'd_mm': 16},
    'DB20': {'A_cm2': 3.142, 'd_mm': 20},
    'DB25': {'A_cm2': 4.909, 'd_mm': 25},
    'DB28': {'A_cm2': 6.158, 'd_mm': 28}
}


def fmt(n, digits=3):
    """Format helper (Source)"""
    try:
        if n is None: return "-"
        val = float(n)
        if math.isnan(val): return "-"
        return f"{val:,.{digits}f}"
    except:
        return "-"


# ==========================================
# 2. CALCULATION ENGINE (ACI 318-19)
# ==========================================
def beta1FromFc(fc_MPa):
    """(Source)"""
    if fc_MPa <= 28: return 0.85
    b1 = 0.85 - 0.05 * ((fc_MPa - 28) / 7)
    return max(0.65, b1)


def phiFlexureFromStrain(eps_t):
    """(Source)"""
    if eps_t <= 0.002: return 0.65
    if eps_t >= 0.005: return 0.90
    return 0.65 + (eps_t - 0.002) * (0.25 / 0.003)


def flexureSectionResponse(As_mm2, fc, fy, bw, d, Es=200000, eps_cu=0.003):
    """Iteration for equilibrium (Source)"""
    beta1 = beta1FromFc(fc)
    fs = fy

    # Initial guess
    a = (As_mm2 * fs) / (0.85 * fc * bw) if fc > 0 else 0
    c = a / beta1 if beta1 > 0 else 0

    # Iterate
    for i in range(50):
        if c <= 0.1: c = 0.1  # Prevent div/0

        eps_t = eps_cu * (d - c) / c
        fs_new = min(fy, Es * eps_t)
        fs_new = max(fs_new, -fy)

        a_new = (As_mm2 * fs_new) / (0.85 * fc * bw)

        if abs(fs_new - fs) < 0.1 and abs(a_new - a) < 0.1:
            fs = fs_new
            a = a_new
            break
        fs = fs_new
        a = a_new
        c = a / beta1

    c = a / beta1
    eps_t = eps_cu * (d - c) / c if c > 0 else 0.005
    phi = phiFlexureFromStrain(eps_t)

    T = As_mm2 * fs
    Mn = T * (d - a / 2.0)
    phiMn = phi * Mn

    return {'beta1': beta1, 'a': a, 'c': c, 'eps_t': eps_t, 'fs': fs, 'phi': phi, 'Mn': Mn, 'phiMn': phiMn}


def process_calculation(inputs):
    """Main Logic (Source)"""
    calc_rows = []

    def sec(title):
        calc_rows.append(["SECTION", title, "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        calc_rows.append([item, formula, subs, result, unit, status])

    # 1. Parse Inputs
    b_cm = inputs['b']
    h_cm = inputs['h']
    cover_cm = inputs['cover']
    agg_mm = inputs.get('agg', 20)

    # Convert ksc -> MPa (Source)
    ksc_to_MPa = 0.0980665
    fc = inputs['fc'] * ksc_to_MPa
    fy = inputs['fy'] * ksc_to_MPa
    fyt = inputs['fyt'] * ksc_to_MPa

    bw = b_cm * 10
    h = h_cm * 10
    cover = cover_cm * 10

    barKey = inputs['mainBar']
    stirKey = inputs['stirrupBar']

    db_st = BAR_INFO[stirKey]['d_mm']
    db_main = BAR_INFO[barKey]['d_mm']
    d = h - cover - db_st - db_main / 2.0

    # 2. Header Rows
    sec("1. MATERIAL & SECTION PARAMETERS")
    beta1 = beta1FromFc(fc)

    # As Min (Source)
    rho1 = 0.25 * math.sqrt(fc) / fy
    rho2 = 1.4 / fy
    As_min = max(rho1, rho2) * bw * d

    # As Max (Source)
    Es = 200000
    eps_cu = 0.003
    eps_y = fy / Es
    rho_bal = 0.85 * beta1 * (fc / fy) * (eps_cu / (eps_cu + eps_y))
    As_max = 0.75 * rho_bal * bw * d

    row("Materials", "-", f"fc'={fmt(fc, 2)} MPa, fy={fmt(fy, 2)} MPa", "-", "-")
    row("Section", "-", f"{fmt(bw, 0)} x {fmt(h, 0)} mm", "-", "mm")
    row("As,min", "max(0.25√fc'/fy, 1.4/fy)bd", "-", f"{fmt(As_min, 0)}", "mm²")

    # 3. FLEXURE
    sec("2. FLEXURE DESIGN")

    MuCases = [
        {'key': "L_TOP", 'title': "Left (Top) Mu(-)", 'Mu_tfm': inputs['mu_L_n']},
        {'key': "L_BOT", 'title': "Left (Bot) Mu(+)", 'Mu_tfm': inputs['mu_L_p']},
        {'key': "M_TOP", 'title': "Mid (Top) Mu(-)", 'Mu_tfm': inputs['mu_M_n']},
        {'key': "M_BOT", 'title': "Mid (Bot) Mu(+)", 'Mu_tfm': inputs['mu_M_p']},
        {'key': "R_TOP", 'title': "Right (Top) Mu(-)", 'Mu_tfm': inputs['mu_R_n']},
        {'key': "R_BOT", 'title': "Right (Bot) Mu(+)", 'Mu_tfm': inputs['mu_R_p']},
    ]

    bar_counts = {}
    flex_ok = True

    for case in MuCases:
        title = case['title']
        Mu_tfm = case['Mu_tfm']
        key = case['key']

        # --- FIX: Handle Mu=0 correctly (Source) ---
        if Mu_tfm <= 0.001:
            bar_counts[key] = 2  # Min bars
            # Skip row output but KEEP data in bar_counts
            continue

        Mu_Nmm = Mu_tfm * 9806650.0

        # Binary Search for As (Source)
        As_lo = As_min
        As_hi = As_lo

        # Expand
        for _ in