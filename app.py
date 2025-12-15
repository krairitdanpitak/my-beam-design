import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import numpy as np

# ==========================================
# 1. DATABASE & HELPER FUNCTIONS
# ==========================================

# Rebar Database (Based on GAS Code Source 10-11)
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
    """Format number with commas and decimals"""
    if n is None or pd.isna(n): return "-"
    return f"{n:,.{digits}f}"


def beta1FromFc(fc_MPa):
    """ACI 318-19 Beta1 (Source 26-28)"""
    if fc_MPa <= 28: return 0.85
    b1 = 0.85 - 0.05 * ((fc_MPa - 28) / 7)
    return max(0.65, b1)


def phiFlexureFromStrain(eps_t):
    """ACI 318-19 Phi from Strain (Source 28-30)"""
    if eps_t <= 0.002: return 0.65
    if eps_t >= 0.005: return 0.90
    return 0.65 + (eps_t - 0.002) * (0.25 / 0.003)


def flexureSectionResponse(As_mm2, fc, fy, bw, d, Es=200000, eps_cu=0.003):
    """
    Iteration for strain compatibility (Source 31-44)
    Returns: phiMn, phi, eps_t, etc.
    """
    beta1 = beta1FromFc(fc)
    fs = fy
    # Initial guess
    a = (As_mm2 * fs) / (0.85 * fc * bw) if fc > 0 else 0
    c = a / beta1 if beta1 > 0 else 0

    # Iterate to find equilibrium
    for i in range(30):
        if c <= 0: break
        eps_t = eps_cu * (d - c) / c
        fs_new = min(fy, Es * eps_t)
        # Avoid tension stiffening issues in iteration logic, keep simple as source
        fs_new = max(fs_new, -fy)  # Compression check? Source doesn't handle doubly explicitly here

        a_new = (As_mm2 * fs_new) / (0.85 * fc * bw)

        if abs(fs_new - fs) < 0.1 and abs(a_new - a) < 0.1:
            fs = fs_new
            a = a_new
            break
        fs = fs_new
        a = a_new
        c = a / beta1

    c = a / beta1
    eps_t = eps_cu * (d - c) / c if c > 0 else 999
    phi = phiFlexureFromStrain(eps_t)

    T = As_mm2 * fs
    Mn = T * (d - a / 2.0)
    phiMn = phi * Mn

    return {'beta1': beta1, 'a': a, 'c': c, 'eps_t': eps_t, 'fs': fs, 'phi': phi, 'Mn': Mn, 'phiMn': phiMn}


# ==========================================
# 2. MAIN CALCULATION LOGIC (PORTED)
# ==========================================
def process_calculation(inputs):
    calc_rows = []

    def sec(title): calc_rows.append(["SECTION", title, "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        calc_rows.append([item, formula, subs, result, unit, status])

    # 1. Parse Inputs (Source 3-9)
    b_cm = inputs['b']
    h_cm = inputs['h']
    cover_cm = inputs['cover']
    agg_mm = inputs.get('agg', 20)
    fc_ksc = inputs['fc']
    fy_ksc = inputs['fy']
    fyt_ksc = inputs['fyt']

    mainBarKey = inputs['mainBar']
    stirrupBarKey = inputs['stirrupBar']

    # Moments (tf-m)
    MuCases = [
        {'key': "L_TOP", 'title': "Left (Top) Mu(-)", 'Mu_tfm': inputs['mu_L_n']},
        {'key': "L_BOT", 'title': "Left (Bot) Mu(+)", 'Mu_tfm': inputs['mu_L_p']},
        {'key': "M_TOP", 'title': "Mid (Top) Mu(-)", 'Mu_tfm': inputs['mu_M_n']},
        {'key': "M_BOT", 'title': "Mid (Bot) Mu(+)", 'Mu_tfm': inputs['mu_M_p']},
        {'key': "R_TOP", 'title': "Right (Top) Mu(-)", 'Mu_tfm': inputs['mu_R_n']},
        {'key': "R_BOT", 'title': "Right (Bot) Mu(+)", 'Mu_tfm': inputs['mu_R_p']},
    ]
    # Shear (tf)
    VuCases = [
        {'key': "V_L", 'title': "Left", 'Vu_tf': inputs['vu_L']},
        {'key': "V_M", 'title': "Mid", 'Vu_tf': inputs['vu_M']},
        {'key': "V_R", 'title': "Right", 'Vu_tf': inputs['vu_R']},
    ]

    # 3. Unit Conversions (Source 13-19)
    ksc_to_MPa = 0.0980665
    fc = fc_ksc * ksc_to_MPa
    fy = fy_ksc * ksc_to_MPa
    fyt = fyt_ksc * ksc_to_MPa

    bw = b_cm * 10  # mm
    h = h_cm * 10  # mm
    cover = cover_cm * 10  # mm

    db_st = BAR_INFO[stirrupBarKey]['d_mm']
    db_main = BAR_INFO[mainBarKey]['d_mm']
    d = h - cover - db_st - db_main / 2.0

    # 6. Header Rows (Source 53-57)
    sec("1. MATERIAL & SECTION PARAMETERS")
    beta1 = beta1FromFc(fc)

    # As Min (Source 45-46)
    rho1 = 0.25 * math.sqrt(fc) / fy
    rho2 = 1.4 / fy
    As_min = max(rho1, rho2) * bw * d

    # As Max (Source 47-49)
    Es = 200000
    eps_cu = 0.003
    eps_y = fy / Es
    rho_bal = 0.85 * beta1 * (fc / fy) *