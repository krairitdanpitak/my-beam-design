import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import numpy as np

# ==========================================
# 1. DATABASE & HELPER FUNCTIONS
# ==========================================

# Rebar Database
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
    """Format number with commas and decimals (Safe version)"""
    try:
        if n is None: return "-"
        val = float(n)
        if math.isnan(val): return "-"
        return f"{val:,.{digits}f}"
    except:
        return "-"


def beta1FromFc(fc_MPa):
    """ACI 318-19 Beta1"""
    if fc_MPa <= 28: return 0.85
    b1 = 0.85 - 0.05 * ((fc_MPa - 28) / 7)
    return max(0.65, b1)


def phiFlexureFromStrain(eps_t):
    """ACI 318-19 Phi from Strain"""
    if eps_t <= 0.002: return 0.65
    if eps_t >= 0.005: return 0.90
    return 0.65 + (eps_t - 0.002) * (0.25 / 0.003)


def flexureSectionResponse(As_mm2, fc, fy, bw, d, Es=200000, eps_cu=0.003):
    """
    Iteration for strain compatibility
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
        # Simple convergence check
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

    def sec(title):
        calc_rows.append(["SECTION", title, "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        calc_rows.append([item, formula, subs, result, unit, status])

    # 1. Parse Inputs
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

    # 3. Unit Conversions
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

    # 6. Header Rows
    sec("1. MATERIAL & SECTION PARAMETERS")
    beta1 = beta1FromFc(fc)

    # As Min
    rho1 = 0.25 * math.sqrt(fc) / fy
    rho2 = 1.4 / fy
    As_min = max(rho1, rho2) * bw * d

    # As Max
    Es = 200000
    eps_cu = 0.003
    eps_y = fy / Es
    rho_bal = 0.85 * beta1 * (fc / fy) * (eps_cu / (eps_cu + eps_y))
    As_max = 0.75 * rho_bal * bw * d

    row("Materials", "-", f"fc'={fmt(fc, 2)} MPa, fy={fmt(fy, 2)} MPa", "-", "-")
    row("Section", "-", f"{fmt(bw, 0)} x {fmt(h, 0)} mm", "-", "mm")
    row("Effective d", "d=h-cv-dst-db/2", f"{h}-{cover}-{db_st}-{db_main / 2}", f"{fmt(d, 1)}", "mm")
    row("As,min", "max(0.25√fc'/fy, 1.4/fy)bd", f"max({rho1:.5f}, {rho2:.5f})*{bw}*{d}", f"{fmt(As_min, 0)}", "mm²")
    row("As,max", "0.75*rho_bal*bd", f"rho_bal={rho_bal:.4f}", f"{fmt(As_max, 0)}", "mm²")

    # 7. FLEXURE DESIGN
    sec("2. FLEXURE DESIGN (STRAIN-BASED φ)")

    bar_counts = {}
    flex_all_ok = True

    for case in MuCases:
        title = case['title']
        Mu_tfm = case['Mu_tfm']

        if Mu_tfm == 0:
            bar_counts[case['key']] = 2  # Minimum bars
            continue

        Mu_Nmm = Mu_tfm * 9806650.0

        # Binary Search for As Required
        As_lo = As_min
        As_hi = As_lo

        # Expand Hi
        for _ in range(30):
            r = flexureSectionResponse(As_hi, fc, fy, bw, d)
            if r['phiMn'] >= Mu_Nmm: break
            As_hi *= 1.3
            if As_hi > As_max:
                As_hi = As_max
                break

        # Binary Search
        As_req = As_hi
        for _ in range(50):
            As_mid = 0.5 * (As_lo + As_hi)
            r = flexureSectionResponse(As_mid, fc, fy, bw, d)
            if r['phiMn'] >= Mu_Nmm:
                As_hi = As_mid
            else:
                As_lo = As_mid
        As_req = As_hi

        # Provide Bars
        bar_area = BAR_INFO[mainBarKey]['A_cm2'] * 100  # mm2
        n = math.ceil(As_req / bar_area)
        if n < 2: n = 2
        As_prov = n * bar_area

        # Check Provided
        rProv = flexureSectionResponse(As_prov, fc, fy, bw, d)
        passStrength = rProv['phiMn'] >= Mu_Nmm
        passAsMax = As_req <= As_max + 1e-6

        # Clear Spacing Check
        usable = bw - 2.0 * (cover + db_st)
        clear = (usable - n * db_main) / (n - 1) if n > 1 else usable - db_main
        clear_req = max(db_main, 25.0, 4.0 * agg_mm / 3.0)
        passClear = clear >= clear_req - 1

        overallPass = passStrength and passAsMax and passClear
        if not overallPass: flex_all_ok = False

        bar_counts[case['key']] = n

        # Rows
        row(f"{title}: Mu", "-", "-", f"{fmt(Mu_tfm, 3)}", "tf-m")
        row(f"{title}: As,req", "Solve φMn ≥ Mu", f"As_min={fmt(As_min, 0)}", f"{fmt(As_req, 0)}", "mm²")
        row(f"{title}: Provide", f"Use {mainBarKey}", f"{n}-{mainBarKey} (As={fmt(As_prov, 0)})",
            "OK" if overallPass else "NO", "-")
        row(f"{title}: Check", "φMn ≥ Mu", f"φMn={fmt(rProv['phiMn'] / 9.8e6, 3)}", "PASS" if passStrength else "FAIL",
            "tf-m")
        if n > 1:
            row(f"{title}: Spacing", f"clr ≥ {fmt(clear_req, 1)}", f"clr={fmt(clear, 1)}",
                "PASS" if passClear else "FAIL", "mm")

    # 8. SHEAR DESIGN
    sec("3. SHEAR DESIGN")

    Vc_N = 0.17 * math.sqrt(fc) * bw * d
    phi_v = 0.75
    phiVc_N = phi_v * Vc_N
    Vn_max_N = 0.66 * math.sqrt(fc) * bw * d
    phiVnmax_N = phi_v * Vn_max_N

    row("Vc", "0.17√fc' bd", f"0.17√{fmt(fc, 1)}*{bw}*{fmt(d, 1)}", f"{fmt(Vc_N / 9806.65, 3)}", "tf")
    row("φVc", "0.75 * Vc", "-", f"{fmt(phiVc_N / 9806.65, 3)}", "tf")

    Av = 2.0 * BAR_INFO[stirrupBarKey]['A_cm2'] * 100  # mm2 (2 legs)

    # Av/s Min Limit
    req1 = 0.062 * math.sqrt(fc) * bw / fyt
    req2 = 0.35 * bw / fyt
    s_avmin = Av / max(req1, req2)