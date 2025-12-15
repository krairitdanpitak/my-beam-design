import streamlit as st
import matplotlib

# ใช้ Agg backend เพื่อป้องกัน Error เรื่อง GUI บน Server/Cloud
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import numpy as np

# ==========================================
# 1. SETUP & CONFIGURATION
# ==========================================
st.set_page_config(page_title="RC Beam Designer Pro", layout="wide")

# Database เหล็กเสริม (อ้างอิงพื้นที่หน้าตัดตามมาตรฐาน)
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
    """ฟังก์ชันจัดรูปแบบตัวเลขให้สวยงาม"""
    try:
        if n is None: return "-"
        val = float(n)
        if math.isnan(val): return "-"
        return f"{val:,.{digits}f}"
    except:
        return "-"


# ==========================================
# 2. CALCULATION LOGIC (ACI 318-19 Metric)
# ==========================================
def beta1FromFc(fc_MPa):
    """คำนวณค่า Beta1"""
    if fc_MPa <= 28: return 0.85
    b1 = 0.85 - 0.05 * ((fc_MPa - 28) / 7)
    return max(0.65, b1)


def phiFlexureFromStrain(eps_t):
    """คำนวณค่า Phi จาก Strain"""
    if eps_t <= 0.002: return 0.65
    if eps_t >= 0.005: return 0.90
    return 0.65 + (eps_t - 0.002) * (0.25 / 0.003)


def flexureSectionResponse(As_mm2, fc, fy, bw, d, Es=200000, eps_cu=0.003):
    """คำนวณกำลังรับโมเมนต์ดัด (Iterative Method)"""
    beta1 = beta1FromFc(fc)
    fs = fy

    # เริ่มต้นสมมติค่า a
    a = (As_mm2 * fs) / (0.85 * fc * bw) if fc > 0 else 0
    c = a / beta1 if beta1 > 0 else 0

    # วนลูปหาจุดสมดุลแรง (Equilibrium)
    for i in range(50):
        if c <= 0.1: c = 0.1
        eps_t = eps_cu * (d - c) / c

        # คำนวณ fs ใหม่ตาม Strain จริง
        fs_new = min(fy, Es * eps_t)
        fs_new = max(fs_new, -fy)

        a_new = (As_mm2 * fs_new) / (0.85 * fc * bw)

        # เช็คว่าลู่เข้าหรือยัง
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

    return {'phi': phi, 'phiMn': phiMn, 'eps_t': eps_t}


def solve_required_as(Mu_Nmm, As_min, As_max, fc, fy, bw, d):
    """Binary Search เพื่อหาปริมาณเหล็กที่ต้องการ (As Required)"""
    As_lo = As_min
    As_hi = As_lo

    # 1. ขยายขอบเขตบน (Expand Hi)
    for _ in range(30):
        r = flexureSectionResponse(As_hi, fc, fy, bw, d)
        if r['phiMn'] >= Mu_Nmm: break
        As_hi *= 1.3
        if As_hi > As_max:
            As_hi = As_max
            break

    # 2. ค้นหาค่าละเอียด (Binary Search)
    As_req = As_hi
    for _ in range(50):
        As_mid = 0.5 * (As_lo + As_hi)
        r = flexureSectionResponse(As_mid, fc, fy, bw, d)
        if r['phiMn'] >= Mu_Nmm:
            As_hi = As_mid
        else:
            As_lo = As_mid

    return As_hi


def process_calculation(inputs):
    """ฟังก์ชันหลักสำหรับคำนวณทั้งระบบ"""
    calc_rows = []

    # Helper ในการสร้างแถวตาราง
    def sec(title):
        calc_rows.append(["SECTION", title, "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        calc_rows.append([item, formula, subs, result, unit, status])

    # 1. ดึงค่า Input
    b_cm = inputs['b']
    h_cm = inputs['h']
    cover_cm = inputs['cover']
    agg_mm = inputs.get('agg', 20)

    # แปลงหน่วย ksc -> MPa
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

    # Effective Depth (d)
    d = h - cover - db_st - db_main / 2.0

    # 2. ส่วน Header ข้อมูลวัสดุ
    sec("1. MATERIAL & SECTION PARAMETERS")
    beta1 = beta1FromFc(fc)

    # As Min / As Max
    rho1 = 0.25 * math.sqrt(fc) / fy
    rho2 = 1.4 / fy
    As_min = max(rho1, rho2) * bw * d

    Es = 200000
    eps_cu = 0.003
    eps_y = fy / Es
    rho_bal = 0.85 * beta1 * (fc / fy) * (eps_cu / (eps_cu + eps_y))
    As_max = 0.75 * rho_bal * bw * d

    row("Materials", "-", f"fc'={fmt(fc, 2)} MPa", "-", "-")
    row("Section", "-", f"{fmt(bw, 0)} x {fmt(h, 0)} mm", "-", "mm")
    row("As,min", "max(0.25√fc'/fy, 1.4/fy)bd", "-", f"{fmt(As_min, 0)}", "mm²")

    # 3. ออกแบบรับแรงดัด (FLEXURE)
    sec("2. FLEXURE DESIGN")

    # รายการโมเมนต์ (เขียนแยกบรรทัดเพื่อป้องกัน Syntax Error)
