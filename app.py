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
    MuCases = [
        {'key': "L_TOP", 't': "Left (Top) Mu(-)", 'v': inputs['mu_L_n']},
        {'key': "L_BOT", 't': "Left (Bot) Mu(+)", 'v': inputs['mu_L_p']},
        {'key': "M_TOP", 't': "Mid (Top) Mu(-)", 'v': inputs['mu_M_n']},
        {'key': "M_BOT", 't': "Mid (Bot) Mu(+)", 'v': inputs['mu_M_p']},
        {'key': "R_TOP", 't': "Right (Top) Mu(-)", 'v': inputs['mu_R_n']},
        {'key': "R_BOT", 't': "Right (Bot) Mu(+)", 'v': inputs['mu_R_p']}
    ]

    bar_counts = {}
    flex_ok = True

    for case in MuCases:
        title = case['t']
        Mu_tfm = case['v']
        key = case['key']

        # กรณีโมเมนต์เป็น 0 ให้ใส่เหล็กขั้นต่ำ 2 เส้น
        if Mu_tfm <= 0.001:
            bar_counts[key] = 2
            continue

        Mu_Nmm = Mu_tfm * 9806650.0

        # คำนวณ As Required
        As_req = solve_required_as(Mu_Nmm, As_min, As_max, fc, fy, bw, d)

        # เลือกจำนวนเส้น
        bar_area = BAR_INFO[barKey]['A_cm2'] * 100
        n = math.ceil(As_req / bar_area)
        if n < 2: n = 2
        As_prov = n * bar_area

        # ตรวจสอบกำลังที่ได้จริง
        rProv = flexureSectionResponse(As_prov, fc, fy, bw, d)
        passStr = rProv['phiMn'] >= Mu_Nmm
        passMax = As_req <= As_max + 1

        # ตรวจสอบระยะเรียง (Clear Spacing)
        usable = bw - 2.0 * (cover + db_st)
        clear = (usable - n * db_main) / (n - 1) if n > 1 else usable - db_main
        req_clr = max(db_main, 25.0, 4.0 * agg_mm / 3.0)
        passClr = clear >= req_clr - 1

        overall = passStr and passMax and passClr
        if not overall: flex_ok = False
        bar_counts[key] = n

        # เพิ่มข้อมูลลงตาราง
        row(f"{title}: Mu", "-", "-", f"{fmt(Mu_tfm, 3)}", "tf-m")
        row(f"{title}: Provide", f"Use {barKey}", f"{n}-{barKey}", "OK" if overall else "NO", "-")
        row(f"{title}: Check", "φMn ≥ Mu", f"{fmt(rProv['phiMn'] / 9.8e6, 3)}", "PASS" if passStr else "FAIL", "tf-m")

    # 4. ออกแบบรับแรงเฉือน (SHEAR)
    sec("3. SHEAR DESIGN")

    Vc_N = 0.17 * math.sqrt(fc) * bw * d
    phi_v = 0.75
    phiVc_N = phi_v * Vc_N

    Av = 2.0 * BAR_INFO[stirKey]['A_cm2'] * 100
    req1 = 0.062 * math.sqrt(fc) * bw / fyt
    req2 = 0.35 * bw / fyt
    s_avmin = Av / max(req1, req2)

    VuCases = [
        {'key': "V_L", 't': "Left", 'v': inputs['vu_L']},
        {'key': "V_M", 't': "Mid", 'v': inputs['vu_M']},
        {'key': "V_R", 't': "Right", 'v': inputs['vu_R']}
    ]

    shear_res = {}
    shear_ok = True

    row("Vc", "0.17√fc' bd", "-", f"{fmt(Vc_N / 9806.65, 2)}", "tf")
    row("φVc", "0.75 * Vc", "-", f"{fmt(phiVc_N / 9806.65, 2)}", "tf")

    for case in VuCases:
        loc = case['t']
        Vu_tf = case['v']
        Vu_N = Vu_tf * 9806.65

        needVs = Vu_N > phiVc_N
        s1 = (d / 2.0) if needVs else (0.75 * d)
        smax = min(s1, 600.0)

        s_req = 9999
        if needVs:
            Vs_req_N = (Vu_N / phi_v) - Vc_N
            if Vs_req_N > 0:
                s_req = (Av * fyt * d) / Vs_req_N

        s_sel = min(s_avmin, smax)
        if needVs: s_sel = min(s_sel, s_req)

        # ปัดระยะห่างลงทีละ 25 มม. (ตามมาตรฐานก่อสร้าง)
        s_sel = math.floor(s_sel / 25.0) * 25.0
        s_sel = max(50.0, s_sel)

        Vs_prov = (Av * fyt * d) / s_sel
        phiVn = phi_v * (Vc_N + Vs_prov)
        passStr = Vu_N <= phiVn + 1

        if not passStr: shear_ok = False
        shear_res[case['key']] = s_sel

        row(f"{loc}: Vu", "-", "-", f"{fmt(Vu_tf, 3)}", "tf")
        row(f"{loc}: Provide", f"min(req, {fmt(s_avmin, 0)})", "-", f"{stirKey} @ {fmt(s_sel / 10, 0)} cm", "-",
            "OK" if passStr else "NO")

    # 5. สรุปผล
    sec("4. FINAL STATUS")
    row("Overall", "-", "-", "OK" if (flex_ok and shear_ok) else "NOT OK", "-", "")

    return calc_rows, bar_counts, shear_res


# ==========================================
# 3. PLOTTING FUNCTION
# ==========================================
def create_beam_section(b, h, cover, top_n, bot_n, stir_txt, m_db, s_db, title, bar_name):
    """ฟังก์ชันวาดรูปหน้าตัดคาน"""
    fig, ax = plt.subplots(figsize=(4, 5))

    # วาดคอนกรีต
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='#333', facecolor='#FFF')
    ax.add_patch(rect)

    # วาดเหล็กปลอก
    margin = cover + s_db / 20
    rect_s = patches.Rectangle((margin, margin), b - 2 * margin, h - 2 * margin,
                               linewidth=2, edgecolor='#2E7D32', facecolor='none', linestyle='-')
    ax.add_patch(rect_s)

    # ฟังก์ชันย่อยวาดเหล็กแกน
    def draw_row(n, y, color):
        if n < 1: return
        dia = m_db / 10
        width = b - 2 * margin - dia
        xs = [b / 2] if n == 1 else np.linspace(margin + dia / 2, b - margin - dia / 2, n)
        for x in xs:
            circle = patches.Circle((x, y), radius=dia / 2, edgecolor='black', facecolor=color)
            ax.add_patch(circle)

    # ป้องกัน Error กรณีไม่มีค่า (ให้ Default = 2 เส้น)
    top_n = top_n if top_n else 2
    bot_n = bot_n if bot_n else 2

    draw_row(top_n, h - margin - m_db / 20, '#1976D2')  # เหล็กบน (สีน้ำเงิน)
    draw_row(bot_n, margin + m_db / 20, '#D32F2F')  # เหล็กล่าง (สีแดง)

    ax.set_xlim(-5, b + 5)
    ax.set_ylim(-h * 0.2, h * 1.2)
    ax.axis('off')
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10, fontweight='bold')

    # ใส่ Text ระบุจำนวน
    ax.text(b / 2, h * 1.05, f"{top_n}-{bar_name}", ha='center', color='#1976D2', fontsize=9)
    ax.text(b / 2, -h * 0.05, f"{bot_n}-{bar_name}", ha='center', color='#D32F2F', fontsize=9)
    ax.text(b / 2, -h * 0.15, f"Stir: {stir_txt}", ha='center', color='#2E7D32', fontsize=9)

    return fig


# ==========================================
# 4. STREAMLIT USER INTERFACE
# ==========================================
# CSS สำหรับตาราง
st.markdown("""
<style>
    .report-table {width: 100%; border-collapse: collapse;}
    .report-table th, .report-table td {border: 1px solid #ddd; padding: 8px; font-size: 14px;}
    .report-table th {background-color: #f2f2f2; text-align: left;}
    .pass-ok {color: green; font-weight: bold;}
    .pass-no {color: red; font-weight: bold;}
    .sec-row {background-color: #e0e0e0; font-weight: bold;}
</style>
""", unsafe_allow_html=True)

st.title("RC Beam Designer Pro (ACI 318-19)")

# --- Input Form ---
with st.sidebar.form("inputs"):
    st.header("1. Parameters")
    c1, c2, c3 = st.columns(3)
    fc = c1.number_input("fc' (ksc)", value=240)
    fy = c2.number_input("fy (ksc)", value=4000)
    fyt = c3.number_input("fyt (ksc)", value=2400)

    c1, c2 = st.columns(2)
    mainBarKey = c1.selectbox("Main Bar", list(BAR_INFO.keys()), index=4)
    stirrupBarKey = c2.selectbox("Stirrup", list(BAR_INFO.keys()), index=0)

    st.header("2. Geometry")
    c1, c2, c3 = st.columns(3)
    b = c1.number_input("b (cm)", value=25)
    h = c2.number_input("h (cm)", value=50)
    cover = c3.number_input("cover (cm)", value=3.0)
    agg = st.number_input("Agg (mm)", value=20)

    st.header("3. Loads")
    st.markdown("**Left Support**")
    c1, c2 = st.columns(2)
    mu_L_n = c1.number_input("Mu- Top", value=8.0, key='mln')
    mu_L_p = c2.number_input("Mu+ Bot", value=4.0, key='mlp')
    vu_L = st.number_input("Vu Left", value=12.0)

    st.markdown("**Mid Span**")
    c1, c2 = st.columns(2)
    mu_M_n = c1.number_input("Mu- Top", value=0.0, key='mmn')
    mu_M_p = c2.number_input("Mu+ Bot", value=8.0, key='mmp')
    vu_M = st.number_input("Vu Mid", value=8.0)

    st.markdown("**Right Support**")
    c1, c2 = st.columns(2)
    mu_R_n = c1.number_input("Mu- Top", value=8.0, key='mrn')
    mu_R_p = c2.number_input("Mu+ Bot", value=4.0, key='mrp')
    vu_R = st.number_input("Vu Right", value=12.0)

    run_btn = st.form_submit_button("Run Calculation")

# --- Processing ---
if run_btn:
    try:
        inputs = {
            'b': b, 'h': h, 'cover': cover, 'agg': agg,
            'fc': fc, 'fy': fy, 'fyt': fyt,
            'mainBar': mainBarKey, 'stirrupBar': stirrupBarKey,
            'mu_L_n': mu_L_n, 'mu_L_p': mu_L_p,
            'mu_M_n': mu_M_n, 'mu_M_p': mu_M_p,
            'mu_R_n': mu_R_n, 'mu_R_p': mu_R_p,
            'vu_L': vu_L, 'vu_M': vu_M, 'vu_R': vu_R
        }

        rows, bars, shears = process_calculation(inputs)

        # แสดงรูปกราฟิก 3 ช่วง (ซ้าย-กลาง-ขวา)
        col1, col2, col3 = st.columns(3)
        m_db = BAR_INFO[mainBarKey]['d_mm']
        s_db = BAR_INFO[stirrupBarKey]['d_mm']

        with col1:
            st.pyplot(create_beam_section(b, h, cover, bars.get('L_TOP', 2), bars.get('L_BOT', 2),
                                          f"@{shears['V_L'] / 10:.0f}cm", m_db, s_db, "Left Support", mainBarKey))
        with col2:
            st.pyplot(create_beam_section(b, h, cover, bars.get('M_TOP', 2), bars.get('M_BOT', 2),
                                          f"@{shears['V_M'] / 10:.0f}cm", m_db, s_db, "Mid Span", mainBarKey))
        with col3:
            st.pyplot(create_beam_section(b, h, cover, bars.get('R_TOP', 2), bars.get('R_BOT', 2),
                                          f"@{shears['V_R'] / 10:.0f}cm", m_db, s_db, "Right Support", mainBarKey))

        # แสดงตารางรายการคำนวณ HTML
        st.markdown("### Calculation Report")
        html = "<table class='report-table'>"
        html += "<tr><th>Item</th><th>Formula / Ref</th><th>Substitution</th><th>Result</th><th>Unit</th><th>Status</th></tr>"
        for r in rows:
            if r[0] == "SECTION":
                html += f"<tr class='sec-row'><td colspan='6'>{r[1]}</td></tr>"
            else:
                cls = "pass-ok" if r[5] in ["OK", "PASS"] else "pass-no"
                html += f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td>{r[3]}</td><td>{r[4]}</td><td class='{cls}'>{r[5]}</td></tr>"
        html += "</table>"
        st.markdown(html, unsafe_allow_html=True)

    except Exception as e:
        st.error(f"เกิดข้อผิดพลาด: {e}")
        st.info("ลองตรวจสอบค่า Input หรือติดต่อผู้พัฒนา")
else:
    st.info("👈 กรุณากรอกข้อมูลทางด้านซ้าย แล้วกดปุ่ม 'Run Calculation'")
