import streamlit as st
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import numpy as np
import os
import requests
from fpdf import FPDF

# ==========================================
# 1. SETUP & ROBUST FONT LOADER
# ==========================================
st.set_page_config(page_title="RC Beam Designer Pro", layout="wide")


def check_font():
    """ดาวน์โหลดฟอนต์ THSarabunNew (มี Link สำรอง)"""
    font_name = "THSarabunNew.ttf"

    # ถ้ามีไฟล์แล้ว ไม่ต้องโหลดใหม่
    if os.path.exists(font_name) and os.path.getsize(font_name) > 1000:
        return font_name

    urls = [
        "https://github.com/nutjunkie/thaifonts/raw/master/THSarabunNew.ttf",
        "https://raw.githubusercontent.com/fpdf2/fpdf2/master/test/fonts/THSarabunNew.ttf"
    ]

    for url in urls:
        try:
            r = requests.get(url, allow_redirects=True, timeout=10)
            if r.status_code == 200:
                with open(font_name, 'wb') as f:
                    f.write(r.content)
                return font_name
        except:
            continue

    return None


# ==========================================
# 2. DATABASE & HELPER
# ==========================================
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
    try:
        if n is None: return "-"
        val = float(n)
        if math.isnan(val): return "-"
        return f"{val:,.{digits}f}"
    except:
        return "-"


# ==========================================
# 3. CALCULATION LOGIC
# ==========================================
def beta1FromFc(fc_MPa):
    if fc_MPa <= 28: return 0.85
    b1 = 0.85 - 0.05 * ((fc_MPa - 28) / 7)
    return max(0.65, b1)


def phiFlexureFromStrain(eps_t):
    if eps_t <= 0.002: return 0.65
    if eps_t >= 0.005: return 0.90
    return 0.65 + (eps_t - 0.002) * (0.25 / 0.003)


def flexureSectionResponse(As_mm2, fc, fy, bw, d, Es=200000, eps_cu=0.003):
    beta1 = beta1FromFc(fc)
    fs = fy
    a = (As_mm2 * fs) / (0.85 * fc * bw) if fc > 0 else 0
    c = a / beta1 if beta1 > 0 else 0

    for i in range(50):
        if c <= 0.1: c = 0.1
        eps_t = eps_cu * (d - c) / c
        fs_new = min(fy, Es * eps_t)
        fs_new = max(fs_new, -fy)
        a_new = (As_mm2 * fs_new) / (0.85 * fc * bw)

        if abs(fs_new - fs) < 0.1 and abs(a_new - a) < 0.1:
            fs = fs_new;
            a = a_new;
            break
        fs = fs_new;
        a = a_new;
        c = a / beta1

    c = a / beta1
    eps_t = eps_cu * (d - c) / c if c > 0 else 0.005
    phi = phiFlexureFromStrain(eps_t)
    T = As_mm2 * fs
    Mn = T * (d - a / 2.0)
    phiMn = phi * Mn

    return {'phi': phi, 'phiMn': phiMn, 'eps_t': eps_t}


def solve_required_as(Mu_Nmm, As_min, As_max, fc, fy, bw, d):
    As_lo = As_min
    As_hi = As_lo
    for _ in range(30):
        r = flexureSectionResponse(As_hi, fc, fy, bw, d)
        if r['phiMn'] >= Mu_Nmm: break
        As_hi *= 1.3
        if As_hi > As_max: As_hi = As_max; break

    for _ in range(50):
        As_mid = 0.5 * (As_lo + As_hi)
        r = flexureSectionResponse(As_mid, fc, fy, bw, d)
        if r['phiMn'] >= Mu_Nmm:
            As_hi = As_mid
        else:
            As_lo = As_mid
    return As_hi


def process_calculation(inputs):
    calc_rows = []

    def sec(title):
        calc_rows.append(["SECTION", title, "", "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        calc_rows.append([item, formula, subs, result, unit, status])

    b_cm = inputs['b']
    h_cm = inputs['h']
    cover_cm = inputs['cover']
    agg_mm = inputs.get('agg', 20)
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

    sec("1. MATERIAL & SECTION PARAMETERS")
    As_min = max(0.25 * math.sqrt(fc) / fy, 1.4 / fy) * bw * d

    row("Materials", "-", f"fc'={fmt(fc, 2)} MPa", "-", "-")
    row("Section", "-", f"{fmt(bw, 0)} x {fmt(h, 0)} mm", "-", "mm")
    row("As,min", "max(0.25√fc'/fy, 1.4/fy)bd", "-", f"{fmt(As_min, 0)}", "mm²")

    sec("2. FLEXURE DESIGN")
    # ใช้ชื่อสั้นลงเพื่อป้องกันตารางล้น
    MuCases = [
        {'key': "L_TOP", 't': "Left (Top)", 'v': inputs['mu_L_n']},
        {'key': "L_BOT", 't': "Left (Bot)", 'v': inputs['mu_L_p']},
        {'key': "M_TOP", 't': "Mid (Top)", 'v': inputs['mu_M_n']},
        {'key': "M_BOT", 't': "Mid (Bot)", 'v': inputs['mu_M_p']},
        {'key': "R_TOP", 't': "Right (Top)", 'v': inputs['mu_R_n']},
        {'key': "R_BOT", 't': "Right (Bot)", 'v': inputs['mu_R_p']}
    ]

    bar_counts = {}
    flex_ok = True

    # คำนวณ As Max
    beta1 = beta1FromFc(fc)
    Es = 200000
    eps_cu = 0.003
    eps_y = fy / Es
    rho_bal = 0.85 * beta1 * (fc / fy) * (eps_cu / (eps_cu + eps_y))
    As_max = 0.75 * rho_bal * bw * d

    for case in MuCases:
        title = case['t']
        Mu_tfm = case['v']
        key = case['key']
        if Mu_tfm <= 0.001:
            bar_counts[key] = 2
            continue

        Mu_Nmm = Mu_tfm * 9806650.0
        As_req = solve_required_as(Mu_Nmm, As_min, As_max, fc, fy, bw, d)

        bar_area = BAR_INFO[barKey]['A_cm2'] * 100
        n = math.ceil(As_req / bar_area)
        if n < 2: n = 2
        As_prov = n * bar_area

        rProv = flexureSectionResponse(As_prov, fc, fy, bw, d)
        passStr = rProv['phiMn'] >= Mu_Nmm
        passMax = As_req <= As_max + 1

        usable = bw - 2.0 * (cover + db_st)
        clear = (usable - n * db_main) / (n - 1) if n > 1 else usable - db_main
        req_clr = max(db_main, 25.0, 4.0 * agg_mm / 3.0)
        passClr = clear >= req_clr - 1

        overall = passStr and passMax and passClr
        if not overall: flex_ok = False
        bar_counts[key] = n

        row(f"{title} Mu", "-", "-", f"{fmt(Mu_tfm, 3)}", "tf-m", "")
        row(f"{title} Prov", f"Use {barKey}", f"{n}-{barKey}", "OK" if overall else "NO", "-", "OK")
        row(f"{title} Check", "φMn ≥ Mu", f"{fmt(rProv['phiMn'] / 9.8e6, 3)}", "PASS" if passStr else "FAIL", "tf-m",
            "PASS")

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
    shear_ok = True
    shear_res = {}

    row("Vc", "0.17√fc' bd", "-", f"{fmt(Vc_N / 9806.65, 2)}", "tf", "")
    row("φVc", "0.75 * Vc", "-", f"{fmt(phiVc_N / 9806.65, 2)}", "tf", "")

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
        s_sel = math.floor(s_sel / 25.0) * 25.0
        s_sel = max(50.0, s_sel)
        Vs_prov = (Av * fyt * d) / s_sel
        phiVn = phi_v * (Vc_N + Vs_prov)
        passStr = Vu_N <= phiVn + 1
        if not passStr: shear_ok = False
        shear_res[case['key']] = s_sel

        row(f"{loc} Vu", "-", "-", f"{fmt(Vu_tf, 3)}", "tf", "")
        status_shear = "OK" if passStr else "NO"
        row(f"{loc} Provide", f"min(req, {fmt(s_avmin, 0)})", "-", f"{stirKey}@{fmt(s_sel / 10, 0)}cm", "-",
            status_shear)

    sec("4. FINAL STATUS")
    final_status = "OK" if (flex_ok and shear_ok) else "NOT OK"
    row("Overall", "-", "-", final_status, "-", final_status)
    return calc_rows, bar_counts, shear_res


# ==========================================
# 3. PDF GENERATION (FIXED FONT)
# ==========================================
class PDF(FPDF):
    def header(self):
        pass  # Handle manually


def create_pdf_bytes(inputs, rows, img_files):
    pdf = PDF()

    # 1. โหลดฟอนต์ (บังคับโหลดให้ได้)
    font_path = check_font()
    if not font_path:
        # ถ้าโหลดไม่ได้จริงๆ ให้ใช้ Arial (แต่จะแจ้งเตือนใน UI)
        has_thai = False
        st.error("⚠️ ไม่สามารถโหลดฟอนต์ภาษาไทยได้ รายงานจะแสดงเป็นภาษาอังกฤษ")
    else:
        has_thai = True
        pdf.add_font('THSarabunNew', '', font_path, uni=True)
        pdf.add_font('THSarabunNew', 'B', font_path, uni=True)

    def txt(s):
        s = str(s)
        if has_thai: return s
        # Safe Mode: ตัดคำไทยทิ้งถ้าไม่มีฟอนต์
        return s.encode('latin-1', 'ignore').decode('latin-1')

    pdf.add_page()

    # ตั้งค่าฟอนต์
    if has_thai:
        header_font = ('THSarabunNew', 'B', 16)
        body_font = ('THSarabunNew', '', 14)
        table_font = ('THSarabunNew', '', 12)
    else:
        header_font = ('Arial', 'B', 12)
        body_font = ('Arial', '', 10)
        table_font = ('Arial', '', 9)

    # --- Header ---
    pdf.set_font(*header_font)
    pdf.cell(0, 10, txt('ENGINEERING DESIGN REPORT'), 0, 1, 'C')
    pdf.set_font(*body_font)
    pdf.cell(0, 8, txt('Reinforced Concrete Beam Design (ACI 318-19)'), 0, 1, 'C')
    pdf.ln(5)

    # --- Project Info ---
    pdf.cell(30, 8, txt("Project:"), 0)
    pdf.cell(90, 8, txt(inputs['project']), 0)
    pdf.cell(20, 8, txt("Date:"), 0)
    pdf.cell(0, 8, "15/12/2568", 0, 1)

    pdf.cell(30, 8, txt("Engineer:"), 0)
    pdf.cell(90, 8, txt(inputs['engineer']), 0)
    pdf.cell(20, 8, txt("Code:"), 0)
    pdf.cell(0, 8, "ACI 318-19", 0, 1)
    pdf.ln(5)

    # --- Material Box ---
    pdf.set_fill_color(245, 245, 245)
    pdf.rect(10, pdf.get_y(), 90, 30, 'F')
    pdf.rect(105, pdf.get_y(), 95, 30, 'F')

    pdf.set_xy(12, pdf.get_y() + 2)
    pdf.set_font(header_font[0], 'B', 14 if has_thai else 11)
    pdf.cell(80, 8, txt("Materials"), 0, 2)
    pdf.set_font(*body_font)
    pdf.cell(80, 6, txt(f"Concrete (fc') = {inputs['fc']} ksc"), 0, 2)
    pdf.cell(80, 6, txt(f"Main Steel (fy) = {inputs['fy']} ksc"), 0, 0)

    pdf.set_xy(107, pdf.get_y() - 14)
    pdf.set_font(header_font[0], 'B', 14 if has_thai else 11)
    pdf.cell(80, 8, txt("Section"), 0, 2)
    pdf.set_font(*body_font)
    pdf.cell(80, 6, txt(f"Size = {inputs['b']} x {inputs['h']} cm"), 0, 2)
    pdf.cell(80, 6, txt(f"Cover = {inputs['cover']} cm"), 0, 0)
    pdf.ln(15)

    # --- Images ---
    pdf.ln(5)
    pdf.set_font(header_font[0], 'B', 14 if has_thai else 12)
    pdf.cell(0, 10, txt("Design Summary"), 0, 1)
    y_img = pdf.get_y()
    w_img = 60
    if len(img_files) >= 3:
        try:
            pdf.image(img_files[0], x=10, y=y_img, w=w_img)
            pdf.image(img_files[1], x=75, y=y_img, w=w_img)
            pdf.image(img_files[2], x=140, y=y_img, w=w_img)
        except:
            pass
    pdf.ln(85)

    # --- Table ---
    pdf.add_page()
    pdf.set_font(header_font[0], 'B', 14 if has_thai else 12)
    pdf.cell(0, 10, txt("Calculation Details"), 0, 1)

    # Table Header (Adjusted Widths)
    pdf.set_fill_color(220, 220, 220)
    pdf.set_font(header_font[0], 'B', 12 if has_thai else 10)
    # Col Widths: Item, Formula, Sub, Res, Unit, Status
    cols = [35, 45, 45, 20, 15, 30]
    headers = ["Item", "Formula", "Substitution", "Result", "Unit", "Status"]

    for i, h in enumerate(headers):
        pdf.cell(cols[i], 8, txt(h), 1, 0, 'C', True)
    pdf.ln()

    # Table Rows
    pdf.set_font(*table_font)
    for r in rows:
        if r[0] == "SECTION":
            pdf.set_fill_color(240, 240, 240)
            pdf.cell(sum(cols), 7, txt(r[1]), 1, 1, 'L', True)
        else:
            pdf.cell(cols[0], 7, txt(str(r[0])), 1)
            pdf.cell(cols[1], 7, txt(str(r[1])), 1)
            pdf.cell(cols[2], 7, txt(str(r[2])), 1)
            pdf.cell(cols[3], 7, txt(str(r[3])), 1)
            pdf.cell(cols[4], 7, txt(str(r[4])), 1)

            # Status
            status_txt = txt(str(r[5]))
            if "OK" in status_txt or "PASS" in status_txt:
                pdf.set_text_color(0, 100, 0)
            elif "NO" in status_txt or "FAIL" in status_txt:
                pdf.set_text_color(200, 0, 0)
            pdf.cell(cols[5], 7, status_txt, 1, 1, 'C')
            pdf.set_text_color(0, 0, 0)

    # Save
    temp = "report.pdf"
    pdf.output(temp)
    with open(temp, "rb") as f:
        data = f.read()
    try:
        os.remove(temp)
    except:
        pass
    return data


def create_beam_section(b, h, cover, top_n, bot_n, stir_txt, m_db, s_db, title, bar_name):
    fig, ax = plt.subplots(figsize=(4, 5))
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='#333', facecolor='#FFF')
    ax.add_patch(rect)
    margin = cover + s_db / 20
    rect_s = patches.Rectangle((margin, margin), b - 2 * margin, h - 2 * margin, linewidth=2, edgecolor='#2E7D32',
                               facecolor='none', linestyle='-')
    ax.add_patch(rect_s)

    def draw_row(n, y, color):
        if n < 1: return
        dia = m_db / 10
        xs = [b / 2] if n == 1 else np.linspace(margin + dia / 2, b - margin - dia / 2, n)
        for x in xs:
            circle = patches.Circle((x, y), radius=dia / 2, edgecolor='black', facecolor=color)
            ax.add_patch(circle)

    top_n = top_n if top_n else 2
    bot_n = bot_n if bot_n else 2
    draw_row(top_n, h - margin - m_db / 20, '#1976D2')
    draw_row(bot_n, margin + m_db / 20, '#D32F2F')

    ax.set_xlim(-5, b + 5)
    ax.set_ylim(-h * 0.2, h * 1.2)
    ax.axis('off')
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.text(b / 2, h * 1.05, f"{top_n}-{bar_name}", ha='center', color='#1976D2', fontsize=9)
    ax.text(b / 2, -h * 0.05, f"{bot_n}-{bar_name}", ha='center', color='#D32F2F', fontsize=9)
    ax.text(b / 2, -h * 0.15, f"Stir: {stir_txt}", ha='center', color='#2E7D32', fontsize=9)
    return fig


# ==========================================
# 4. UI MAIN
# ==========================================
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

if 'calc_done' not in st.session_state:
    st.session_state['calc_done'] = False

with st.sidebar.form("inputs"):
    st.header("Project Info")
    # --- UPDATED DEFAULT VALUES ---
    project_name = st.text_input("Project Name", value="อาคารสำนักงาน 2 ชั้น")
    engineer_name = st.text_input("Engineer Name", value="นายไกรฤทธิ์ ด่านพิทักษ์")

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
    st.markdown("**Left Support (tf-m, tf)**")
    c1, c2 = st.columns(2)
    mu_L_n = c1.number_input("Mu- Top (tf-m)", value=8.0, key='mln')
    mu_L_p = c2.number_input("Mu+ Bot (tf-m)", value=4.0, key='mlp')
    vu_L = st.number_input("Vu Left (tf)", value=12.0)

    st.markdown("**Mid Span (tf-m, tf)**")
    c1, c2 = st.columns(2)
    mu_M_n = c1.number_input("Mu- Top (tf-m)", value=0.0, key='mmn')
    mu_M_p = c2.number_input("Mu+ Bot (tf-m)", value=8.0, key='mmp')
    vu_M = st.number_input("Vu Mid (tf)", value=8.0)

    st.markdown("**Right Support (tf-m, tf)**")
    c1, c2 = st.columns(2)
    mu_R_n = c1.number_input("Mu- Top (tf-m)", value=8.0, key='mrn')
    mu_R_p = c2.number_input("Mu+ Bot (tf-m)", value=4.0, key='mrp')
    vu_R = st.number_input("Vu Right (tf)", value=12.0)

    run_btn = st.form_submit_button("Run Calculation")

if run_btn:
    inputs = {
        'project': project_name, 'engineer': engineer_name,
        'b': b, 'h': h, 'cover': cover, 'agg': agg,
        'fc': fc, 'fy': fy, 'fyt': fyt,
        'mainBar': mainBarKey, 'stirrupBar': stirrupBarKey,
        'mu_L_n': mu_L_n, 'mu_L_p': mu_L_p,
        'mu_M_n': mu_M_n, 'mu_M_p': mu_M_p,
        'mu_R_n': mu_R_n, 'mu_R_p': mu_R_p,
        'vu_L': vu_L, 'vu_M': vu_M, 'vu_R': vu_R
    }

    rows, bars, shears = process_calculation(inputs)

    st.session_state['data'] = inputs
    st.session_state['rows'] = rows
    st.session_state['bars'] = bars
    st.session_state['shears'] = shears
    st.session_state['calc_done'] = True

if st.session_state.get('calc_done'):
    data = st.session_state['data']
    rows = st.session_state['rows']
    bars = st.session_state['bars']
    shears = st.session_state['shears']

    img_files = []
    m_db = BAR_INFO[data['mainBar']]['d_mm']
    s_db = BAR_INFO[data['stirrupBar']]['d_mm']

    col1, col2, col3 = st.columns(3)
    with col1:
        fig1 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('L_TOP', 2), bars.get('L_BOT', 2),
                                   f"@{shears['V_L'] / 10:.0f}cm", m_db, s_db, "Left Support", data['mainBar'])
        st.pyplot(fig1)
        fig1.savefig("temp_L.png", dpi=100, bbox_inches='tight')
        img_files.append("temp_L.png")
        plt.close(fig1)

    with col2:
        fig2 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('M_TOP', 2), bars.get('M_BOT', 2),
                                   f"@{shears['V_M'] / 10:.0f}cm", m_db, s_db, "Mid Span", data['mainBar'])
        st.pyplot(fig2)
        fig2.savefig("temp_M.png", dpi=100, bbox_inches='tight')
        img_files.append("temp_M.png")
        plt.close(fig2)

    with col3:
        fig3 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('R_TOP', 2), bars.get('R_BOT', 2),
                                   f"@{shears['V_R'] / 10:.0f}cm", m_db, s_db, "Right Support", data['mainBar'])
        st.pyplot(fig3)
        fig3.savefig("temp_R.png", dpi=100, bbox_inches='tight')
        img_files.append("temp_R.png")
        plt.close(fig3)

    st.write("---")

    pdf_bytes = create_pdf_bytes(data, rows, img_files)

    st.download_button(
        label="🖨️ Print / Download Report (PDF)",
        data=pdf_bytes,
        file_name="Beam_Design_Report.pdf",
        mime="application/pdf"
    )

    st.markdown("### Calculation Report")
    html = "<table class='report-table'>"
    html += "<tr><th>Item</th><th>Formula</th><th>Substitution</th><th>Result</th><th>Unit</th><th>Status</th></tr>"
    for r in rows:
        if r[0] == "SECTION":
            html += f"<tr class='sec-row'><td colspan='6'>{r[1]}</td></tr>"
        else:
            cls = "pass-ok" if "OK" in r[5] or "PASS" in r[5] else "pass-no"
            html += f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td>{r[3]}</td><td>{r[4]}</td><td class='{cls}'>{r[5]}</td></tr>"
    html += "</table>"
    st.markdown(html, unsafe_allow_html=True)

else:
    st.info("👈 Please enter parameters and click 'Run Calculation'")
