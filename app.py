import streamlit as st
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import numpy as np
import io
import base64
import streamlit.components.v1 as components

# ==========================================
# 1. SETUP & CSS
# ==========================================
st.set_page_config(page_title="RC Beam Designer Pro", layout="wide")

# CSS สำหรับปุ่ม Link และตาราง
st.markdown("""
<style>
    .print-btn {
        background-color: #008CBA;
        border: none;
        color: white !important;
        padding: 12px 28px;
        text-align: center;
        text-decoration: none;
        display: inline-block;
        font-size: 18px;
        margin: 10px 0px;
        cursor: pointer;
        border-radius: 5px;
        font-weight: bold;
        box-shadow: 0 2px 5px rgba(0,0,0,0.2);
        transition: background-color 0.3s;
    }
    .print-btn:hover {
        background-color: #005f7f;
        color: white !important;
    }
    .report-table {width: 100%; border-collapse: collapse; font-family: sans-serif;}
    .report-table th, .report-table td {border: 1px solid #ddd; padding: 8px; font-size: 14px;}
    .report-table th {background-color: #f2f2f2; text-align: left; font-weight: bold;}
    .pass-ok {color: green; font-weight: bold;}
    .pass-no {color: red; font-weight: bold;}
    .sec-row {background-color: #e0e0e0; font-weight: bold; font-size: 15px;}
</style>
""", unsafe_allow_html=True)

# ==========================================
# 2. DATABASE
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
        overall = passStr and passMax
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
# 4. PLOTTING & REPORT GEN
# ==========================================
def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}"


def create_beam_section(b, h, cover, top_n, bot_n, stir_txt, m_db, s_db, title, bar_name):
    fig, ax = plt.subplots(figsize=(3, 4))
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='#333', facecolor='#FFF')
    ax.add_patch(rect)
    margin = cover + s_db / 20
    rect_s = patches.Rectangle((margin, margin), b - 2 * margin, h - 2 * margin, linewidth=2, edgecolor='#2E7D32',
                               facecolor='none', linestyle='-')
    ax.add_patch(rect_s)

    def draw_row(n, y, color):
        n = int(n)
        if n < 1: return
        dia = m_db / 10
        if n == 1:
            xs = [b / 2]
        else:
            xs = np.linspace(margin + dia / 2, b - margin - dia / 2, n)
        for x in xs:
            circle = patches.Circle((x, y), radius=dia / 2, edgecolor='black', facecolor=color)
            ax.add_patch(circle)

    top_n = int(top_n) if top_n else 2
    bot_n = int(bot_n) if bot_n else 2
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


def generate_full_html_report(inputs, rows, img_b64_list, auto_print=False):
    """สร้าง HTML Template พร้อมปุ่มกดในตัว"""

    # 1. Table
    table_rows = ""
    for r in rows:
        if r[0] == "SECTION":
            table_rows += f"<tr class='sec-row'><td colspan='6'>{r[1]}</td></tr>"
        else:
            status_cls = "pass-ok" if "OK" in r[5] or "PASS" in r[5] else "pass-no"
            table_rows += f"""
            <tr>
                <td>{r[0]}</td>
                <td>{r[1]}</td>
                <td>{r[2]}</td>
                <td>{r[3]}</td>
                <td>{r[4]}</td>
                <td class='{status_cls}'>{r[5]}</td>
            </tr>
            """

    # 2. Script: Auto print only if requested
    print_script = ""
    if auto_print:
        print_script = """
        <script>
            window.onload = function() {
                setTimeout(function(){ window.print(); }, 800);
            }
        </script>
        """

    # 3. HTML Content
    html_content = f"""
    <!DOCTYPE html>
    <html lang="th">
    <head>
        <meta charset="UTF-8">
        <title>Engineering Design Report</title>
        <link href="https://fonts.googleapis.com/css2?family=Sarabun:wght@400;700&display=swap" rel="stylesheet">
        <style>
            body {{ font-family: 'Sarabun', sans-serif; padding: 20px; -webkit-print-color-adjust: exact; }}
            h1, h3 {{ text-align: center; margin: 5px; }}
            .header {{ margin-bottom: 20px; border-bottom: 2px solid #333; padding-bottom: 10px; }}
            .info-container {{ display: flex; justify-content: space-between; margin-bottom: 20px; }}
            .info-box {{ width: 48%; border: 1px solid #ddd; padding: 10px; border-radius: 5px; background: #f9f9f9; }}
            .images {{ display: flex; justify-content: space-around; margin: 20px 0; }}
            .images img {{ width: 30%; border: 1px solid #ddd; padding: 5px; }}
            table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 12px; }}
            th, td {{ border: 1px solid #444; padding: 6px; }}
            th {{ background-color: #eee; }}
            .sec-row {{ background-color: #ddd; font-weight: bold; }}
            .pass-ok {{ color: green; font-weight: bold; text-align: center; }}
            .pass-no {{ color: red; font-weight: bold; text-align: center; }}

            @media print {{
                .no-print {{ display: none !important; }}
                body {{ padding: 0; }}
            }}
            .print-btn-internal {{
                background-color: #4CAF50; color: white; padding: 10px 20px;
                border: none; border-radius: 5px; cursor: pointer; font-size: 16px;
                margin-bottom: 20px;
            }}
        </style>
    </head>
    <body>
        <div class="no-print" style="text-align: center;">
            <button onclick="window.print()" class="print-btn-internal">🖨️ Print This Page / พิมพ์หน้านี้</button>
        </div>

        <div class="header">
            <h1>ENGINEERING DESIGN REPORT</h1>
            <h3>Reinforced Concrete Beam Design (ACI 318-19)</h3>
        </div>

        <div class="info-container">
            <div class="info-box">
                <strong>Project:</strong> {inputs['project']}<br>
                <strong>Engineer:</strong> {inputs['engineer']}<br>
                <strong>Date:</strong> 15/12/2568
            </div>
            <div class="info-box">
                <strong>Materials:</strong> fc'={inputs['fc']} ksc, fy={inputs['fy']} ksc<br>
                <strong>Section:</strong> {inputs['b']} x {inputs['h']} cm, Cover={inputs['cover']} cm
            </div>
        </div>

        <h3>Design Summary</h3>
        <div class="images">
            <img src="{img_b64_list[0]}" />
            <img src="{img_b64_list[1]}" />
            <img src="{img_b64_list[2]}" />
        </div>

        <h3>Calculation Details</h3>
        <table>
            <thead>
                <tr>
                    <th width="25%">Item</th>
                    <th width="30%">Formula</th>
                    <th width="20%">Substitution</th>
                    <th>Result</th>
                    <th>Unit</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>

        {print_script}
    </body>
    </html>
    """
    return html_content


# ==========================================
# 5. UI MAIN
# ==========================================
st.title("RC Beam Designer Pro (ACI 318-19)")

if 'calc_done' not in st.session_state:
    st.session_state['calc_done'] = False

with st.sidebar.form("inputs"):
    st.header("Project Info")
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
    c1, c2 = st.columns(2)
    mu_L_n = c1.number_input("Mu- Top (Left)", value=8.0)
    mu_L_p = c2.number_input("Mu+ Bot (Left)", value=4.0)
    vu_L = st.number_input("Vu Left", value=12.0)

    c1, c2 = st.columns(2)
    mu_M_n = c1.number_input("Mu- Top (Mid)", value=0.0)
    mu_M_p = c2.number_input("Mu+ Bot (Mid)", value=8.0)
    vu_M = st.number_input("Vu Mid", value=8.0)

    c1, c2 = st.columns(2)
    mu_R_n = c1.number_input("Mu- Top (Right)", value=8.0)
    mu_R_p = c2.number_input("Mu+ Bot (Right)", value=4.0)
    vu_R = st.number_input("Vu Right", value=12.0)

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

if st.session_state['calc_done']:
    data = st.session_state['data']
    rows = st.session_state['rows']
    bars = st.session_state['bars']
    shears = st.session_state['shears']
    m_db = BAR_INFO[data['mainBar']]['d_mm']
    s_db = BAR_INFO[data['stirrupBar']]['d_mm']

    # 1. GRAPHICS
    fig1 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('L_TOP', 2), bars.get('L_BOT', 2),
                               f"@{shears['V_L'] / 10:.0f}cm", m_db, s_db, "Left Support", data['mainBar'])
    img1_b64 = fig_to_base64(fig1)

    fig2 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('M_TOP', 2), bars.get('M_BOT', 2),
                               f"@{shears['V_M'] / 10:.0f}cm", m_db, s_db, "Mid Span", data['mainBar'])
    img2_b64 = fig_to_base64(fig2)

    fig3 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('R_TOP', 2), bars.get('R_BOT', 2),
                               f"@{shears['V_R'] / 10:.0f}cm", m_db, s_db, "Right Support", data['mainBar'])
    img3_b64 = fig_to_base64(fig3)

    # 2. GENERATE TWO VERSIONS OF HTML
    # Version A: Preview (No auto print)
    html_preview = generate_full_html_report(data, rows, [img1_b64, img2_b64, img3_b64], auto_print=False)

    # Version B: Link (Auto print)
    html_link = generate_full_html_report(data, rows, [img1_b64, img2_b64, img3_b64], auto_print=True)

    # Encode Version B for Link
    b64_html = base64.b64encode(html_link.encode('utf-8')).decode('utf-8')
    href = f'data:text/html;charset=utf-8;base64,{b64_html}'

    # 3. DISPLAY
    st.markdown("### ✅ คำนวณเสร็จสิ้น (Calculation Finished)")

    # PRINT BUTTON (เปิดแท็บใหม่ + สั่งพิมพ์อัตโนมัติ)
    st.markdown(f"""
        <div style="text-align: center; margin: 20px;">
            <a href="{href}" target="_blank" class="print-btn">
                🖨️ เปิดหน้ารายงานเพื่อพิมพ์ (Open Print Report)
            </a>
        </div>
    """, unsafe_allow_html=True)

    st.write("---")

    # SHOW PREVIEW (แบบไม่เด้ง Print)
    st.subheader("Preview (ตัวอย่างรายงาน)")
    st.components.v1.html(html_preview, height=600, scrolling=True)

else:
    st.info("👈 กรุณากรอกข้อมูลทางด้านซ้าย แล้วกดปุ่ม 'Run Calculation'")
