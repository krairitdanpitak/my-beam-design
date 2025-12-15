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
# --- CHANGED: Page Title ---
st.set_page_config(page_title="RC Beam Design SDM", layout="wide")

st.markdown("""
<style>
    /* CSS สำหรับปุ่มพิมพ์ */
    .print-btn-internal {
        background-color: #008CBA;
        border: none;
        color: white !important;
        padding: 12px 28px;
        text-align: center;
        text-decoration: none;
        display: inline-block;
        font-size: 16px;
        margin: 10px 0px;
        cursor: pointer;
        border-radius: 5px;
        font-family: 'Sarabun', sans-serif;
        font-weight: bold;
        box-shadow: 0 2px 5px rgba(0,0,0,0.2);
    }
    .print-btn-internal:hover { background-color: #005f7f; }

    /* CSS สำหรับตารางในหน้าเว็บ */
    .report-table {width: 100%; border-collapse: collapse; font-family: sans-serif; font-size: 14px;}
    .report-table th, .report-table td {border: 1px solid #ddd; padding: 8px;}
    .report-table th {background-color: #f2f2f2; text-align: center; font-weight: bold;}

    .pass-ok {color: green; font-weight: bold;}
    .pass-no {color: red; font-weight: bold;}
    .sec-row {background-color: #e0e0e0; font-weight: bold; font-size: 15px;}

    /* สีแดงสำหรับ Load Input */
    .load-value {color: #D32F2F !important; font-weight: bold;}
</style>
""", unsafe_allow_html=True)

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
        val = float(n)
        if math.isnan(val): return "-"
        return f"{val:,.{digits}f}"
    except:
        return "-"


# ==========================================
# 3. CALCULATION LOGIC (DETAILED)
# ==========================================
def beta1FromFc(fc_MPa):
    if fc_MPa <= 28: return 0.85
    b1 = 0.85 - 0.05 * ((fc_MPa - 28) / 7)
    return max(0.65, b1)


def phiFlexureFromStrain(eps_t):
    if eps_t <= 0.002: return 0.65
    if eps_t >= 0.005: return 0.90
    return 0.65 + (eps_t - 0.002) * (0.25 / 0.003)


def flexureSectionResponse(As_mm2, fc_MPa, fy_MPa, bw_mm, d_mm, Es=200000, eps_cu=0.003):
    beta1 = beta1FromFc(fc_MPa)
    fs = fy_MPa
    a = (As_mm2 * fs) / (0.85 * fc_MPa * bw_mm) if fc_MPa > 0 else 0
    c = a / beta1 if beta1 > 0 else 0

    for i in range(50):
        if c <= 0.1: c = 0.1
        eps_t = eps_cu * (d_mm - c) / c
        fs_new = min(fy_MPa, Es * eps_t)
        fs_new = max(fs_new, -fy_MPa)
        a_new = (As_mm2 * fs_new) / (0.85 * fc_MPa * bw_mm)
        if abs(fs_new - fs) < 0.1 and abs(a_new - a) < 0.1:
            fs = fs_new;
            a = a_new;
            break
        fs = fs_new;
        a = a_new;
        c = a / beta1

    c = a / beta1
    eps_t = eps_cu * (d_mm - c) / c if c > 0 else 0.005
    phi = phiFlexureFromStrain(eps_t)

    T = As_mm2 * fs
    Mn = T * (d_mm - a / 2.0)
    phiMn = phi * Mn

    return {'phi': phi, 'phiMn': phiMn, 'eps_t': eps_t, 'Mn': Mn}


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
    rows = []

    def sec(title):
        rows.append(["SECTION", title, "", "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        rows.append([item, formula, subs, result, unit, status])

    b_cm = inputs['b'];
    h_cm = inputs['h'];
    cover_cm = inputs['cover']
    bw = b_cm * 10;
    h = h_cm * 10;
    cover = cover_cm * 10
    agg_mm = inputs['agg']

    fc_ksc = inputs['fc'];
    fy_ksc = inputs['fy'];
    fyt_ksc = inputs['fyt']
    ksc_to_MPa = 0.0980665
    fc_MPa = fc_ksc * ksc_to_MPa
    fy_MPa = fy_ksc * ksc_to_MPa
    fyt_MPa = fyt_ksc * ksc_to_MPa

    barKey = inputs['mainBar']
    stirKey = inputs['stirrupBar']
    db_main = BAR_INFO[barKey]['d_mm']
    db_st = BAR_INFO[stirKey]['d_mm']

    d = h - cover - db_st - (db_main / 2.0)

    # 1. MATERIALS
    sec("1. MATERIAL & SECTION PARAMETERS")
    row("Concrete & Steel", "fc', fy, fyt", f"fc'={fmt(fc_MPa, 2)} MPa, fy={fmt(fy_MPa, 0)} MPa", "-", "-")
    row("Section (b x h)", "-", f"{fmt(bw, 0)} x {fmt(h, 0)}", "-", "mm")
    row("Effective depth d", "d = h - cover - db(st) - db/2",
        f"{fmt(h, 0)} - {fmt(cover, 0)} - {fmt(db_st, 0)} - {fmt(db_main / 2, 1)}", f"{fmt(d, 1)}", "mm")

    beta1 = beta1FromFc(fc_MPa)
    row("β1", "ACI 318-19", f"fc'={fmt(fc_MPa, 2)} MPa", f"{fmt(beta1, 2)}", "-")

    term1 = 0.25 * math.sqrt(fc_MPa) / fy_MPa
    term2 = 1.4 / fy_MPa
    rho_min = max(term1, term2)
    As_min = rho_min * bw * d
    row("As,min", "max(0.25√fc'/fy, 1.4/fy) bw·d", f"max({fmt(term1, 5)}, {fmt(term2, 5)})·{fmt(bw, 0)}·{fmt(d, 0)}",
        f"{fmt(As_min, 0)}", "mm²")

    Es = 200000;
    eps_cu = 0.003;
    eps_y = fy_MPa / Es
    rho_bal = 0.85 * beta1 * (fc_MPa / fy_MPa) * (eps_cu / (eps_cu + eps_y))
    As_max = 0.75 * rho_bal * bw * d

    # 2. FLEXURE
    sec("2. FLEXURE DESIGN")
    MuCases = [
        {'key': "L_TOP", 't': "Left (Top) Mu(-)", 'v': inputs['mu_L_n']},
        {'key': "L_BOT", 't': "Left (Bot) Mu(+)", 'v': inputs['mu_L_p']},
        {'key': "M_TOP", 't': "Mid (Top) Mu(-)", 'v': inputs['mu_M_n']},
        {'key': "M_BOT", 't': "Mid (Bot) Mu(+)", 'v': inputs['mu_M_p']},
        {'key': "R_TOP", 't': "Right (Top) Mu(-)", 'v': inputs['mu_R_n']},
        {'key': "R_BOT", 't': "Right (Bot) Mu(+)", 'v': inputs['mu_R_p']}
    ]
    bar_counts = {}

    for case in MuCases:
        title = case['t']
        Mu_tfm = case['v']
        key = case['key']

        row(f"{title}: Mu", "-", "-", f"{fmt(Mu_tfm, 3)}", "tf-m", "")

        if Mu_tfm <= 0.001:
            bar_counts[key] = 2
            continue

        Mu_Nmm = Mu_tfm * 9806650.0
        As_req = solve_required_as(Mu_Nmm, As_min, As_max, fc_MPa, fy_MPa, bw, d)
        row(f"{title}: As,req", "Solve φMn(As) ≥ Mu", f"As_min={fmt(As_min, 0)}", f"{fmt(As_req, 0)}", "mm²")

        bar_area = BAR_INFO[barKey]['A_cm2'] * 100
        n = math.ceil(As_req / bar_area)
        if n < 2: n = 2
        As_prov = n * bar_area
        bar_counts[key] = n

        row(f"{title}: Provide", f"Use {barKey} (As={fmt(bar_area, 0)})", f"Req: {fmt(As_req, 0)}",
            f"{n}-{barKey} ({fmt(As_prov, 0)})", "mm²", "OK")

        res = flexureSectionResponse(As_prov, fc_MPa, fy_MPa, bw, d)
        phiMn_kNm = res['phiMn'] / 1e6
        Mu_kNm = Mu_Nmm / 1e6

        row(f"{title}: φ (Strain)", "φ = f(εt)", f"εt = {fmt(res['eps_t'], 4)}", f"{fmt(res['phi'], 2)}", "-")
        status_str = "PASS" if res['phiMn'] >= Mu_Nmm else "FAIL"
        row(f"{title}: Check Strength", "φMn ≥ Mu", f"{fmt(phiMn_kNm, 2)} ≥ {fmt(Mu_kNm, 2)} kNm", status_str, "-",
            status_str)

        usable_width = bw - 2 * (cover + db_st)
        clear_spacing = (usable_width - n * db_main) / (n - 1) if n > 1 else usable_width - db_main
        req_clr = max(db_main, 25.0, 4 / 3 * agg_mm)
        status_clr = "PASS" if clear_spacing >= req_clr - 1 else "FAIL"
        row(f"{title}: Clear Spacing", "s_clr ≥ max(db, 25, 4/3 agg)", f"{fmt(clear_spacing, 1)} ≥ {fmt(req_clr, 1)}",
            status_clr, "mm", status_clr)

    # 3. SHEAR
    sec("3. SHEAR DESIGN")
    Vc_N = 0.17 * math.sqrt(fc_MPa) * bw * d
    Vc_tf = Vc_N / 9806.65
    row("Vc", "0.17√fc' bw·d", f"0.17·√{fmt(fc_MPa, 1)}·{fmt(bw, 0)}·{fmt(d, 0)}", f"{fmt(Vc_tf, 3)}", "tf")

    phi_v = 0.75
    phiVc_tf = phi_v * Vc_tf
    row("φVc", "0.75 · Vc", f"0.75 · {fmt(Vc_tf, 3)}", f"{fmt(phiVc_tf, 3)}", "tf")

    Av = 2 * BAR_INFO[stirKey]['A_cm2'] * 100
    row("Stirrup Av", f"2-leg {stirKey}", "-", f"{fmt(Av, 1)}", "mm²")

    term1_v = 0.062 * math.sqrt(fc_MPa) * bw / fyt_MPa
    term2_v = 0.35 * bw / fyt_MPa
    Av_per_s_min = max(term1_v, term2_v)
    s_max_Av_min = Av / Av_per_s_min
    row("s(Av,min)", "Av / max(0.062√fc' bw/fyt, 0.35 bw/fyt)", f"{fmt(Av, 1)} / {fmt(Av_per_s_min, 2)}",
        f"{fmt(s_max_Av_min, 0)}", "mm")

    VuCases = [
        {'key': "V_L", 't': "Left", 'v': inputs['vu_L']},
        {'key': "V_M", 't': "Mid", 'v': inputs['vu_M']},
        {'key': "V_R", 't': "Right", 'v': inputs['vu_R']}
    ]
    shear_res = {}

    for case in VuCases:
        loc = case['t']
        Vu_tf = case['v']
        row(f"{loc}: Vu", "-", "-", f"{fmt(Vu_tf, 3)}", "tf", "")

        Vu_N = Vu_tf * 9806.65
        phiVc_N = phiVc_tf * 9806.65
        needVs = Vu_N > phiVc_N

        Vs_req_tf = 0
        if needVs:
            Vs_req_N = (Vu_N / phi_v) - Vc_N
            Vs_req_tf = Vs_req_N / 9806.65
            row(f"{loc}: Vs,req", "Vu/φ - Vc", f"{fmt(Vu_tf, 3)}/0.75 - {fmt(Vc_tf, 3)}", f"{fmt(Vs_req_tf, 3)}", "tf")
        else:
            row(f"{loc}: Vs,req", "Vu ≤ φVc", "-", "0", "tf")

        s_req = 9999
        if needVs and Vs_req_tf > 0:
            s_req = (Av * fyt_MPa * d) / (Vs_req_N)
            row(f"{loc}: s,req", "Av·fyt·d / Vs,req",
                f"{fmt(Av, 1)}·{fmt(fyt_MPa, 0)}·{fmt(d, 0)} / {fmt(Vs_req_N, 0)}", f"{fmt(s_req, 0)}", "mm")

        limit_high_shear = (Vs_req_tf * 9806.65) > (2 * Vc_N)
        s_geom_max = min(d / 4, 300) if limit_high_shear else min(d / 2, 600)
        s_limit = min(s_max_Av_min, s_geom_max)
        if needVs: s_limit = min(s_limit, s_req)

        s_prov = math.floor(s_limit / 10.0) * 10.0
        if s_prov < 50: s_prov = 50
        shear_res[case['key']] = s_prov

        row(f"{loc}: s limit", "min(s_req, s(Av,min), d/2)",
            f"min({fmt(s_req, 0)}, {fmt(s_max_Av_min, 0)}, {fmt(s_geom_max, 0)})", f"{fmt(s_limit, 0)}", "mm")
        row(f"{loc}: Provide", f"Use {stirKey}", "-", f"{stirKey} @ {fmt(s_prov / 10, 0)} cm", "-", "OK")

        Vs_prov_N = (Av * fyt_MPa * d) / s_prov
        phiVn_N = phi_v * (Vc_N + Vs_prov_N)
        phiVn_tf = phiVn_N / 9806.65
        status_shear = "PASS" if phiVn_tf >= Vu_tf - 0.05 else "FAIL"
        row(f"{loc}: Check", "φ(Vc+Vs) ≥ Vu", f"{fmt(phiVn_tf, 2)} ≥ {fmt(Vu_tf, 2)}", status_shear, "tf", status_shear)

    sec("4. FINAL STATUS")
    row("Overall", "-", "-", "DESIGN COMPLETE", "-", "OK")
    return rows, bar_counts, shear_res


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
        n = int(n);
        dia = m_db / 10
        if n == 1:
            xs = [b / 2]
        else:
            xs = np.linspace(margin + dia / 2, b - margin - dia / 2, n)
        for x in xs:
            ax.add_patch(patches.Circle((x, y), radius=dia / 2, edgecolor='black', facecolor=color))

    draw_row(int(top_n or 2), h - margin - m_db / 20, '#1976D2')
    draw_row(int(bot_n or 2), margin + m_db / 20, '#D32F2F')

    ax.set_xlim(-5, b + 5);
    ax.set_ylim(-h * 0.2, h * 1.2);
    ax.axis('off');
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.text(b / 2, h * 1.05, f"{top_n}-{bar_name}", ha='center', color='#1976D2', fontsize=9)
    ax.text(b / 2, -h * 0.05, f"{bot_n}-{bar_name}", ha='center', color='#D32F2F', fontsize=9)
    ax.text(b / 2, -h * 0.15, f"Stir: {stir_txt}", ha='center', color='#2E7D32', fontsize=9)
    return fig


def generate_full_html_report(inputs, rows, img_b64_list):
    table_rows = ""
    for r in rows:
        if r[0] == "SECTION":
            table_rows += f"<tr class='sec-row'><td colspan='6'>{r[1]}</td></tr>"
        else:
            status_cls = "pass-ok" if "OK" in r[5] or "PASS" in r[5] else "pass-no"
            val_cls = "load-value" if r[1] == "-" and (r[4] == "tf-m" or r[4] == "tf") else ""
            table_rows += f"""
            <tr>
                <td>{r[0]}</td>
                <td>{r[1]}</td>
                <td>{r[2]}</td>
                <td class='{val_cls}'>{r[3]}</td>
                <td>{r[4]}</td>
                <td class='{status_cls}'>{r[5]}</td>
            </tr>
            """

    html = f"""
    <!DOCTYPE html>
    <html lang="th">
    <head>
        <meta charset="UTF-8">
        <title>Engineering Design Report</title>
        <link href="https://fonts.googleapis.com/css2?family=Sarabun:wght@400;700&display=swap" rel="stylesheet">
        <style>
            body {{ font-family: 'Sarabun', sans-serif; padding: 20px; -webkit-print-color-adjust: exact; color: black; background: white; }}
            h1, h3 {{ text-align: center; margin: 5px; }}
            .header {{ position: relative; margin-bottom: 20px; border-bottom: 2px solid #333; padding-bottom: 10px; }}

            .beam-box {{
                position: absolute;
                top: 0;
                right: 0;
                border: 2px solid #333;
                padding: 8px 20px;
                font-size: 20px;
                font-weight: bold;
                background: #fff;
            }}

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
            .load-value {{ color: #D32F2F !important; font-weight: bold; }}

            .footer-section {{
                margin-top: 50px;
                page-break-inside: avoid;
            }}
            .signature-block {{
                width: 300px;
                text-align: center;
            }}
            .sign-line {{
                border-bottom: 1px solid #000;
                margin: 40px 0 10px 0;
                width: 100%;
            }}

            @media print {{
                .no-print {{ display: none !important; }}
                body {{ padding: 0; margin: 0; }}
                .info-box {{ border: 1px solid #333; }}
                .load-value {{ color: #D32F2F !important; -webkit-print-color-adjust: exact; }}
            }}

            .print-btn-internal {{
                background-color: #4CAF50; color: white; padding: 12px 24px;
                border: none; border-radius: 5px; cursor: pointer; font-size: 16px;
                margin-bottom: 20px; font-weight: bold; font-family: 'Sarabun', sans-serif;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            }}
            .print-btn-internal:hover {{ background-color: #45a049; }}
        </style>
    </head>
    <body>
        <div class="no-print" style="text-align: center;">
            <button onclick="window.print()" class="print-btn-internal">🖨️ Print This Page / พิมพ์หน้านี้</button>
        </div>

        <div class="header">
            <div class="beam-box">{inputs['beam_id']}</div>
            <h1>ENGINEERING DESIGN REPORT</h1>
            <h3>RC Beam Design SDM</h3>
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

        <br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

        <h3>Calculation Details</h3>
        <table>
            <thead>
                <tr>
                    <th width="20%">Item</th>
                    <th width="30%">Formula</th>
                    <th width="25%">Substitution</th>
                    <th>Result</th>
                    <th>Unit</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>

        <div class="footer-section">
            <div class="signature-block">
                <div style="text-align: left; font-weight: bold;">Designed by:</div>
                <div class="sign-line"></div>
                <div>({inputs['engineer']})</div>
                <div>วิศวกรโครงสร้าง</div>
            </div>
        </div>
    </body>
    </html>
    """
    return html


# ==========================================
# 5. UI MAIN
# ==========================================
# --- CHANGED: Main Title ---
st.title("RC Beam Design SDM")

if 'calc_done' not in st.session_state:
    st.session_state['calc_done'] = False

with st.sidebar.form("inputs"):
    st.header("Project Info")
    project_name = st.text_input("Project Name", value="อาคารสำนักงาน 2 ชั้น")
    beam_id = st.text_input("Beam Number", value="B-01")
    engineer_name = st.text_input("Engineer Name", value="นายไกรฤทธิ์ ด่านพิทักษ์")

    st.header("1. Parameters")
    c1, c2, c3 = st.columns(3)
    fc = c1.number_input("fc' (ksc)", value=240, key='fc')
    fy = c2.number_input("fy (ksc)", value=4000, key='fy')
    fyt = c3.number_input("fyt (ksc)", value=2400, key='fyt')

    c1, c2 = st.columns(2)
    mainBarKey = c1.selectbox("Main Bar", list(BAR_INFO.keys()), index=4)
    stirrupBarKey = c2.selectbox("Stirrup", list(BAR_INFO.keys()), index=0)

    st.header("2. Geometry")
    c1, c2, c3 = st.columns(3)
    b = c1.number_input("b (cm)", value=25, key='b')
    h = c2.number_input("h (cm)", value=50, key='h')
    cover = c3.number_input("cover (cm)", value=3.0, key='cover')
    agg = st.number_input("Agg (mm)", value=20, key='agg')

    st.header("3. Loads")
    st.markdown("**Left Support (Mu: tf-m, Vu: tf)**")
    c1, c2 = st.columns(2)
    mu_L_n = c1.number_input("Mu- Top (tf-m)", value=8.0, key='L_mn')
    mu_L_p = c2.number_input("Mu+ Bot (tf-m)", value=4.0, key='L_mp')
    vu_L = st.number_input("Vu Left (tf)", value=12.0, key='L_v')

    st.markdown("**Mid Span (Mu: tf-m, Vu: tf)**")
    c1, c2 = st.columns(2)
    mu_M_n = c1.number_input("Mu- Top (tf-m)", value=0.0, key='M_mn')
    mu_M_p = c2.number_input("Mu+ Bot (tf-m)", value=8.0, key='M_mp')
    vu_M = st.number_input("Vu Mid (tf)", value=8.0, key='M_v')

    st.markdown("**Right Support (Mu: tf-m, Vu: tf)**")
    c1, c2 = st.columns(2)
    mu_R_n = c1.number_input("Mu- Top (tf-m)", value=8.0, key='R_mn')
    mu_R_p = c2.number_input("Mu+ Bot (tf-m)", value=4.0, key='R_mp')
    vu_R = st.number_input("Vu Right (tf)", value=12.0, key='R_v')

    run_btn = st.form_submit_button("Run Calculation")

if run_btn:
    inputs = {
        'project': project_name, 'beam_id': beam_id,
        'engineer': engineer_name,
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

    # 1. Generate Images
    fig1 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('L_TOP', 2), bars.get('L_BOT', 2),
                               f"@{shears['V_L'] / 10:.0f}cm", m_db, s_db, "Left Support", data['mainBar'])
    img1_b64 = fig_to_base64(fig1)
    fig2 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('M_TOP', 2), bars.get('M_BOT', 2),
                               f"@{shears['V_M'] / 10:.0f}cm", m_db, s_db, "Mid Span", data['mainBar'])
    img2_b64 = fig_to_base64(fig2)
    fig3 = create_beam_section(data['b'], data['h'], data['cover'], bars.get('R_TOP', 2), bars.get('R_BOT', 2),
                               f"@{shears['V_R'] / 10:.0f}cm", m_db, s_db, "Right Support", data['mainBar'])
    img3_b64 = fig_to_base64(fig3)

    # 2. Generate HTML Report
    html_report = generate_full_html_report(data, rows, [img1_b64, img2_b64, img3_b64])

    st.success("✅ คำนวณเสร็จสิ้น (Calculation Finished)")
    st.components.v1.html(html_report, height=800, scrolling=True)

else:
    st.info("👈 กรุณากรอกข้อมูลทางด้านซ้าย แล้วกดปุ่ม 'Run Calculation'")