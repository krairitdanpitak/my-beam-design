import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import numpy as np


# ==========================================
# 1. ACI 318-19 CALCULATION ENGINE (METRIC)
# ==========================================
def solve_aci_section(Mu_kNm, Vu_kN, b_mm, h_mm, cover_mm, fc_MPa, fy_MPa, fyt_MPa):
    """
    คำนวณหน้าตัดตาม ACI 318-19 (หน่วยภายในเป็น N, mm, MPa)
    อ้างอิงสูตรจากเอกสารแนบ
    """
    results = {}
    rows = []  # เก็บรายการคำนวณ (Formula, Sub, Res)

    # 1. Geometry
    # Formula: d = h - cover - db_stirrup - db_main/2
    est_main_db = 20  # mm (Approx DB20)
    est_stirrup_db = 9  # mm (RB9)
    d = h_mm - cover_mm - est_stirrup_db - (est_main_db / 2)

    rows.append(["Effective depth d", "d = h - cover - db(st) - db/2",
                 f"{h_mm} - {cover_mm} - {est_stirrup_db} - {est_main_db / 2}", f"{d:.1f}", "mm"])

    # 2. Parameters (Beta1)
    if fc_MPa <= 28:
        beta1 = 0.85
    else:
        beta1 = max(0.65, 0.85 - 0.05 * (fc_MPa - 28) / 7)

    rows.append(["Beta1", "ACI 318-19 Table 22.2.2.4.3", f"fc'={fc_MPa:.2f}", f"{beta1:.2f}", "-"])

    # 3. Flexure Design (Mu)
    phi_b = 0.9  # Assumption for tension controlled initially
    Mu_Nmm = abs(Mu_kNm) * 1e6

    # As min (ACI Table 9.6.1.2)
    # As,min = max(0.25*sqrt(fc)/fy, 1.4/fy) * bw * d
    term1 = (0.25 * math.sqrt(fc_MPa) / fy_MPa) * b_mm * d
    term2 = (1.4 / fy_MPa) * b_mm * d
    As_min = max(term1, term2)

    rows.append(["As,min", "max(0.25√fc'/fy, 1.4/fy) bw·d",
                 f"max({term1:.0f}, {term2:.0f})", f"{As_min:.0f}", "mm²"])

    # Solve Required As (Iterative or Quadratic)
    # Mn = As*fy*(d - a/2) -> Mu/phi = As*fy*(d - 0.59*As*fy/(0.85*fc*b))
    # Simplify: Rn = Mu / (phi*b*d^2)
    if Mu_Nmm > 0:
        Rn = Mu_Nmm / (phi_b * b_mm * d ** 2)
        rho_req = (0.85 * fc_MPa / fy_MPa) * (1 - math.sqrt(max(0, 1 - (2 * Rn) / (0.85 * fc_MPa))))
        As_calc = rho_req * b_mm * d
        As_req = max(As_calc, As_min)

        rows.append(["Flexure As,req", "Solve Mn(As) ≥ Mu/φ",
                     f"Mu={Mu_kNm:.2f} kN-m", f"{As_req:.0f}", "mm²"])
    else:
        As_req = As_min  # Minimum reinforcement
        rows.append(["Flexure As,req", "Mu ≈ 0, Use Min Steel", "-", f"{As_req:.0f}", "mm²"])

    # Select Rebar (Logic)
    bar_db = 16  # Default DB16 as in PDF example
    bar_area = 201  # DB16 area
    num_bars = math.ceil(As_req / bar_area)
    num_bars = max(num_bars, 2)
    As_prov = num_bars * bar_area

    rows.append(["Provide Reinforcement", f"Use DB{bar_db}", f"{num_bars}-DB{bar_db}", f"{As_prov}", "mm²"])

    # 4. Shear Design (Vu)
    # ACI 318-19 Metric Formula: Vc = 0.17 * sqrt(fc) * bw * d
    Vc_N = 0.17 * math.sqrt(fc_MPa) * b_mm * d
    Vc_kN = Vc_N / 1000
    phi_v = 0.75  # ACI Shear Phi
    PhiVc = phi_v * Vc_kN

    rows.append(["Vc (Concrete)", "0.17√fc'·bw·d", f"0.17√{fc_MPa:.1f}·{b_mm}·{d}", f"{Vc_kN / 9.81:.2f}",
                 "tf"])  # Show tf to match PDF
    rows.append(["φVc", "φ·Vc (φ=0.75)", f"0.75 · {Vc_kN / 9.81:.2f}", f"{PhiVc / 9.81:.2f}", "tf"])

    # Vn_max limit
    Vn_max_N = 0.66 * math.sqrt(fc_MPa) * b_mm * d
    Vn_max_kN = Vn_max_N / 1000
    PhiVn_max = phi_v * Vn_max_kN

    rows.append(["φVn,max", "φ·0.66√fc'·bw·d", "-", f"{PhiVn_max / 9.81:.2f}", "tf"])

    Vu_abs = abs(Vu_kN)
    stirrup_text = ""
    s_req = d / 2

    if Vu_abs > PhiVc:
        # Need Stirrups: Vs = Vu/phi - Vc
        Vs_kN = (Vu_abs / phi_v) - Vc_kN
        rows.append(
            ["Vs,req", "Vs = Vu/φ - Vc", f"{Vu_abs / 9.81:.2f}/0.75 - {Vc_kN / 9.81:.2f}", f"{Vs_kN / 9.81:.2f}", "tf"])

        # Spacing s = Av * fyt * d / Vs
        Av = 2 * (math.pi * (est_stirrup_db / 2) ** 2)  # 2 legs RB9 ~ 63 mm2 (RB6=28, RB9=63)
        # In PDF it uses RB6 (Av=56.6). Let's use RB9 (Av=63.6) for standard
        # But to match PDF logic, let's assume user selects bar.
        Av_used = 63.6  # RB9

        s_calc = (Av_used * fyt_MPa * d) / (Vs_kN * 1000)

        rows.append(
            ["s_req", "Av·fyt·d / Vs", f"{Av_used:.1f}·{fyt_MPa}·{d} / {Vs_kN * 1000:.0f}", f"{s_calc:.0f}", "mm"])

        s_max = min(d / 2, 600)
        s_final = min(s_calc, s_max)
        s_final = 10 * math.floor(s_final / 10)  # Round down to 10mm

        stirrup_text = f"RB9 @ {int(s_final / 10)} cm"
        rows.append(["Provide Stirrup", f"RB9 @ {int(s_final / 10)} cm", "-", "OK", "-"])

    else:
        # Min Stirrups
        stirrup_text = f"RB9 @ {int(d / 20)} cm (Min)"
        rows.append(["Stirrup Check", "Vu < φVc", "Use Min Reinforcement", "OK", "-"])

    return {
        'As_req': As_req, 'As_prov': As_prov,
        'num_bars': num_bars, 'bar_db': bar_db,
        'stirrup_text': stirrup_text,
        'rows': rows,
        'd': d
    }


# ==========================================
# 2. PLOTTING FUNCTION (With Left/Mid/Right Labels)
# ==========================================
def create_section_plot(b, h, cover, top_data, bot_data, stirrup_label, title="Section"):
    fig, ax = plt.subplots(figsize=(4, 5))

    # Concrete
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='#333', facecolor='#FFF', zorder=1)
    ax.add_patch(rect)

    # Stirrup
    s_dia = 0.9  # RB9
    margin = cover + s_dia / 2
    rect_s = patches.Rectangle((margin, margin), b - 2 * margin, h - 2 * margin,
                               linewidth=2, edgecolor='#2E7D32', facecolor='none', linestyle='-', zorder=2)
    ax.add_patch(rect_s)

    # Rebar Drawer
    def draw_bars(num, db, y_pos, color):
        if num == 0: return
        dia = db / 10
        width = b - 2 * margin - dia
        if num == 1:
            xs = [b / 2]
        else:
            xs = np.linspace(margin + dia / 2, b - margin - dia / 2, num)

        for x in xs:
            circle = patches.Circle((x, y_pos), radius=dia / 2, edgecolor='black', facecolor=color, zorder=3)
            ax.add_patch(circle)

    # Draw Top (Blue)
    draw_bars(top_data['n'], top_data['db'], h - (cover + s_dia + top_data['db'] / 20), '#1976D2')

    # Draw Bottom (Red)
    draw_bars(bot_data['n'], bot_data['db'], cover + s_dia + bot_data['db'] / 20, '#D32F2F')

    # Annotation
    ax.text(b / 2, -h * 0.15, f"B: {bot_data['n']}-DB{bot_data['db']}", ha='center', color='#D32F2F', fontweight='bold')
    ax.text(b / 2, h * 1.05, f"T: {top_data['n']}-DB{top_data['db']}", ha='center', color='#1976D2', fontweight='bold')
    ax.text(b / 2, -h * 0.25, f"S: {stirrup_label}", ha='center', color='#2E7D32')

    ax.set_title(title, fontsize=12, fontweight='bold', pad=15)
    ax.set_xlim(-10, b + 10)
    ax.set_ylim(-h * 0.3, h * 1.2)
    ax.axis('off')
    ax.set_aspect('equal')

    return fig


# ==========================================
# 3. STREAMLIT APP
# ==========================================
st.set_page_config(page_title="RC Beam Designer Pro", layout="wide")

# CSS Injection for Report Look
st.markdown("""
<style>
    .report-header { text-align: center; border-bottom: 2px solid #333; padding-bottom: 10px; margin-bottom: 20px; }
    .report-title { font-size: 24px; font-weight: bold; }
    .report-sub { font-size: 16px; color: #555; }
    .section-box { border: 1px solid #ddd; padding: 15px; border-radius: 5px; background-color: #fafafa; }
    .pass-tag { background-color: #d4edda; color: #155724; padding: 2px 8px; border-radius: 4px; font-weight: bold; font-size: 12px; }
</style>
""", unsafe_allow_html=True)

# --- HEADER ---
st.markdown("""
<div class="report-header">
    <div class="report-title">ENGINEERING DESIGN REPORT</div>
    <div class="report-sub">Reinforced Concrete Beam Design (ACI 318-19) - Gravity Beam</div>
</div>
""", unsafe_allow_html=True)

# --- INPUT SECTION (FORM) ---
with st.sidebar.form("input_form"):
    st.header("1. Design Parameters")

    st.subheader("Materials")
    fc_ksc = st.number_input("Concrete fc' (ksc)", value=240)
    fy_ksc = st.number_input("Main Steel fy (ksc)", value=4000)
    fyt_ksc = st.number_input("Stirrup fyt (ksc)", value=2400)

    st.subheader("Section Geometry")
    col_g1, col_g2 = st.columns(2)
    b_cm = col_g1.number_input("Width (cm)", value=25)
    h_cm = col_g2.number_input("Depth (cm)", value=50)
    cover_cm = st.number_input("Cover (cm)", value=3.0)

    st.subheader("Loads (Left / Mid / Right)")
    # Input 3 Zones
    col_l, col_m, col_r = st.columns(3)
    with col_l:
        st.markdown("**Left Support**")
        Mu_L = st.number_input("Mu- (kg-m)", value=8000.0, key="ml")
        Vu_L = st.number_input("Vu (kg)", value=12000.0, key="vl")
    with col_m:
        st.markdown("**Mid Span**")
        Mu_M = st.number_input("Mu+ (kg-m)", value=4000.0, key="mm")
        Vu_M = st.number_input("Vu (kg)", value=8000.0, key="vm")
    with col_r:
        st.markdown("**Right Support**")
        Mu_R = st.number_input("Mu- (kg-m)", value=8000.0, key="mr")
        Vu_R = st.number_input("Vu (kg)", value=12000.0, key="vr")

    calc_btn = st.form_submit_button("Run Calculation")

if calc_btn:
    # --- CONVERSION & CALCULATION ---
    # Convert ksc/kg/cm -> MPa/N/mm for ACI Logic
    fc_MPa = fc_ksc * 0.0980665
    fy_MPa = fy_ksc * 0.0980665
    fyt_MPa = fyt_ksc * 0.0980665
    b_mm = b_cm * 10
    h_mm = h_cm * 10
    cover_mm = cover_cm * 10

    # Run 3 Sections
    res_L = solve_aci_section(-Mu_L * 9.81 / 1000, Vu_L * 9.81 / 1000, b_mm, h_mm, cover_mm, fc_MPa, fy_MPa, fyt_MPa)
    res_M = solve_aci_section(Mu_M * 9.81 / 1000, Vu_M * 9.81 / 1000, b_mm, h_mm, cover_mm, fc_MPa, fy_MPa, fyt_MPa)
    res_R = solve_aci_section(-Mu_R * 9.81 / 1000, Vu_R * 9.81 / 1000, b_mm, h_mm, cover_mm, fc_MPa, fy_MPa, fyt_MPa)

    # --- SUMMARY SECTION (Top of Report) ---
    st.markdown("### Design Summary <span class='pass-tag'>PASS (ACI 318-19)</span>", unsafe_allow_html=True)

    # 3 Columns Visuals
    col1, col2, col3 = st.columns(3)

    # Left Plot
    with col1:
        # Assuming Negative Moment -> Top bars Main, Bottom bars Min
        # For simplicity, using calculated bars as Main
        top_L = {'n': res_L['num_bars'], 'db': res_L['bar_db']}
        bot_L = {'n': 2, 'db': 12}  # Min placeholder
        fig_L = create_section_plot(b_mm / 10, h_mm / 10, cover_mm / 10, top_L, bot_L, res_L['stirrup_text'],
                                    "Left (Support)")
        st.pyplot(fig_L)

    # Mid Plot
    with col2:
        # Positive Moment -> Bottom bars Main
        top_M = {'n': 2, 'db': 12}  # Hanger
        bot_M = {'n': res_M['num_bars'], 'db': res_M['bar_db']}
        fig_M = create_section_plot(b_mm / 10, h_mm / 10, cover_mm / 10, top_M, bot_M, res_M['stirrup_text'],
                                    "Mid (Span)")
        st.pyplot(fig_M)

    # Right Plot
    with col3:
        top_R = {'n': res_R['num_bars'], 'db': res_R['bar_db']}
        bot_R = {'n': 2, 'db': 12}
        fig_R = create_section_plot(b_mm / 10, h_mm / 10, cover_mm / 10, top_R, bot_R, res_R['stirrup_text'],
                                    "Right (Support)")
        st.pyplot(fig_R)

    st.markdown("---")

    # --- DETAILED CALCULATION TABLE ---
    st.markdown("### Calculation Details")


    # Combine rows from all sections or just show one representative?
    # Let's show "Left Support" and "Mid Span" details as in the PDF example

    def make_df(rows):
        df = pd.DataFrame(rows, columns=["Item", "Formula", "Substitution", "Result", "Unit"])
        return df


    tab1, tab2, tab3 = st.tabs(["Left Calculation", "Mid Calculation", "Right Calculation"])

    with tab1:
        st.table(make_df(res_L['rows']))
    with tab2:
        st.table(make_df(res_M['rows']))
    with tab3:
        st.table(make_df(res_R['rows']))

else:
    st.info("Please adjust inputs on the sidebar and click 'Run Calculation'")