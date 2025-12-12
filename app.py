import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from fpdf import FPDF
import math
import os


# -----------------------------------------------------------
# 1. Calculation Logic
# -----------------------------------------------------------
def design_beam(Mu, Vu, b, h, cover, fc, fy, fy_v):
    phi_b = 0.90
    phi_v = 0.85

    # Estimate d
    est_main_db = 2.0  # DB20
    est_stirrup_db = 0.9  # RB9
    d = h - cover - est_stirrup_db - (est_main_db / 2)

    # Flexure
    Mu_kgcm = Mu * 100
    beta1 = 0.85 if fc <= 280 else max(0.65, 0.85 - 0.05 * ((fc - 280) / 70))
    rho_b = (0.85 * beta1 * fc / fy) * (6120 / (6120 + fy))
    rho_min = 14 / fy

    Rn = Mu_kgcm / (phi_b * b * d ** 2)
    check_val = 1 - (2 * Rn / (0.85 * fc))

    if check_val < 0:
        return None  # Section too small

    rho_req = (0.85 * fc / fy) * (1 - math.sqrt(check_val))
    rho_final = max(rho_req, rho_min)
    As_req = rho_final * b * d

    # Shear
    Vc = 0.53 * math.sqrt(fc) * b * d
    phi_Vc = phi_v * Vc

    stirrup_info = "Min. Requirement"
    s_req = d / 2
    if Vu > phi_Vc:
        Vs_req = (Vu / phi_v) - Vc
        Av = 2 * (math.pi * (est_stirrup_db / 2) ** 2)
        s_calc = (Av * fy_v * d) / Vs_req
        s_req = min(s_calc, d / 2)
        stirrup_info = "Shear Reinforcement Needed"

    return {
        'd': d, 'As_req': As_req, 'rho': rho_final,
        'Vc': Vc, 'PhiVc': phi_Vc,
        'stirrup_status': stirrup_info,
        'stirrup_spacing': s_req
    }


def select_rebar_logic(As_req, b, cover):
    # Logic เลือกเหล็กแบบง่าย
    size = 20  # DB20
    area = 3.14
    num = math.ceil(As_req / area)
    num = max(num, 2)
    return num, size


# -----------------------------------------------------------
# 2. Plotting Logic
# -----------------------------------------------------------
def create_plot(b, h, cover, num_bars, main_db, stirrup_db=9):
    fig, ax = plt.subplots(figsize=(5, 5 * (h / b)))

    # Concrete
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='#333', facecolor='#F0F2F6')
    ax.add_patch(rect)

    # Stirrup
    s_dia = stirrup_db / 10
    sx, sy = cover + s_dia / 2, cover + s_dia / 2
    sw, sh = b - 2 * cover - s_dia, h - 2 * cover - s_dia
    rect_s = patches.Rectangle((sx, sy), sw, sh, linewidth=1.5, edgecolor='#0068C9', facecolor='none', linestyle='--')
    ax.add_patch(rect_s)

    # Main Bars
    m_dia = main_db / 10
    inner_w = b - 2 * cover - 2 * s_dia
    min_gap = max(2.5, m_dia)
    max_per_layer = int((inner_w + min_gap) / (m_dia + min_gap))
    if max_per_layer < 2: max_per_layer = 2

    layers = []
    rem = num_bars
    while rem > 0:
        take = min(rem, max_per_layer)
        layers.append(take)
        rem -= take

    cy = cover + s_dia + m_dia / 2
    for n in layers:
        if n == 1:
            xs = [b / 2]
        else:
            w = b - 2 * cover - 2 * s_dia - m_dia
            gap = w / (n - 1) if n > 1 else 0
            start = cover + s_dia + m_dia / 2
            xs = [start + j * gap for j in range(n)]
        for x in xs:
            c = patches.Circle((x, cy), radius=m_dia / 2, edgecolor='black', facecolor='#FF4B4B')
            ax.add_patch(c)
        cy += (m_dia + min_gap)

    plt.xlim(-5, b + 5)
    plt.ylim(-5, h + 5)
    plt.axis('off')
    ax.set_aspect('equal')
    return fig


# -----------------------------------------------------------
# 3. PDF Generator
# -----------------------------------------------------------
def create_pdf(inputs, res, num_bars, size_mm):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(0, 10, "RC Beam Design Report", 0, 1, 'C')
    pdf.line(10, 20, 200, 20)
    pdf.ln(10)

    pdf.set_font("Arial", '', 12)
    # Inputs
    pdf.set_fill_color(240, 240, 240)
    pdf.cell(0, 8, "1. Design Parameters", 0, 1, 'L', fill=True)
    pdf.cell(100, 8, f"Concrete (fc'): {inputs['fc']} ksc", 0, 0)
    pdf.cell(0, 8, f"Steel (fy): {inputs['fy']} ksc", 0, 1)
    pdf.cell(100, 8, f"Section: {inputs['b']} x {inputs['h']} cm", 0, 1)
    pdf.cell(100, 8, f"Moment (Mu): {inputs['Mu']:,} kg-m", 0, 0)
    pdf.cell(0, 8, f"Shear (Vu): {inputs['Vu']:,} kg", 0, 1)
    pdf.ln(5)

    # Results
    pdf.cell(0, 8, "2. Results", 0, 1, 'L', fill=True)
    pdf.cell(100, 8, f"Required As: {res['As_req']:.2f} cm2", 0, 0)
    pdf.cell(0, 8, f"Phi*Vc: {res['PhiVc']:.2f} kg", 0, 1)
    pdf.set_font("Arial", 'B', 12)
    pdf.set_text_color(0, 100, 0)
    pdf.cell(0, 10, f"SELECTED: {num_bars}-DB{size_mm} (Main) + RB9@{int(res['stirrup_spacing'])}cm", 0, 1)

    # Plot Image
    temp_img = "temp_plot_pdf.png"
    fig = create_plot(inputs['b'], inputs['h'], inputs['cover'], num_bars, size_mm)
    fig.savefig(temp_img, bbox_inches='tight', dpi=100)

    pdf.image(temp_img, x=60, y=None, w=90)
    if os.path.exists(temp_img): os.remove(temp_img)

    # Output logic for Streamlit
    return pdf.output(dest='S').encode('latin-1')


# -----------------------------------------------------------
# 4. Streamlit UI
# -----------------------------------------------------------
st.set_page_config(page_title="RC Beam Designer", page_icon="🏗️")

st.title("🏗️ RC Beam Design App")
st.caption("Designed according to EIT Standard (SDM)")

with st.sidebar:
    st.header("1. Material & Section")
    fc = st.number_input("Concrete fc' (ksc)", 240, 500, 240, 10)
    fy = st.number_input("Rebar fy (ksc)", 3000, 5000, 4000, 100)
    b = st.number_input("Width b (cm)", 15, 100, 25)
    h = st.number_input("Depth h (cm)", 30, 200, 50)
    cover = st.slider("Covering (cm)", 2.0, 5.0, 3.0)

    st.header("2. Loads")
    Mu = st.number_input("Moment Mu (kg-m)", 0, 100000, 12000)
    Vu = st.number_input("Shear Vu (kg)", 0, 100000, 8000)

inputs = {'fc': fc, 'fy': fy, 'b': b, 'h': h, 'cover': cover, 'Mu': Mu, 'Vu': Vu}
res = design_beam(Mu, Vu, b, h, cover, fc, fy, 2400)

if res is None:
    st.error("❌ Section too small! Please increase dimensions.")
else:
    num_bars, size_mm = select_rebar_logic(res['As_req'], b, cover)

    tab1, tab2 = st.tabs(["📊 Design Result", "📝 Full Report"])

    with tab1:
        col1, col2 = st.columns([1, 1.5])
        with col1:
            st.metric("Required As", f"{res['As_req']:.2f} cm²")
            st.metric("Main Rebar", f"{num_bars} - DB{size_mm}")
            st.metric("Stirrup", f"RB9 @ {int(res['stirrup_spacing'])} cm")
        with col2:
            fig = create_plot(b, h, cover, num_bars, size_mm)
            st.pyplot(fig)

    with tab2:
        st.write("### Ready to Export?")
        pdf_data = create_pdf(inputs, res, num_bars, size_mm)
        st.download_button(
            label="📄 Download PDF Report",
            data=pdf_data,
            file_name="Beam_Calculation.pdf",
            mime="application/pdf"
        )
