import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from fpdf import FPDF
import math
import os


# -----------------------------------------------------------
# 1. Calculation Logic (เพิ่มการเก็บ Step และเหล็กบน)
# -----------------------------------------------------------
def design_beam(Mu, Vu, b, h, cover, fc, fy, fy_v):
    steps = []  # เก็บข้อความรายการคำนวณ

    # 1. Constants & Properties
    phi_b = 0.90
    phi_v = 0.85
    steps.append(f"**1. Properties:** fc'={fc} ksc, fy={fy} ksc, Size={b}x{h} cm")

    # 2. Estimate d
    est_main_db = 2.0  # DB20
    est_stirrup_db = 0.9  # RB9
    d = h - cover - est_stirrup_db - (est_main_db / 2)
    steps.append(f"**2. Effective Depth (d):**")
    steps.append(f"   d = h - cover - stirrup - (db/2)")
    steps.append(f"   d = {h} - {cover} - {est_stirrup_db} - ({est_main_db}/2) = {d:.2f} cm")

    # 3. Flexure Design (Singly Reinforced)
    Mu_kgcm = Mu * 100
    steps.append(f"**3. Flexural Design:**")

    # Beta1
    beta1 = 0.85 if fc <= 280 else max(0.65, 0.85 - 0.05 * ((fc - 280) / 70))

    # Rho Balanced
    rho_b = (0.85 * beta1 * fc / fy) * (6120 / (6120 + fy))
    rho_max = 0.75 * rho_b
    rho_min = 14 / fy
    steps.append(f"   rho_min = {rho_min:.5f}, rho_max = {rho_max:.5f}")

    # Rn
    Rn = Mu_kgcm / (phi_b * b * d ** 2)
    steps.append(f"   Rn = Mu / (phi * b * d^2) = {Mu_kgcm} / ({phi_b}*{b}*{d:.2f}^2) = {Rn:.2f} ksc")

    # Check Section
    check_val = 1 - (2 * Rn / (0.85 * fc))
    if check_val < 0:
        return None, "Error: Section too small (Rn too high)"

    # Rho Required
    rho_req = (0.85 * fc / fy) * (1 - math.sqrt(check_val))
    rho_final = max(rho_req, rho_min)
    As_req = rho_final * b * d

    steps.append(f"   rho_req = {rho_req:.5f} -> Use rho = {rho_final:.5f}")
    steps.append(f"   As_req = rho * b * d = {rho_final:.5f} * {b} * {d:.2f} = **{As_req:.2f} cm²**")

    # 4. Shear Design
    steps.append(f"**4. Shear Design:**")
    Vc = 0.53 * math.sqrt(fc) * b * d
    phi_Vc = phi_v * Vc
    steps.append(f"   Vc = 0.53 * sqrt(fc') * b * d = {Vc:.2f} kg")
    steps.append(f"   Phi*Vc = {phi_Vc:.2f} kg (Vu = {Vu} kg)")

    stirrup_info = ""
    s_req = d / 2

    if Vu > phi_Vc:
        Vs_req = (Vu / phi_v) - Vc
        Av = 2 * (math.pi * (est_stirrup_db / 2) ** 2)  # 2 legs RB9
        s_calc = (Av * fy_v * d) / Vs_req
        s_req = min(s_calc, d / 2)
        stirrup_info = f"RB9 @ {int(s_req)} cm (Vs Required)"
        steps.append(f"   result: Vu > Phi*Vc -> Need Shear Rebar")
        steps.append(f"   Vs = (Vu/phi) - Vc = {Vs_req:.2f} kg")
        steps.append(f"   Spacing (s) = (Av*fy*d)/Vs = {s_calc:.2f} cm -> Use {int(s_req)} cm")
    else:
        s_req = d / 2
        stirrup_info = f"RB9 @ {int(s_req)} cm (Min Control)"
        steps.append(f"   result: Vu <= Phi*Vc -> Min Stirrup")
        steps.append(f"   Max Spacing = d/2 = {d / 2:.2f} cm -> Use {int(s_req)} cm")

    return {
        'd': d, 'As_req': As_req, 'rho': rho_final,
        'Vc': Vc, 'PhiVc': phi_Vc,
        'stirrup_info': stirrup_info,
        'stirrup_spacing': s_req
    }, steps


def select_rebar_logic(As_req, b):
    # 1. Bottom Bars (Main)
    main_db = 20  # สมมติเลือก DB20 เป็นหลัก
    area_main = 3.14
    num_bottom = math.ceil(As_req / area_main)
    num_bottom = max(num_bottom, 2)  # ขั้นต่ำ 2 เส้น

    # 2. Top Bars (Hanger) - ใส่เพื่อยึดเหล็กปลอก
    # กฎง่ายๆ: ถ้าคานเล็กใช้ 2-DB12, ถ้าคานใหญ่ใช้ 2-DB16
    top_db = 12
    if b > 30: top_db = 16
    num_top = 2

    return {
        'bottom_num': num_bottom, 'bottom_db': main_db,
        'top_num': num_top, 'top_db': top_db
    }


# -----------------------------------------------------------
# 2. Plotting Logic (เพิ่มเหล็กบน)
# -----------------------------------------------------------
def create_plot(b, h, cover, rebar_data, stirrup_db=9):
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

    # Helper to draw bars
    def draw_layer(num, db, y_center, color):
        dia = db / 10
        inner_w = b - 2 * cover - 2 * s_dia
        min_gap = max(2.5, dia)
        max_per_layer = int((inner_w + min_gap) / (dia + min_gap))
        if max_per_layer < 2: max_per_layer = 2

        layers = []
        rem = num
        while rem > 0:
            take = min(rem, max_per_layer)
            layers.append(take)
            rem -= take

        curr_y = y_center
        # Note: This simple logic draws bottom-up for bottom bars.
        # For top bars, we might need simple single layer for now.

        for n in layers:
            if n == 1:
                xs = [b / 2]
            else:
                w = b - 2 * cover - 2 * s_dia - dia
                gap = w / (n - 1) if n > 1 else 0
                start = cover + s_dia + dia / 2
                xs = [start + j * gap for j in range(n)]

            for x in xs:
                c = patches.Circle((x, curr_y), radius=dia / 2, edgecolor='black', facecolor=color)
                ax.add_patch(c)
            # If multiple layers, move 'curr_y' (Logic applies mostly to bottom)
            if y_center < h / 2:  # Bottom bars
                curr_y += (dia + min_gap)
            else:  # Top bars (Move down if multiple layers - simplified here)
                curr_y -= (dia + min_gap)

    # Draw Bottom Bars (Red)
    bot_dia = rebar_data['bottom_db'] / 10
    y_bot = cover + s_dia + bot_dia / 2
    draw_layer(rebar_data['bottom_num'], rebar_data['bottom_db'], y_bot, '#FF4B4B')

    # Draw Top Bars (Blue/Green)
    top_dia = rebar_data['top_db'] / 10
    y_top = h - (cover + s_dia + top_dia / 2)
    draw_layer(rebar_data['top_num'], rebar_data['top_db'], y_top, '#2E7D32')

    plt.xlim(-5, b + 5)
    plt.ylim(-5, h + 5)
    plt.axis('off')
    ax.set_aspect('equal')
    return fig


# -----------------------------------------------------------
# 3. PDF Generator (เพิ่ม Calculation Steps)
# -----------------------------------------------------------
def create_pdf(inputs, res, rebar_data, calc_steps):
    pdf = FPDF()
    pdf.add_page()

    # Header
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(0, 10, "RC Beam Design Calculation", 0, 1, 'C')
    pdf.set_font("Arial", '', 10)
    pdf.cell(0, 5, "Based on EIT Standard (SDM)", 0, 1, 'C')
    pdf.line(10, 25, 200, 25)
    pdf.ln(10)

    # 1. Summary Box
    pdf.set_fill_color(240, 240, 240)
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 8, "  Design Summary", 0, 1, 'L', fill=True)
    pdf.set_font("Arial", '', 11)

    summary_text = (
        f"Beam Size: {inputs['b']} x {inputs['h']} cm\n"
        f"Bottom Rebar: {rebar_data['bottom_num']}-DB{rebar_data['bottom_db']}\n"
        f"Top Rebar: {rebar_data['top_num']}-DB{rebar_data['top_db']}\n"
        f"Stirrup: {res['stirrup_info']}"
    )
    pdf.multi_cell(0, 6, summary_text)
    pdf.ln(5)

    # 2. Detailed Calculation
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 8, "  Detailed Calculation Steps", 0, 1, 'L', fill=True)
    pdf.ln(2)
    pdf.set_font("Arial", '', 10)

    for line in calc_steps:
        # Clean markdown bold symbols for PDF
        clean_line = line.replace("**", "")
        pdf.cell(0, 6, clean_line, 0, 1)

    # 3. Section Image
    pdf.ln(5)
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 8, "  Section Detail", 0, 1, 'L', fill=True)

    temp_img = "temp_plot.png"
    fig = create_plot(inputs['b'], inputs['h'], inputs['cover'], rebar_data)
    fig.savefig(temp_img, bbox_inches='tight', dpi=100)

    # Center image
    pdf.image(temp_img, x=70, y=None, w=70)
    if os.path.exists(temp_img): os.remove(temp_img)

    return pdf.output(dest='S').encode('latin-1')


# -----------------------------------------------------------
# 4. Streamlit UI
# -----------------------------------------------------------
st.set_page_config(page_title="RC Beam Designer Pro", page_icon="🏗️", layout="wide")

st.title("🏗️ RC Beam Design Pro")
st.caption("Complete Calculation with Steps & Detailing")

# Sidebar
with st.sidebar:
    st.header("Materials & Geometry")
    fc = st.number_input("Concrete fc' (ksc)", 240, 500, 240, 10)
    fy = st.number_input("Main Steel fy (ksc)", 3000, 5000, 4000, 100)
    b = st.number_input("Width b (cm)", 15, 100, 25)
    h = st.number_input("Depth h (cm)", 30, 200, 50)
    cover = st.slider("Covering (cm)", 2.0, 5.0, 3.0)

    st.header("Design Loads")
    Mu = st.number_input("Factored Moment Mu (kg-m)", 0, 100000, 12000)
    Vu = st.number_input("Factored Shear Vu (kg)", 0, 100000, 8000)

# Process
inputs = {'fc': fc, 'fy': fy, 'b': b, 'h': h, 'cover': cover, 'Mu': Mu, 'Vu': Vu}
res, steps = design_beam(Mu, Vu, b, h, cover, fc, fy, 2400)

if res is None:
    st.error(steps)  # Show error message
else:
    # Select Rebar
    rebar_data = select_rebar_logic(res['As_req'], b)

    col_main, col_detail = st.columns([1.2, 1])

    with col_main:
        st.subheader("✅ Design Results")
        c1, c2 = st.columns(2)
        c1.info(f"**Bottom (Main):** {rebar_data['bottom_num']} - DB{rebar_data['bottom_db']}")
        c2.success(f"**Top (Hanger):** {rebar_data['top_num']} - DB{rebar_data['top_db']}")
        st.warning(f"**Stirrup:** {res['stirrup_info']}")

        # Calculation Steps Expander
        with st.expander("📝 ดูรายการคำนวณละเอียด (Show Calculation Steps)", expanded=True):
            for line in steps:
                st.markdown(line)

    with col_detail:
        st.subheader("📐 Section Detailing")
        fig = create_plot(b, h, cover, rebar_data)
        st.pyplot(fig)

        # Download
        pdf_bytes = create_pdf(inputs, res, rebar_data, steps)
        st.download_button(
            label="📄 Download Calculation Report (PDF)",
            data=pdf_bytes,
            file_name="Beam_Design_Full.pdf",
            mime="application/pdf",
            use_container_width=True
        )
