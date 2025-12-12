import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from fpdf import FPDF
import math
import os


# -----------------------------------------------------------
# 1. Calculation Logic (ส่วนคำนวณ)
# -----------------------------------------------------------
def design_beam(Mu, Vu, b, h, cover, fc, fy, fy_v):
    steps = []

    # 1. Properties
    phi_b = 0.90
    phi_v = 0.85
    steps.append(f"**1. Properties:** fc'={fc} ksc, fy={fy} ksc, Size={b}x{h} cm")

    # 2. Estimate d
    est_main_db = 2.0  # DB20
    est_stirrup_db = 0.9  # RB9
    d = h - cover - est_stirrup_db - (est_main_db / 2)
    steps.append(f"**2. Effective Depth (d):** {d:.2f} cm")

    # 3. Flexure
    Mu_kgcm = Mu * 100
    beta1 = 0.85 if fc <= 280 else max(0.65, 0.85 - 0.05 * ((fc - 280) / 70))
    rho_b = (0.85 * beta1 * fc / fy) * (6120 / (6120 + fy))
    rho_min = 14 / fy

    Rn = Mu_kgcm / (phi_b * b * d ** 2)
    check_val = 1 - (2 * Rn / (0.85 * fc))

    if check_val < 0:
        return None, "Error: หน้าตัดเล็กเกินไป (Section too small)"

    rho_req = (0.85 * fc / fy) * (1 - math.sqrt(check_val))
    rho_final = max(rho_req, rho_min)
    As_req = rho_final * b * d

    steps.append(f"**3. Flexure:** As required = **{As_req:.2f} cm²** (Rho={rho_final:.4f})")

    # 4. Shear
    Vc = 0.53 * math.sqrt(fc) * b * d
    phi_Vc = phi_v * Vc

    stirrup_info = ""
    s_req = d / 2
    if Vu > phi_Vc:
        Vs_req = (Vu / phi_v) - Vc
        Av = 2 * (math.pi * (est_stirrup_db / 2) ** 2)
        s_calc = (Av * fy_v * d) / Vs_req
        s_req = min(s_calc, d / 2)
        stirrup_info = f"RB9 @ {int(s_req)} cm (Shear Reinf.)"
    else:
        s_req = d / 2
        stirrup_info = f"RB9 @ {int(s_req)} cm (Min Control)"

    steps.append(f"**4. Shear:** {stirrup_info}")

    return {
        'd': d, 'As_req': As_req, 'rho': rho_final,
        'Vc': Vc, 'PhiVc': phi_Vc,
        'stirrup_info': stirrup_info,
        'stirrup_spacing': s_req
    }, steps


def select_rebar_logic(As_req, b):
    # Bottom Bars
    main_db = 20  # DB20
    area_main = 3.14
    num_bottom = math.ceil(As_req / area_main)
    num_bottom = max(num_bottom, 2)

    # Top Bars (Hanger)
    top_db = 12 if b <= 30 else 16
    num_top = 2

    return {
        'bottom_num': num_bottom, 'bottom_db': main_db,
        'top_num': num_top, 'top_db': top_db
    }


# -----------------------------------------------------------
# 2. Plotting Logic (เพิ่ม Text บอกขนาดเหล็กในรูป)
# -----------------------------------------------------------
def create_plot(b, h, cover, rebar_data, stirrup_db=9):
    # ปรับขนาดรูปให้มีที่ว่างขอบๆ สำหรับใส่ตัวหนังสือ
    fig, ax = plt.subplots(figsize=(6, 6 * (h / b)))

    # Concrete
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='#333', facecolor='#F5F5F5')
    ax.add_patch(rect)

    # Stirrup
    s_dia = stirrup_db / 10
    sx, sy = cover + s_dia / 2, cover + s_dia / 2
    sw, sh = b - 2 * cover - s_dia, h - 2 * cover - s_dia
    rect_s = patches.Rectangle((sx, sy), sw, sh, linewidth=1.5, edgecolor='#0068C9', facecolor='none', linestyle='--')
    ax.add_patch(rect_s)

    # Helper draw function
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

            # Move Y for next layer
            if y_center < h / 2:
                curr_y += (dia + min_gap)  # Bottom move up
            else:
                curr_y -= (dia + min_gap)  # Top move down

    # Draw Bottom
    bot_db = rebar_data['bottom_db']
    y_bot = cover + s_dia + (bot_db / 10) / 2
    draw_layer(rebar_data['bottom_num'], bot_db, y_bot, '#D32F2F')  # Red

    # Draw Top
    top_db = rebar_data['top_db']
    y_top = h - (cover + s_dia + (top_db / 10) / 2)
    draw_layer(rebar_data['top_num'], top_db, y_top, '#2E7D32')  # Green

    # --- เพิ่ม Text บอกขนาดเหล็ก (Labels) ---
    # Top Label
    label_top = f"{rebar_data['top_num']}-DB{top_db}"
    ax.text(b / 2, h + 3, label_top, ha='center', va='bottom', fontsize=14, color='#2E7D32', fontweight='bold')

    # Bottom Label
    label_bot = f"{rebar_data['bottom_num']}-DB{bot_db}"
    ax.text(b / 2, -3, label_bot, ha='center', va='top', fontsize=14, color='#D32F2F', fontweight='bold')

    # Stirrup Label (Side)
    label_stirrup = f"RB{stirrup_db}"
    ax.text(-2, h / 2, label_stirrup, ha='right', va='center', fontsize=10, color='#0068C9', rotation=90)

    # Set Plot Limits (เผื่อที่ให้ Text)
    plt.xlim(-8, b + 8)
    plt.ylim(-8, h + 8)
    plt.axis('off')
    ax.set_aspect('equal')
    return fig


# -----------------------------------------------------------
# 3. PDF Generator
# -----------------------------------------------------------
def create_pdf(inputs, res, rebar_data, calc_steps):
    pdf = FPDF()
    pdf.add_page()

    pdf.set_font("Arial", 'B', 16)
    pdf.cell(0, 10, "RC Beam Design Report", 0, 1, 'C')
    pdf.set_font("Arial", '', 10)
    pdf.cell(0, 5, "EIT Standard (SDM)", 0, 1, 'C')
    pdf.line(10, 25, 200, 25)
    pdf.ln(10)

    # Summary
    pdf.set_font("Arial", 'B', 12)
    pdf.set_fill_color(230, 230, 230)
    pdf.cell(0, 8, "  1. Design Summary", 0, 1, 'L', fill=True)
    pdf.set_font("Arial", '', 11)

    summary = (
        f"Size: {inputs['b']} x {inputs['h']} cm | Cover: {inputs['cover']} cm\n"
        f"Material: fc' {inputs['fc']} / fy {inputs['fy']} ksc\n"
        f"Main Rebar (Bottom): {rebar_data['bottom_num']}-DB{rebar_data['bottom_db']}\n"
        f"Hanger Rebar (Top): {rebar_data['top_num']}-DB{rebar_data['top_db']}\n"
        f"Stirrup: {res['stirrup_info']}"
    )
    pdf.multi_cell(0, 7, summary)
    pdf.ln(5)

    # Plot
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 8, "  2. Section Detail", 0, 1, 'L', fill=True)
    temp_img = "temp_plot_final.png"
    fig = create_plot(inputs['b'], inputs['h'], inputs['cover'], rebar_data)
    fig.savefig(temp_img, bbox_inches='tight', dpi=100)
    pdf.image(temp_img, x=65, y=None, w=80)
    if os.path.exists(temp_img): os.remove(temp_img)
    pdf.ln(5)

    # Steps
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 8, "  3. Calculation Note", 0, 1, 'L', fill=True)
    pdf.set_font("Arial", '', 10)
    for line in calc_steps:
        pdf.cell(0, 6, line.replace("**", ""), 0, 1)

    return pdf.output(dest='S').encode('latin-1')


# -----------------------------------------------------------
# 4. Streamlit UI (เพิ่ม Form และปุ่มคำนวณ)
# -----------------------------------------------------------
st.set_page_config(page_title="RC Beam Designer", page_icon="🏗️", layout="wide")

st.title("🏗️ RC Beam Design")
st.markdown("---")

# --- Input Form (Sidebar) ---
with st.sidebar.form("design_form"):
    st.header("1. Input Parameters")

    st.subheader("Materials")
    fc = st.number_input("Concrete fc' (ksc)", value=240, step=10)
    fy = st.number_input("Rebar fy (ksc)", value=4000, step=100)

    st.subheader("Geometry")
    col1, col2 = st.columns(2)
    b = col1.number_input("Width (cm)", value=25)
    h = col2.number_input("Depth (cm)", value=50)
    cover = st.slider("Covering (cm)", 2.0, 5.0, 3.0)

    st.subheader("Loads")
    Mu = st.number_input("Moment Mu (kg-m)", value=12000)
    Vu = st.number_input("Shear Vu (kg)", value=8000)

    # ปุ่มคำนวณอยู่ภายใน Form
    submit_btn = st.form_submit_button("🚀 คำนวณออกแบบ (Calculate)")

# --- Main Area ---
if submit_btn:
    # 1. Run Calculation
    inputs = {'fc': fc, 'fy': fy, 'b': b, 'h': h, 'cover': cover, 'Mu': Mu, 'Vu': Vu}
    res, steps = design_beam(Mu, Vu, b, h, cover, fc, fy, 2400)

    if res is None:
        st.error(steps)
    else:
        rebar_data = select_rebar_logic(res['As_req'], b)

        # 2. Display Results
        col_img, col_data = st.columns([1, 1.2])

        with col_img:
            st.subheader("รูปตัดคาน (Section)")
            fig = create_plot(b, h, cover, rebar_data)
            st.pyplot(fig)

        with col_data:
            st.subheader("ผลการออกแบบ (Results)")
            c1, c2 = st.columns(2)
            c1.success(f"**เหล็กล่าง:** {rebar_data['bottom_num']}-DB{rebar_data['bottom_db']}")
            c2.info(f"**เหล็กบน:** {rebar_data['top_num']}-DB{rebar_data['top_db']}")
            st.warning(f"**เหล็กปลอก:** {res['stirrup_info']}")
            st.metric("ปริมาณเหล็กที่ต้องการ (As req)", f"{res['As_req']:.2f} cm²")

            with st.expander("ดูรายการคำนวณ"):
                for s in steps:
                    st.write(s)

            # Download Button
            pdf_bytes = create_pdf(inputs, res, rebar_data, steps)
            st.download_button(
                label="📄 ดาวน์โหลดรายงาน PDF",
                data=pdf_bytes,
                file_name="Beam_Design.pdf",
                mime="application/pdf",
                use_container_width=True
            )
else:
    # หน้าจอตอนแรก (ยังไม่กดปุ่ม)
    st.info("👈 กรุณากรอกข้อมูลทางด้านซ้าย แล้วกดปุ่ม 'คำนวณออกแบบ' เพื่อเริ่มต้น")

    # แสดงรูปตัวอย่างเปล่าๆ ให้ดูสวยงาม
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
        ### วิธีใช้งาน
        1. กรอกค่ากำลังวัสดุ และขนาดหน้าตัด
        2. ใส่ค่าแรง (Moment และ Shear)
        3. กดปุ่ม **คำนวณออกแบบ** ด้านล่างสุด
        """)