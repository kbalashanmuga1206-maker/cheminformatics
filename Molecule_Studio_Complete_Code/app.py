"""
Molecule Studio — mobile-friendly chemistry explorer
Run locally:
    pip install -r requirements.txt
    streamlit run app.py
"""
from __future__ import annotations

import io
import random
from datetime import datetime

import streamlit as st
import pandas as pd

from rdkit import Chem
from rdkit.Chem import (
    AllChem, Descriptors, Draw, Lipinski, rdMolDescriptors,
    Crippen, QED
)
from rdkit.Chem.Draw import rdMolDraw2D

try:
    import py3Dmol
    HAS_3DMOL = True
except Exception:
    HAS_3DMOL = False

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image as RLImage
)

# --------------------------- Page setup ------------------------------------
st.set_page_config(
    page_title="Molecule Studio",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Poppins:wght@400;500;600;700;800&display=swap');
html, body, [class*="css"] { font-family: 'Poppins', sans-serif; }
.stApp { background: linear-gradient(135deg,#f4f7ff 0%,#f8f4ff 55%,#eefcff 100%); }
.block-container { padding-top: 1.2rem; padding-bottom: 2.5rem; max-width: 1200px; }
.hero {
  padding: 26px 28px; border-radius: 22px; color: white;
  background: linear-gradient(115deg,#4338ca,#7c3aed 55%,#0891b2);
  box-shadow: 0 12px 30px rgba(67,56,202,.18); margin-bottom: 18px;
}
.hero h1 {font-size: clamp(2rem,5vw,3.1rem); font-weight:800; margin:0;}
.hero p {font-size:1.05rem; opacity:.95; margin:.35rem 0 0;}
div[data-testid="stMetric"] {
  background: white; padding: 14px 16px; border-radius: 16px;
  border: 1px solid #e8e9ff; box-shadow: 0 5px 18px rgba(30,41,59,.05);
}
.stTabs [data-baseweb="tab"] {font-size:1rem; font-weight:700;}
.stButton>button, .stDownloadButton>button {
  border-radius: 12px; font-weight:700; min-height:2.8rem;
}
div[data-testid="stSidebar"] {background: #f0efff;}
.small-note {color:#64748b; font-size:.88rem;}
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hero">
  <h1>🧪 Molecule Studio</h1>
  <p>Explore molecular structures, properties and chemistry — all in one place.</p>
</div>
""", unsafe_allow_html=True)

# --------------------------- Helpers ---------------------------------------
EXAMPLES = {
    "Caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Nicotine": "CN1CCC[C@H]1c1cccnc1",
    "Benzene": "c1ccccc1",
    "Ethanol": "CCO",
    "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
}

def mol_from_smiles(smiles: str):
    if not smiles or not smiles.strip():
        return None, "Enter a SMILES notation first."
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None, "Invalid SMILES. Check atom symbols, brackets and bond notation."
    return mol, "Valid SMILES — molecular structure generated successfully."

def property_dict(mol):
    return {
        "Molecular formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight (g/mol)": round(Descriptors.MolWt(mol), 3),
        "Exact mass (Da)": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "TPSA (Å²)": round( rdMolDescriptors.CalcTPSA(mol), 3),
        "H-bond donors": Lipinski.NumHDonors(mol),
        "H-bond acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable bonds": Lipinski.NumRotatableBonds(mol),
        "Heavy atoms": mol.GetNumHeavyAtoms(),
        "Aromatic rings": Lipinski.NumAromaticRings(mol),
        "Fraction Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "QED score": round(QED.qed(mol), 3),
    }

def draw_2d(mol, size=(650, 430)):
    copy = Chem.Mol(mol)
    AllChem.Compute2DCoords(copy)
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.bondLineWidth = 2.2
    drawer.DrawMolecule(copy)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def make_3d_html(mol, style="stick"):
    if not HAS_3DMOL:
        return None
    copy = Chem.AddHs(Chem.Mol(mol))
    status = AllChem.EmbedMolecule(copy, randomSeed=42)
    if status != 0:
        return None
    AllChem.MMFFOptimizeMolecule(copy, maxIters=200)
    block = Chem.MolToMolBlock(copy)
    viewer = py3Dmol.view(width=760, height=460)
    viewer.addModel(block, "mol")
    if style == "sphere":
        viewer.setStyle({"stick": {}, "sphere": {"scale": 0.28}})
    elif style == "line":
        viewer.setStyle({"line": {"linewidth": 3}})
    else:
        viewer.setStyle({"stick": {"radius": 0.16}})
    viewer.setBackgroundColor("#fbfdff")
    viewer.zoomTo()
    return viewer._make_html()

def build_pdf(name, smiles, props, png_bytes=None):
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer, pagesize=A4, rightMargin=18*mm, leftMargin=18*mm,
        topMargin=16*mm, bottomMargin=16*mm
    )
    styles = getSampleStyleSheet()
    title = ParagraphStyle(
        "StudioTitle", parent=styles["Title"], textColor=colors.HexColor("#4338ca"),
        fontName="Helvetica-Bold", fontSize=22, spaceAfter=8
    )
    story = [
        Paragraph("Molecule Studio — Molecular Report", title),
        Paragraph(f"<b>Molecule:</b> {name}", styles["Normal"]),
        Paragraph(f"<b>Generated:</b> {datetime.now().strftime('%d %b %Y, %I:%M %p')}", styles["Normal"]),
        Spacer(1, 8),
        Paragraph("<b>SMILES notation</b>", styles["Heading2"]),
        Paragraph(smiles.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"), styles["Code"]),
        Spacer(1, 8),
    ]
    if png_bytes:
        story.append(RLImage(io.BytesIO(png_bytes), width=125*mm, height=82*mm))
        story.append(Spacer(1, 8))
    story.append(Paragraph("Molecular properties", styles["Heading2"]))
    data = [["Property", "Value"]] + [[str(k), str(v)] for k, v in props.items()]
    table = Table(data, colWidths=[85*mm, 75*mm], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#4338ca")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("GRID", (0,0), (-1,-1), .4, colors.HexColor("#cbd5e1")),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f1f5f9")]),
        ("PADDING", (0,0), (-1,-1), 7),
    ]))
    story.append(table)
    story.append(Spacer(1, 10))
    story.append(Paragraph("Note: Calculated descriptors are educational cheminformatics estimates, not medical or regulatory advice.", styles["Italic"]))
    doc.build(story)
    return buffer.getvalue()

# --------------------------- Sidebar ---------------------------------------
st.sidebar.title("🧭 Navigation")
page = st.sidebar.radio(
    "Choose a workspace",
    ["Molecule Explorer", "Compare Molecules", "Chemistry Quiz", "About"]
)
st.sidebar.markdown("---")
st.sidebar.caption("Molecule Studio • Learning edition")
st.sidebar.caption("Built with Streamlit + RDKit")

# --------------------------- Molecule Explorer -----------------------------
if page == "Molecule Explorer":
    st.subheader("🔬 Molecule Explorer")
    left, right = st.columns([1.15, 1])
    with left:
        example_name = st.selectbox("Load an example molecule", ["Custom SMILES"] + list(EXAMPLES))
        default_smiles = "" if example_name == "Custom SMILES" else EXAMPLES[example_name]
        smiles = st.text_area(
            "Enter SMILES notation",
            value=default_smiles,
            height=95,
            placeholder="Example: CCO",
            help="SMILES is a text representation of a molecular structure."
        )
        molecule_name = st.text_input(
            "Molecule name",
            value=example_name if example_name != "Custom SMILES" else "My molecule"
        )
        validate = st.button("Validate & Explore molecule", type="primary", use_container_width=True)
    with right:
        st.info("**Quick start**\n\nChoose an example or paste a SMILES string, then select **Validate & Explore molecule**.")

    if validate:
        mol, message = mol_from_smiles(smiles)
        st.session_state["active_smiles"] = smiles.strip()
        st.session_state["active_name"] = molecule_name.strip() or "My molecule"
        st.session_state["active_valid"] = mol is not None
        if mol is None:
            st.error(message)
        else:
            st.success(message)
            st.session_state["active_props"] = property_dict(mol)
            st.session_state["active_mol_png"] = draw_2d(mol)

    active_smiles = st.session_state.get("active_smiles", "")
    if active_smiles and st.session_state.get("active_valid", False):
        mol, msg = mol_from_smiles(active_smiles)
        props = property_dict(mol)
        st.markdown("---")
        st.subheader(f"🧬 {st.session_state.get('active_name', 'Molecule')}")
        tab1, tab2, tab3, tab4 = st.tabs(["2D Structure", "Interactive 3D", "Properties", "PDF Report"])

        with tab1:
            st.markdown("**2D molecular structure**")
            png = draw_2d(mol)
            st.image(png, use_container_width=True)
            st.download_button(
                "⬇️ Download 2D structure (PNG)", data=png,
                file_name="molecule_structure.png", mime="image/png"
            )

        with tab2:
            st.markdown("**Rotate, zoom and inspect the 3D molecule**")
            style = st.selectbox("3D visualization style", ["stick", "sphere", "line"])
            html = make_3d_html(mol, style)
            if html:
                import streamlit.components.v1 as components
                components.html(html, height=490, scrolling=False)
            else:
                st.warning("3D viewer is unavailable or this molecule could not be embedded. Install py3Dmol and try again.")

        with tab3:
            st.markdown("**Calculated molecular descriptors**")
            c1, c2, c3 = st.columns(3)
            c1.metric("Molecular weight", f"{props['Molecular weight (g/mol)']} g/mol")
            c2.metric("LogP", props["LogP"])
            c3.metric("TPSA", f"{props['TPSA (Å²)']} Å²")
            df = pd.DataFrame(list(props.items()), columns=["Property", "Value"])
            st.dataframe(df, use_container_width=True, hide_index=True)

        with tab4:
            st.markdown("Create a PDF containing the SMILES, 2D structure and calculated properties.")
            pdf = build_pdf(st.session_state.get("active_name", "Molecule"), active_smiles, props, draw_2d(mol))
            st.download_button(
                "📄 Download molecular PDF report", data=pdf,
                file_name="molecule_report.pdf", mime="application/pdf",
                type="primary", use_container_width=True
            )

# --------------------------- Compare ---------------------------------------
elif page == "Compare Molecules":
    st.subheader("⚖️ Compare Molecules")
    st.write("Enter up to six SMILES strings. Valid molecules are compared using calculated descriptors.")
    count = st.selectbox("Number of molecules", [2, 3, 4, 5, 6], index=0)
    inputs = []
    cols = st.columns(2)
    for i in range(count):
        with cols[i % 2]:
            nm = st.text_input(f"Molecule {i+1} name", value=f"Molecule {i+1}", key=f"cmp_name_{i}")
            sm = st.text_input(
                f"Molecule {i+1} SMILES",
                value=EXAMPLES["Caffeine"] if i == 0 else (EXAMPLES["Aspirin"] if i == 1 else ""),
                key=f"cmp_smiles_{i}",
                placeholder="Paste SMILES"
            )
            inputs.append((nm, sm))
    if st.button("Compare molecules", type="primary", use_container_width=True):
        rows, errors = [], []
        for nm, sm in inputs:
            mol, msg = mol_from_smiles(sm)
            if mol is None:
                errors.append(f"{nm}: {msg}")
            else:
                row = {"Molecule": nm, "SMILES": sm}
                row.update(property_dict(mol))
                rows.append(row)
        for err in errors:
            st.error(err)
        if rows:
            frame = pd.DataFrame(rows).set_index("Molecule")
            st.success(f"Successfully compared {len(rows)} molecule(s).")
            st.dataframe(frame, use_container_width=True)
            st.download_button(
                "⬇️ Download comparison as CSV",
                data=frame.to_csv().encode("utf-8"),
                file_name="molecule_comparison.csv",
                mime="text/csv"
            )

# --------------------------- Quiz ------------------------------------------
elif page == "Chemistry Quiz":
    st.subheader("🎯 Chemistry Quiz")
    QUESTIONS = [
        {"q": "What does SMILES represent?", "opts": ["A text representation of a molecule", "A protein sequence", "A lab temperature", "A type of microscope"], "a": 0, "why": "SMILES encodes molecular structures as text."},
        {"q": "Which element has the symbol O?", "opts": ["Gold", "Oxygen", "Osmium", "Hydrogen"], "a": 1, "why": "O is the chemical symbol for oxygen."},
        {"q": "What does molecular weight commonly express?", "opts": ["Colour", "Mass per mole", "Boiling point", "pH"], "a": 1, "why": "Molar mass is expressed in grams per mole (g/mol)."},
        {"q": "Which bond involves sharing electron pairs?", "opts": ["Ionic bond", "Covalent bond", "Metallic lattice only", "Hydrogen bond only"], "a": 1, "why": "Covalent bonds form when atoms share electron pairs."},
        {"q": "What is the molecular formula of water?", "opts": ["CO₂", "H₂O", "O₂", "H₂"], "a": 1, "why": "A water molecule contains two hydrogen atoms and one oxygen atom."},
        {"q": "Which descriptor estimates lipophilicity?", "opts": ["LogP", "TPSA", "H-bond donors", "Atom count"], "a": 0, "why": "LogP is commonly used as an estimate of partitioning between octanol and water."},
        {"q": "What does TPSA stand for?", "opts": ["Total Polar Surface Area", "Temperature Pressure Scale Amount", "Total Protein Sequence Analysis", "Titration Point Solution Average"], "a": 0, "why": "TPSA means topological polar surface area."},
    ]
    if "quiz_order" not in st.session_state:
        st.session_state.quiz_order = random.sample(range(len(QUESTIONS)), 5)
        st.session_state.quiz_submitted = False
    with st.form("quiz_form"):
        answers = {}
        for pos, qidx in enumerate(st.session_state.quiz_order, start=1):
            q = QUESTIONS[qidx]
            answers[qidx] = st.radio(f"{pos}. {q['q']}", q["opts"], key=f"quiz_{qidx}")
        submitted = st.form_submit_button("Submit answers", type="primary", use_container_width=True)
    if submitted:
        score = 0
        for qidx in st.session_state.quiz_order:
            q = QUESTIONS[qidx]
            selected = answers[qidx]
            correct = q["opts"][q["a"]]
            if selected == correct:
                score += 1
        st.session_state.quiz_score = score
        st.session_state.quiz_answers = answers
        st.session_state.quiz_submitted = True
    if st.session_state.get("quiz_submitted"):
        score = st.session_state.quiz_score
        st.success(f"Your score: {score}/5")
        for qidx in st.session_state.quiz_order:
            q = QUESTIONS[qidx]
            selected = st.session_state.quiz_answers[qidx]
            correct = q["opts"][q["a"]]
            if selected == correct:
                st.markdown(f"✅ **{q['q']}** — Correct")
            else:
                st.markdown(f"❌ **{q['q']}** — Correct answer: **{correct}**")
            st.caption(q["why"])
        if st.button("Try another quiz"):
            st.session_state.quiz_order = random.sample(range(len(QUESTIONS)), 5)
            st.session_state.quiz_submitted = False
            st.rerun()

# --------------------------- About -----------------------------------------
else:
    st.subheader("ℹ️ About Molecule Studio")
    st.markdown("""
    Molecule Studio is an educational cheminformatics application for exploring molecules from SMILES notation.

    **Main features**
    - SMILES validation and molecule generation
    - 2D molecular structure drawing
    - Interactive 3D viewer with stick, sphere and line styles
    - Molecular descriptors and properties
    - Downloadable PDF molecular reports
    - Side-by-side comparison of up to six molecules
    - Chemistry quiz and instant feedback

    **Technology:** Python, Streamlit, RDKit, py3Dmol, Pandas and ReportLab.

    **Developer:** K. Balashanmuga  
    **Programme:** Artificial Intelligence and Machine Learning (AIML)  
    **Institution:** SRM Institute of Science and Technology, Tiruchirappalli
    """)
    st.info("For learning and academic demonstration. Descriptor values are computational estimates.")

st.markdown("---")
st.markdown(
    "<div style='text-align:center;color:#64748b;font-size:.85rem'>"
    "© 2026 Molecule Studio • Educational Chemistry Explorer</div>",
    unsafe_allow_html=True
)
