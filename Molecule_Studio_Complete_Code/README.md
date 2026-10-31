# Molecule Studio

A mobile-friendly Streamlit chemistry application.

## Features
- SMILES validation
- 2D molecular structure generation and PNG download
- Interactive 3D viewer (stick, sphere, line)
- Molecular properties and descriptors
- PDF molecular report download
- Comparison of up to six molecules and CSV export
- Chemistry quiz with score and explanations

## Run on your computer

1. Install Python 3.10 or newer.
2. Open Terminal in this folder.
3. Create and activate a virtual environment (recommended).
4. Install dependencies:

   ```bash
   python -m pip install -r requirements.txt
   ```

5. Start the app:

   ```bash
   streamlit run app.py
   ```

6. Streamlit will show a local URL, usually http://localhost:8501.

## Deploy to Streamlit Community Cloud

1. Upload `app.py` and `requirements.txt` to a GitHub repository.
2. Sign in to Streamlit Community Cloud.
3. Create an app and select the repository, branch and `app.py`.
4. Deploy.

## Notes
- RDKit provides chemical parsing, descriptors and 2D drawings.
- 3D conformer generation may fail for some structures; try another valid SMILES.
- The app calculates educational cheminformatics descriptors. It is not medical advice.
