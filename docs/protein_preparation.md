# Protein Preparation Workflow for MD Simulations (`prepare_protein.py`)

This script automates the preparation of protein structures for Molecular Dynamics (MD) simulations using Amber, Propka, Reduce, and AutoDockTools. The preparation pipeline includes the following key steps:

---

## 1. Initial PDB Preparation and Renumbering
- Loads the input PDB file.
- Saves a clean version (`image.pdb`) with:
  - Renumbered residues.
  - Standardized residue names.

---

## 2. Hydrogen Addition and Histidine Protonation (using Reduce)
- Runs the **Reduce** program to:
  - Add hydrogen atoms.
  - Determine protonation states of **Histidine (HIS)** residues.
    - Assigns `HIP`, `HID`, or `HIE` based on presence of `HD1` and/or `HE2` hydrogens.

---

## 3. Ionizable Residue Protonation State Determination (using Propka)
- Executes **Propka** to predict pKa values for all ionizable residues:
  - HIS, ASP, GLU, LYS, N-terminus, C-terminus.
- Based on pKa values and specified pH (default = **7.0**), assigns Amber-specific protonation states:
  - `ASH` for protonated Aspartate.
  - `GLH` for protonated Glutamate.
  - `LYN` for deprotonated Lysine.

---

## 4. Disulfide Bond Detection
- Identifies **Cysteine (CYS)** pairs within **3.0 Å** to form disulfide bonds.
- Renames bonded CYS residues as `CYX`.
- Records the specific atoms involved in the disulfide bridges.

---

## 5. Generating TLeap-ready PDB
- Applies all determined modifications:
  - Histidine variants (`HIP`, `HID`, `HIE`)
  - Protonated/deprotonated forms (from Propka)
  - Disulfide-linked residues (`CYX`)
- Outputs a new PDB: `tleap.pdb`, ready for Amber `tleap`.

---

## 6. Topology and Coordinate Generation (using TLeap)
- Dynamically creates a `tleap.in` input file with:
  - Protein structure (`tleap.pdb`)
  - Amber force fields (`ff14SB`, `TIP3P`)
  - Commands for defining disulfide bonds.
- Runs **tleap** to generate:
  - Topology file: `.parm7`
  - Coordinate file: `.crd`

---

## 7. PQR and PDB File Generation (using ambpdb)
- Uses **ambpdb** to convert Amber outputs into:
  - **PQR file** (`protein.pqr`): includes atomic charges and radii.
  - **PDB file** (`protein.pdb`): standard format for visualization and compatibility.

---

## 8. PDBQT File Generation (using AutoDockTools)
- Optionally runs **`prepare_receptor4.py`** to generate:
  - **PDBQT file** (`protein.pdbqt`) for use with AutoDock Vina.

---