# Molecular_Dynamics_Simulation
A GROMACS-based workflow to analyze binding and stability of a small molecule with a mutant protein kinase, including system setup, MD simulation, and trajectory analysis.
# Molecular Dynamics Simulations Workflow for Small Molecule Binding to Protein Kinase Mutant

**Project Type:** Research / Computational Biochemistry  
**Objective:** To study the binding interactions and stability of a small-molecule ligand with a mutant protein kinase using Molecular Dynamics simulations in GROMACS.


# -------------- Tools & Software --------------
- GROMACS (manual or apt install)
- Chimera & PyMOL (visualization & ligand prep)
- AutoDock Vina (optional docking)
- VMD (trajectory visualization)
- Grace (graphing)


# -------------- System Setup -------------- 
1. **Prepare Receptor and Ligand**
   - Chimera used to process the target protein structure and the small-molecule ligand.
   - Ligand optimized and converted to proper formats for GROMACS.

2. **Topology Generation**
   - Protein: CHARMM27 force field
   - Ligand: SwissParam-generated topology
   - Position restraints applied to ligand

3. **Solvation and Ion Addition**
   - Solvent: TIP3P water
   - Ion concentration: 0.1 M, neutral system


# -------------- Simulation Workflow (GROMACS) -------------- 
- **Energy Minimization**: EM.mdp → `gmx mdrun -deffnm EM`
- **Equilibration**: NVT → NPT
- **Production MD**: `MD.mdp`, time adjustable
- **Trajectory Postprocessing**:
  - Re-centering & rewrapping: `gmx trjconv`
  - RMSD, RMSF, Hydrogen Bonds, Radius of Gyration, SASA
  - Plots with `xmgrace`


# -------------- How to Run -------------- 
For detailed Molecular Dynamics simulation commands and step-by-step instructions, see [MD_Simulation_Steps.md](MD_Simulation_Steps.md)
