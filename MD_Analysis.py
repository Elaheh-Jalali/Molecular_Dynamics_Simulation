"""
# Molecular Dynamics Analysis for Protein-Ligand Complex Includes:
- Re-centering & rewrapping trajectory
- Extract first frame
- RMSD & RMSF
- Radius of gyration
- Solvent Accessible Surface Area (SASA)
- Hydrogen bonds analysis
All results are plotted and saved as PNG files.

"""

import MDAnalysis as mda
from MDAnalysis.analysis import rms, hbonds, sas
import matplotlib.pyplot as plt
import numpy as np

# -------------- Input Files -------------- 
topology = "topol.tpr"       # GROMACS topology file
trajectory = "MD_center.xtc" # GROMACS trajectory
u = mda.Universe(topology, trajectory)

# -------------- Atom Selections -------------- 
protein = u.select_atoms("protein")       # protein selection
ligand = u.select_atoms("resname LIG")    # ligand selection

# -------------- Recenter & Rewrap Trajectory -------------- 
for ts in u.trajectory:
    ts.positions -= protein.center_of_mass()
    u.atoms.wrap(compound="fragments")

# -------------- Extract First Frame -------------- 
u.select_atoms("all").write("first_frame.pdb")

# -------------- RMSD -------------- 
rmsd_analysis = rms.RMSD(protein).run()
time = rmsd_analysis.rmsd[:,1]       # time (ps)
rmsd_values = rmsd_analysis.rmsd[:,2] # RMSD (nm)

plt.figure()
plt.plot(time, rmsd_values)
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (nm)")
plt.title("Protein RMSD")
plt.savefig("RMSD_protein.png", dpi=300)
plt.close()

# -------------- RMSF -------------- 
rmsf_analysis = rms.RMSF(protein).run()
residues = np.arange(len(rmsf_analysis.rmsf))

plt.figure()
plt.plot(residues, rmsf_analysis.rmsf)
plt.xlabel("Residue")
plt.ylabel("RMSF (nm)")
plt.title("Protein RMSF")
plt.savefig("RMSF_protein.png", dpi=300)
plt.close()

# -------------- Radius of Gyration -------------- 
rg = protein.radius_of_gyration()
print(f"Radius of gyration (final frame): {rg:.3f} nm")

# -------------- SASA -------------- 
sasa_analysis = sas.SASA(protein).run()
plt.figure()
plt.plot(sasa_analysis.times, sasa_analysis.results.sasa)
plt.xlabel("Time (ps)")
plt.ylabel("SASA (nm^2)")
plt.title("Protein SASA")
plt.savefig("SASA_protein.png", dpi=300)
plt.close()

# -------------- Hydrogen Bonds -------------- 
hbond_analysis = hbonds.HydrogenBondAnalysis(u, selection1="protein", selection2="resname LIG")
hbond_analysis.run()
times = hbond_analysis.times
hb_counts = [len(frame) for frame in hbond_analysis.hbonds]

plt.figure()
plt.plot(times, hb_counts)
plt.xlabel("Time (ps)")
plt.ylabel("Number of H-bonds")
plt.title("Protein-Ligand Hydrogen Bonds")
plt.savefig("HBonds.png", dpi=300)
plt.close()

# -------------- Summary -------------- 
print("MD analysis completed. Plots and first frame saved.")
print(f"Final Radius of Gyration: {rg:.3f} nm")
print(f"Number of H-bonds in last frame: {hb_counts[-1] if hb_counts else 0}")
