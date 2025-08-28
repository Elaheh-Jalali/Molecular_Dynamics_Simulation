# Detailed Molecular Dynamics Simulation Steps for Small Molecule Binding to a Protein Mutant

This file provides step-by-step instructions to perform Molecular Dynamics simulations of a small-molecule ligand binding to a protein mutant using GROMACS. Follow these steps after preparing your system as described in the main README.


# -------------- 1. Update Ubuntu and install build tools --------------
sudo apt update
sudo apt upgrade -y
sudo apt install -y gcc cmake build-essential libfftw3-dev
# This installs essential compilers and libraries needed for GROMACS.


# -------------- 2. Software Installation --------------
# Install GROMACS (quick install)
sudo apt install -y gromacs

# Install PyMOL for visualization
sudo apt-get install -y pymol

# Install Chimera
# Download installer from https://www.cgl.ucsf.edu/chimera/download.html
chmod +x CHIMERA-INSTALLER.bin
./CHIMERA-INSTALLER.bin

# Install AutoDock Vina
# Download from https://vina.scripps.edu/downloads/
tar -xzvf autodock_vina_1_1_2_linux_x86.tgz

# Install GRACE graphing tool
sudo apt-get install -y grace

# Install VMD
# Download .tar.gz from https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD
./configure
cd src
sudo make install
vmd

# Manual GROMACS compilation
# Download from https://manual.gromacs.org/current/download.html
tar xfz gromacs-2025.2.tar.gz
cd gromacs-2025.2


# -------------- 9. Ligand and receptor preparation in Chimera --------------
# Open Protein.pdb, delete protein chains, keep ligand
# Add hydrogens to ligand and save as LIG.mol2
# Ensure the first line of LIG.mol2 file is "@<TRIPOS>MOLECULE" and update the molecule name to "LIG"
perl sort_mol2_bonds.pl LIG.mol2 LIG.mol2
# Upload to SwissParam (http://www.swissparam.ch/) and download topology
# Prepare receptor REC.pdb with DockPrep


# -------------- 10. Load GROMACS environment (if manually compiled) --------------
source /usr/local/gromacs/bin/GMXRC


# -------------- 11. Generate protein topology --------------
gmx pdb2gmx -f REC.pdb -ignh
# Choose force field (e.g., CHARMM27) and water model (e.g., TIP3P)


# -------------- 12. Prepare ligand coordinates --------------
gmx editconf -f LIG.pdb -o LIG.gro
# Copy lines 3 to second last line from LIG.gro to conf.gro
# Adjust atom count in conf.gro
# Visual check in Chimera


# -------------- 13. Edit topol.top --------------
gedit topol.top
# Add 
; Include ligand topology 
#include "LIG.itp"

below- Include forcefield parameters
#include "amberGS.ff/forcefield.itp"

# Add 
LIG 1
below, Protein_chain_E 1


# -------------- 14. Define box and solvate --------------
gmx editconf -f conf.gro -d 1.0 -bt triclinic -o box.gro
gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o box_sol.gro


# -------------- 15. Prepare ions input --------------
gmx grompp -f ions.mdp -c box_sol.gro -p topol.top -o ION.tpr


# -------------- 16. Add ions and neutralize system --------------
gmx genion -s ION.tpr -p topol.top -conc 0.1 -neutral -o box_sol_ion.gro
# Choose solvent group (usually 15)


# -------------- 17. Energy minimization --------------
gmx grompp -f EM.mdp -c box_sol_ion.gro -p topol.top -o EM.tpr
gmx mdrun -v -deffnm EM


# -------------- 18. Create ligand index and position restraints --------------
gmx make_ndx -f LIG.gro -o index_LIG.ndx
# Enter: 0 & ! a H* (select all except hydrogens), then q
gmx genrestr -f LIG.gro -n index_LIG.ndx -o posre_LIG.itp -fc 1000 1000 1000
# Select ligand group number (e.g., 3)


# -------------- 19. Edit topol.top to include ligand position restraints --------------
gedit topol.top
# Add the following ligand position restraints block at the end of the file, 
# immediately after the protein position restraints section:

        ; Ligand position restraints
		#ifdef POSRES
		#include "posre_LIG.itp"
		#endif

# The protein position restraints section looks like this:

		"; Include Position restraint file
		#ifdef POSRES
		#include "posre.itp"
		#endif



# -------------- 20. Create system index file --------------
gmx make_ndx -f EM.gro -o index.ndx
# Enter: 1 | 13 (select protein and ligand), then q


# -------------- 21. NVT equilibration --------------
gmx grompp -f NVT.mdp -c EM.gro -r EM.gro -p topol.top -n index.ndx -maxwarn 2 -o NVT.tpr
gmx mdrun -deffnm NVT


# -------------- 22. NPT equilibration --------------
gmx grompp -f NPT.mdp -c NVT.gro -r NVT.gro -p topol.top -n index.ndx -maxwarn 2 -o NPT.tpr
gmx mdrun -deffnm NPT


# -------------- 23. Production MD run --------------
gmx grompp -f MD.mdp -c NPT.gro -t NPT.cpt -p topol.top -n index.ndx -maxwarn 2 -o MD.tpr
gmx mdrun -deffnm MD
