#!/bin/bash

# 1. Prepare the Protein Topology
# Extract the ligand from the protein-ligand complex
grep UNK1 IUS0_clean.pdb > unk1.pdb

# Extract files containing the force field parameters
tar -zxvf charmm36-jul2022.ff.tgz

# Running the force field on the aldose reductase protein without the ligand
gmx pdb2gmx -f 1US0_clean.pdb -o 3HTB_processed.gro -ter

# 2. Prepare the Ligand Topology
# Convert your ligand to .mol2 file using the Avagrodro's tool. 
# Using the same tool add hygrogen atoms to the molecules
# Sort ligand coordinates to match bond coordinates using the "perl sort_mol2_bonds.pl" script 
# The "perl sort_mol2_bonds.pl" script is located in the script folder in this repository
perl sort_mol2_bonds.pl unk1.mol2 unk1_fix.mol2

# Generate ligand topology with CGenFF
# Visit CGenFF server and upload fixed ligand to generate stream file (.str)
# Convert CHARMM topology to GROMACS format using the cgenff_charmm2gmx_py3_nx2.py  script
# The "cgenff_charmm2gmx_py3_nx2.py" is located in the scripts folder
python cgenff_charmm2gmx_py3_nx2.py UNK1 unk1_fix.mol2 unk1.str charmm36-jul2022.ff


# 3. Build the Protein-Ligand Complex
# Convert ligand to GROMACS format 
gmx editconf -f unk1_ini.pdb -o unk1.gro
# Paste ligand coordinates into protein file to create complex.gro
# Incorporate ligand topology into protein topology


# 4. Defining the Unit Cell & Adding Solvent
# Define unit cell and fill with water
gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

# Adding ions
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# 5. Energy Minimization
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# 6. Equilibration
# Apply restraints to the ligand
gmx make_ndx -f unk1.gro -o index_unk1.ndx
# Execute genrestr module to create position restraints
# Define temperature coupling groups
gmx make_ndx -f em.gro -o index.ndx
# Create temperature coupling groups
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt

# Equilibration, Part 2
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt

# 7. Production of Molecular Dynamics (MD) files
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
gmx mdrun -deffnm md_0_10
