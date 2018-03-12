#! /usr/bin/env python
"""This script prompts a user to enter a PDB filename (and optional arguments 
for custom atom selections and spring constants) to export constraint PDBs."""
import argparse
import mdtraj as md
import numpy as np


def main():
    """Main function
    """

    # Argparse setup
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="Name of the input PDB file")
    parser.add_argument("-s", "--selection", default="protein",
                        help="Provide specific mdtraj-style atom selections"
                        " enclosed in quotes that need to be constrained "
                        "if you want more than the default "
                        "cons_bb/cons_protein PDBs")
    parser.add_argument("-c", "--springconst", default=10,
                        help="Provide a specific spring constant (default: 10)")
    args = parser.parse_args()
    t = md.load(args.pdb)
    custom_selection = args.selection
    spring_constant = args.springconst

    if custom_selection == 'protein':
        cons_bb_array = np.zeros(t.n_atoms)
        cons_protein_array = np.zeros(t.n_atoms)
        cons_protein_array[t.top.select('protein')] = spring_constant
        cons_bb_array[t.top.select('backbone and protein')] = spring_constant
        t.save_pdb('cons_protein.pdb', bfactors=cons_protein_array)
        t.save_pdb('cons_bb.pdb', bfactors=cons_bb_array)
    else:
        cons_array = np.zeros(t.n_atoms)
        cons_array[t.top.select(custom_selection)] = spring_constant
        t.save_pdb('cons.pdb', bfactors=cons_array)


if __name__ == "__main__":
    main()
