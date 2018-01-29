#! /usr/bin/env python
"""This script prompts a user to enter a PDB filename (and optional integers
corresponding to chains) to export constraint PDBs for all protein and
backbone-only atoms."""
import argparse
import mdtraj as md
import numpy as np

def main():
    """Main function
    """

    # Argparse setup
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="Name of the input PDB file")
    parser.add_argument("-c", "--chains", default='protein',
                        help="Provide integer values corresponding"
                        " to chains of interest that need to be constrained "
                        "(default selection: protein)")
    args = parser.parse_args()

    t = md.load(args.pdb)
    chains_of_interest = args.chains.split(',')
    # Pay attention when specifying chains of non-proteins in the cons_bb section

    cons_bb_array = np.zeros(t.n_atoms)
    cons_protein_array = np.zeros(t.n_atoms)

    if chains_of_interest == ['protein']:
        cons_protein_array[t.top.select('protein')] = 10
        cons_bb_array[t.top.select('backbone and protein')] = 10
    else:
        cons_protein_array[t.top.select(
            'chainid '+' or chainid '.join(str(x) for x in chains_of_interest))] = 10
        cons_bb_array[t.top.select(
            'backbone and (chainid '+' or chainid '.join(str(x) for x in chains_of_interest)+')')] = 10

    t.save_pdb('cons_protein.pdb', bfactors=cons_protein_array)
    t.save_pdb('cons_bb.pdb', bfactors=cons_bb_array)

if __name__ == "__main__":
    main()
