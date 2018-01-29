#! /usr/bin/env python
"""This script prompts a user to enter a PDB filename
(and optional padding distance) to calculate min/max and
box dimensions for MD simulation setup"""
import argparse
import mdtraj as md
import numpy as np

def main():
    """Main function
    """

    # Argparse setup
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="Name of the input PDB file")
    parser.add_argument("-p", "--padding", default=12, type=float,
                        help="Distance (in angstroms) to pad the"
                        " protein by when calculating box "
                        "dimensions (default value : 12)")
    args = parser.parse_args()

    # Specify PDB file and size of box
    t = md.load(args.pdb)
    padding_dist = args.padding #angstroms

    #Initial protein dimensions calculations
    max_coords = np.amax(t.xyz[0], axis=0)*10 #Multiplied by 10 to convert into angstroms
    min_coords = np.amin(t.xyz[0], axis=0)*10 #Multiplied by 10 to convert into angstroms
    box_dims = np.ceil(max_coords - min_coords)+(padding_dist*2)
    cellOrigin = min_coords+((max_coords-min_coords)/2)

    # Adjusted box dims and min/max calculations
    max_coords_new = cellOrigin + (box_dims/2)
    min_coords_new = cellOrigin - (box_dims/2)
    box_dims_new = max_coords_new-min_coords_new
    cellOrigin_new = min_coords_new+(box_dims_new/2)
    print('\nBox dims and min/max output:')
    print('Min: '+' '.join(map(str, min_coords_new)))
    print('Max: '+' '.join(map(str, max_coords_new)))
    print('CellOrigin: '+' '.join(map(str, cellOrigin_new)))
    print('Box Dims: '+' '.join(map(str, box_dims_new)) + '\n')

    # NAMD-ready formatting
    print('NAMD-ready formatting:')
    for i in np.arange(3):
        print('cellBasisVector{} {}'.format(i+1, ' '.join(map(str, np.diag(box_dims_new)[i]))))
    print('cellOrigin '+' '.join(map(str, cellOrigin_new)))

if __name__ == "__main__":
    main()
    