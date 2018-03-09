import os as os
import numpy as np
import mdtraj as md
import argparse as ap
from string import ascii_uppercase
from mdtraj.formats import PDBTrajectoryFile

# Methods

# Set up parser
parser = ap.ArgumentParser(description='Rename HIS to HIE, HID, or HIP depending on protonation state')
parser.add_argument('pdb', type=str, nargs=1, help='full path to PDB file')
parser.add_argument('-o', type=str, metavar='prefix', nargs=1, default=None, help='prefix for output PDB')
parser.add_argument('-s', '-split', action='store_true') # Split by chain into separate PDBs
args = parser.parse_args()

pdbfile = args.pdb[0]
if args.o is None:
	prefix = os.path.splitext(args.pdb[0])[0] + '_renamed'
else:
	prefix = args.o[0]
split = args.s

# Process PDBfile
pdb = md.load_pdb(pdbfile)
for residue in pdb.topology.residues:
	is_HIS = (residue.name == 'HIS') or (residue.name == 'HIE') or (residue.name == 'HID') or (residue.name == 'HIP') or (residue.name == 'HSE') or (residue.name == 'HSD') or (residue.name == 'HSP')
	if is_HIS:
		names = [x.name for x in residue.atoms if x.element.symbol == 'H']
		is_empty = (len(names) == 0)
		has_HD1 = (names.count('HD1') == 1)
		has_HE2 = (names.count('HE2') == 1)
		if (is_empty):
			print('WARNING: no hydrogens detected for %s, residue name unchanged' % str(residue))
		else:
			if has_HD1 and has_HE2:
				residue.name = 'HIP'
			if has_HD1 and not has_HE2:
				residue.name = 'HID'
			if not has_HD1 and has_HE2:
				residue.name = 'HIE'
			if not has_HD1 and not has_HE2:
				print('WARNING: no hydrogens detected on ND and NE for %s, residue name unchanged' % str(residue))

if split:
	chains = [np.asarray([atom.index for atom in chain.atoms]) for chain in pdb.top.chains]
	letters = [ascii_uppercase[x.index] for x in pdb.top.chains]
	for i, (sel, letter) in enumerate(zip(chains, letters)):
		with PDBTrajectoryFile(prefix + '_chain%s.pdb' % letter, mode='w') as pdbout:
			subset = pdb.atom_slice(sel)
			pdbout.set_chain_names([letter])
			pdbout.write(subset.xyz[0] * 10., subset.topology)
else:
	pdb.save_pdb(prefix + '.pdb')

print('HIS renaming complete!')
