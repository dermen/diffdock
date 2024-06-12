#!/usr/bin/env python

"""
This script takes a set of raw PDBs and:

    1. "fixes" them
        a. no crystallization reagents
        b. canonical order
        c. no hydrogens

    2. separates protein and LoI

"""

import os
import argparse
from glob import glob

from os.path import join as pjoin

from wharf import pdbutils


def main(glob_string, output_directory):

    print('')
    needed_dirs = ['fixed', 'protein', 'ligand']
    for nd in needed_dirs:
        nd_path = pjoin(output_directory, nd)
        print(nd_path)
        if not os.path.exists(nd_path):
            os.makedirs(nd_path)

    print('results --> %s' % output_directory)
    print('directories correct')

    input_pdbs = glob(glob_string)
    n_input = len(input_pdbs)
    print('Processing %d pdbs...' % n_input)


    for i,pdb in enumerate(input_pdbs):

        bn = os.path.basename(pdb)
        print(' - %04d / %04d %s' % (i+1, n_input, bn))

        with open(pdb, 'r') as f:
            pdb_txt = f.read()


        # 1. fix pdbs

        f_pdb_txt = pdbutils.remove_hydrogens(pdb_txt)
        f_pdb_txt = pdbutils.remove_cryst_reagents(f_pdb_txt)
        f_pdb_txt = pdbutils.enforce_cannonical_pdb_order(f_pdb_txt)
        
        fixed_pdb_path = pjoin(output_directory, 'fixed', bn)
        with open(fixed_pdb_path, 'w') as f:
            f.write(f_pdb_txt)


        # 2. split ligand & protein
        
        split = pdbutils.split_ligands_protein(f_pdb_txt)
        for chain_id in split.keys():
            if len(split[chain_id]) < 2:
                print(f'no ligand found in: {bn}, chain {chain_id}, skipping')
                continue
            elif len(split[chain_id]) > 2:
                print(f'more than one ligand found in: {bn}, chain {chain_id}, skipping')
                continue

            out_bn = os.path.splitext(bn)
            out_name = f'{out_bn[0]}_{chain_id}{out_bn[1]}'
            protein_path = pjoin(output_directory, 'protein', out_name)
            ligand_path  = pjoin(output_directory, 'ligand',  out_name)
            with open(protein_path, 'w') as f:
                f.write(split[chain_id][0])
            with open(ligand_path, 'w') as f:
                f.write(split[chain_id][1])

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get PDBs ready for docking & ML')
    parser.add_argument('pdb_glob',
                        help='a glob string to match PDB files to process')
    parser.add_argument('-o', '--outdir', default='.',
                        help='the directory to write the results')
    args = parser.parse_args()

    main(args.pdb_glob, args.outdir)
