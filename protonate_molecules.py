#!/usr/bin/env python

import copy
import numpy as np
import pandas as pd
import scipy
import sklearn
import rdkit.Chem as Chem


def protonate_molobj(molobj, atom):
    """Find all atom and protonate, return to protonated states

    :molobj: TODO
    :atom: TODO
    :returns: TODO

    """

    atoms = molobj.GetAtoms()
    atoms = [x.GetAtomicNum() for x in atoms]
    atoms = np.asarray(atoms, dtype=int)
    idxs, = np.where(atoms == atom)

    # If no 'atom' in 'atoms', then return none
    if len(idxs) == 0: return

    # Protonate atom with idx
    for idx in idxs:
        idx = int(idx)
        molobj_prime = copy.deepcopy(molobj)
        molobj_prime.GetAtomWithIdx(idx).SetFormalCharge(1)

        yield molobj_prime

    return


def main(args=None):

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version="1.0")
    parser.add_argument('--smiles', action='store', help='', metavar='FILE')

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    # Parse GDB file
    sep = "\t"
    df = pd.read_csv(args.smiles, sep=sep, header=None)

    molecules = df.iloc[:,0].to_list()

    target_atom = 7

    for molsmi in molecules:
        molobj = Chem.MolFromSmiles(molsmi)
        molobjs = protonate_molobj(molobj, target_atom)
        molobjs = list(molobjs)
        if len(molobjs) == 0: continue
        protonated = [Chem.MolToSmiles(x) for x in molobjs]

        rtn = [molsmi] + protonated

        print(", ".join(rtn))


    return


if __name__ == '__main__':
    main()
