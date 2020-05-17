#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy
import sklearn
import copy

import rmsd

import rdkit.Chem as Chem
import rdkit.Chem.rdmolops as rdmolops
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.ChemicalForceFields as ChemicalForceFields
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors


def molobj_to_coordinates(molobj, idx=-1):
    """
    """

    conformer = molobj.GetConformer(id=int(idx))
    coordinates = conformer.GetPositions()
    coordinates = np.array(coordinates)

    return coordinates


def molobj_to_atoms(molobj, atom_type=int):

    atoms = molobj.GetAtoms()

    if atom_type == str:
        atoms = [atom.GetSymbol() for atom in atoms]

    elif atom_type == int:
        atoms = [atom.GetAtomicNum() for atom in atoms]
        atoms = np.array(atoms)

    return atoms


def molobj_to_axyzc(molobj, atom_type=int, idx=-1):
    """
    rdkit molobj to xyz
    """

    atoms = molobj_to_atoms(molobj, atom_type=atom_type)
    coordinates = molobj_to_coordinates(molobj, idx=idx)
    charge = rdmolops.GetFormalCharge(molobj)

    return atoms, coordinates, charge


def genereate_conformers(molsmi,
    max_conf=20,
    min_conf=10,
    max_steps=1000):
    """
    """

    molobj = Chem.MolFromSmiles(molsmi)

    if molobj is None:
        return None

    molobj = Chem.AddHs(molobj)

    status_embed = AllChem.EmbedMolecule(molobj)

    if status_embed != 0:
        return None

    status_optim = AllChem.UFFOptimizeMolecule(molobj, maxIters=max_steps)

    # Keep unconverged uff
    # if status_optim != 0:
    #     return None

    # Check bond lengths
    dist = Chem.rdmolops.Get3DDistanceMatrix(molobj)
    np.fill_diagonal(dist, 10.0)
    min_dist = np.min(dist)

    # For some atom_types in UFF, it will fail
    if min_dist < 0.001:
        print("fail", smilesstr)
        return None

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(molobj)

    confs = min(1 + 3*rot_bond, max_conf)
    confs = max(confs, min_conf)

    status = AllChem.EmbedMultipleConfs(molobj,
        numConfs=confs,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True)

    return molobj


def select_conformer(molobj, method=None):
    """

    Find the conformer with lowest energy

    # TODO Add other energy methods than MMFF
    # TODO Add interface to ppqm.sqm

    """

    res = AllChem.MMFFOptimizeMoleculeConfs(molobj, numThreads=0)
    res = np.array(res)
    energies = res[:,1]
    not_converged = res[:,1]

    # TODO Ignore not_converged != 0

    idx = np.argmin(energies)

    coord = molobj_to_coordinates(molobj, idx=idx)

    conf = Chem.Conformer(len(coord))
    for i, xyz in enumerate(coord):
        conf.SetAtomPosition(i, xyz)

    molobj = copy.deepcopy(molobj)
    molobj.RemoveAllConformers()
    molobj.AddConformer(conf, assignId=True)

    return molobj


def find_charged_atom(molobj):

    atoms = list(molobj.GetAtoms())

    charges = [atom.GetFormalCharge() for atom in atoms]
    charges = np.array(charges)

    idxs, = np.where(charges == 1)

    if len(idxs) == 0:
        return None, None

    idx = idxs[0]

    atom = atoms[idx]
    neighbors = list(atom.GetNeighbors())
    neighbors_types = [a.GetAtomicNum() for a in neighbors]
    neighbors_types = np.array(neighbors_types)
    h_idxs, = np.where(neighbors_types == 1)
    h_idxs = [neighbors[i].GetIdx() for i in h_idxs]

    return idx, h_idxs


def expand_graph(smi):
    """
    Generate 3D coordinates and select most stable conformer

    :molset: TODO
    :returns: XYZ STR

    """

    # Generate conformers
    molobj = genereate_conformers(smi)

    if molobj is None:
        return None

    # Select most stable conformer
    molobj = select_conformer(molobj)

    # Find charged atom for print
    idx, h_idxs = find_charged_atom(molobj)
    if h_idxs is not None:
        h_idxs = [str(x) for x in h_idxs]
        h_idxs = ",".join(h_idxs)

    atoms, xyz, charge = molobj_to_axyzc(molobj, atom_type=str)
    xyzstr = rmsd.set_coordinates(atoms, xyz, title=f" title charge={charge} highlight={idx} hydrogens={h_idxs}")

    return xyzstr


def dump_xyz(molset, name="_tmp_"):
    """TODO: Docstring for dump_xyz.

    :molset: TODO
    :returns: TODO

    """

    neu_smi = molset[0]
    pro_smis = molset[1:]

    xyzstr = expand_graph(neu_smi)
    if xyzstr is None:
        print("error", molset[0])
        return

    with open(name + "_n.xyz", 'w') as f:
        f.write(xyzstr)

    for i, smi in enumerate(pro_smis):
        xyzstr = expand_graph(smi)
        if xyzstr is None:
            print("error", smi)
            continue

        with open(name + f"_{i}.xyz", 'w') as f:
            f.write(xyzstr)

    return


def read_csv(filename, sep=", "):

    with open(filename, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        line = line.strip()
        lines[i] = line.split(sep)

    return lines


def main(args=None):

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version="1.0")
    parser.add_argument('--debug', action='store_true', help='')
    parser.add_argument('--csv', action='store', help='', metavar='FILE')
    parser.add_argument('--scr', action='store', help='', metavar='DIR')

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    molecules = read_csv(args.csv)

    for i, molset in enumerate(molecules):

        print(molset[0])

        name = args.scr + f"{i}"

        dump_xyz(molset, name=name)


    return


if __name__ == '__main__':
    main()


