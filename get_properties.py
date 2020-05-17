#!/usr/bin/env python

import numpy as np
import pandas as pd

import mopac

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

    # TODO collect all mopac information
    # TODO use pandas to store information, for easy CSV dump

    molecules = read_csv(args.csv)

    col_molidx = "MoleculeIdx"
    col_proidx = "ProtonatedIdx"
    col_refsmi = "ReferenceSmiles"
    col_prosmi = "ProtonatedSmiles"
    col_nenergy = "NeutralEnergy"
    col_energy = "ProtonatedEnergy"

    columns = [col_molidx, col_proidx, col_refsmi, col_prosmi, col_nenergy, col_energy]

    cdata = pd.DataFrame(columns=columns)

    for i, molset in enumerate(molecules):

        filename = args.scr + "/" + f"xyz{i}_n.out"
        try:
            properties = mopac.read_properties(filename)
        except:
            print("error", filename)
            continue

        neutral_energy = properties["h"]

        for j, prosmi in enumerate(molset[1:]):

            filename = args.scr + "/" + f"xyz{i}_{j}.out"
            try:
                properties = mopac.read_properties(filename)
            except:
                print("error",  filename)
                continue

            energy = properties["h"]

            row = {
                col_molidx: i,
                col_proidx: j,
                col_refsmi: molset[0],
                col_prosmi: molset[1+j],
                col_nenergy: neutral_energy,
                col_energy: energy,
            }

            cdata = cdata.append(row, ignore_index=True)

    cdata.to_csv("pm3_properties.csv", index=False)

    return


if __name__ == '__main__':
    main()
