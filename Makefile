
# Protonate
data/gdb7_protonated.smi:
	python protonate_molecules.py --smiles gdb11_size07.smi > data/gdb7_protonated.smi


# TODO Conformational search
coord:
	mkdir -p xyz/
	python get_structures.py --csv data/gdb7_protonated.smi --scr xyz/

# TODO Optimize

# TODO Free energy

# TODO Parse to reference.csv



