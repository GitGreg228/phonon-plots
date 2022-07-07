from pymatgen.core.structure import Structure
from pymatgen.symmetry.kpath import KPathSeek, KPathBase, KPathSetyawanCurtarolo
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.plotter import PhononBSPlotter
from pymatgen.io.phonopy import get_ph_bs_symm_line_from_dict
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import json
import argparse
import os
import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='.', help='Path to investigate')
    parser.add_argument('--poscar', type=str, default='POSCAR', help='POSCAR file to open')
    args = parser.parse_args()

    s = Structure.from_file(os.path.join(args.path, args.poscar))
    k = KPathSeek(s, symprec=0.2)
    qpoints = k.kpath['kpoints']
    for qpoint in qpoints.keys():
        qpoints[qpoint] = np.round(qpoints[qpoint], 3).tolist()
    path = k.kpath['path']
    labels = list()
    kpt_list = list()
    for subpath in path:
        _labels = list()
        _kpt_list = list()
        for qpoint in subpath:
            _kpt_list.append(np.round(qpoints[qpoint], 3).tolist())
            _labels.append(qpoint)
        labels.append(_labels)
        kpt_list.append(_kpt_list)
    path = list()
    path_labels = list()
    for i, _kpt_list in enumerate(kpt_list):
        subpath = list()
        subpath_labels = list()
        for j, qpoint in enumerate(_kpt_list):
            subpath.append(' '.join([str(q) for q in qpoint]))
            subpath_labels.append(labels[i][j])
        path.append(' '.join(subpath))
        path_labels.append(' '.join(subpath_labels))

    with open(os.path.join(args.path, 'phonopy_disp.yaml'), 'r') as f:
        phonopy = yaml.load(f, Loader=yaml.FullLoader)
    dim = phonopy['phonopy']['configuration']['dim']
    tol = phonopy['phonopy']['configuration']['symmetry_tolerance']
    atom_name = set([specie.symbol for specie in s.species])
    print('BAND = ' + ', '.join(path))
    print('BAND_LABELS = ' + ' '.join(path_labels))
    with open(os.path.join(args.path, 'band.conf'), 'w') as f:
        f.write(f'DIM = {dim}\n')
        f.write(f'ATOM_NAME = {" ".join([a for a in atom_name])}\n')
        f.write(f'BAND = {", ".join(path)}\n')
        f.write(f'BAND_LABELS = {" ".join(path_labels)}\n')
        f.write(f'SYMMETRY_TOLERANCE = {tol}\n')
    with open(os.path.join(args.path, 'qpoints.yaml'), 'w') as f:
        yaml.dump(k.kpath, f)
    print(f"Please see {os.path.join(args.path, 'band.conf')} and run phonopy.")


if __name__ == '__main__':
    main()