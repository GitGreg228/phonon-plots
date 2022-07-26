import argparse
import os

import numpy as np
import yaml
from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.symmetry.kpath import KPathSeek


def read_structure(path):
    with open(os.path.join(path, 'phonopy_disp.yaml'), 'r') as f:
        struc = yaml.load(f, Loader=yaml.FullLoader)['unit_cell']
    lattice = struc['lattice']
    sites = list()
    for point in struc['points']:
        sites.append(PeriodicSite(species=point['symbol'], coords=point['coordinates'], lattice=lattice))
    return Structure.from_sites(sites)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='.', help='Path to investigate')
    parser.add_argument('--poscar', type=str, default='POSCAR', help='POSCAR file to open')
    args = parser.parse_args()

    if os.path.isfile(os.path.join(args.path, 'phonopy_disp.yaml')):
        print(f"Reading structure from {os.path.join(args.path, 'phonopy_disp.yaml')}")
        s = read_structure(args.path)
    else:
        print(f"Reading structure from {os.path.join(args.path, args.poscar)}")
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
