from pymatgen.core.structure import Structure
from pymatgen.symmetry.kpath import KPathSeek, KPathBase, KPathSetyawanCurtarolo
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.plotter import PhononBSPlotter
from pymatgen.io.phonopy import get_ph_bs_symm_line_from_dict
from tqdm import tqdm
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

    with open(os.path.join(args.path, 'qpoints.yaml'), 'r') as f:
        qpoints = yaml.load(f, Loader=yaml.FullLoader)

    with open(os.path.join(args.path, 'band.yaml'), 'r') as f:
        band = yaml.load(f, Loader=yaml.FullLoader)
    bs = get_ph_bs_symm_line_from_dict(band, labels_dict=qpoints['kpoints'])
    plotter = PhononBSPlotter(bs)
    data = plotter.bs_plot_data()

    new_ticks = dict()
    new_ticks['distance'] = list()
    new_ticks['label'] = list()
    for i, dist in enumerate(data['ticks']['distance']):
        _dist = np.round(dist, 5)
        if _dist not in new_ticks['distance']:
            new_ticks['distance'].append(_dist)
            old_label = data['ticks']['label'][i].replace('$', '')
            if old_label == 'GAMMA':
                new_label = r'$\Gamma$'
            else:
                if r'\mid' in old_label:
                    split = old_label.split(r'\mid')
                    sep = list()
                    for qpoint in split:
                        if '_' in qpoint:
                            new_split = qpoint.split('_')
                            sep.append(r'\mathrm{' + new_split[0] + r'}_{' + new_split[1] + r'}')
                        else:
                            sep.append(r'\mathrm{' + qpoint + r'}')
                    new_label = r'$' + r'\mid'.join(sep) + r'$'
                else:
                    new_label = r'$\mathrm{' + old_label + r'}$'
            new_ticks['label'].append(new_label)

    font = {'weight': 'normal',
            'size': 15}
    mpl.rc('font', **font)

    print('Plotting phonon bands ...')
    plt.figure(figsize=(8, 6))
    for i in tqdm(range(len(data['distances']))):
        for j, freq in enumerate(data['frequency'][i]):
            plt.plot(data['distances'][i], freq, color='steelblue', zorder=0)
    plt.xticks(new_ticks['distance'], labels=new_ticks['label'])
    plt.autoscale(enable=True, axis='x', tight=True)
    ylim = np.max(data['frequency']) * 1.05
    plt.vlines(new_ticks['distance'], 0, ylim, color='k', zorder=1)
    plt.ylim(0, ylim)
    plt.ylabel('Frequency, THz')
    plt.savefig(os.path.join(args.path, 'band.pdf'))


if __name__ == '__main__':
    main()