import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib import gridspec
from pymatgen.io.phonopy import get_ph_bs_symm_line_from_dict
from pymatgen.phonon.plotter import PhononBSPlotter
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='.', help='Path to investigate')
    parser.add_argument('--poscar', type=str, default='POSCAR', help='POSCAR file to open')
    parser.add_argument('--fs', type=int, default=20, help='fontsize')
    parser.add_argument('--width', type=int, default=12, help='image width')
    parser.add_argument('--height', type=int, default=6, help='image height')
    args = parser.parse_args()

    with open(os.path.join(args.path, 'qpoints.yaml'), 'r') as f:
        qpoints = yaml.load(f, Loader=yaml.FullLoader)

    with open(os.path.join(args.path, 'band.yaml'), 'r') as f:
        band = yaml.load(f, Loader=yaml.FullLoader)

    dos = np.loadtxt(os.path.join(args.path, 'total_dos.dat'), unpack=True)

    bs = get_ph_bs_symm_line_from_dict(band, labels_dict=qpoints['kpoints'])
    plotter = PhononBSPlotter(bs)
    data = plotter.bs_plot_data()

    new_ticks = dict()
    new_ticks['distance'] = list()
    new_ticks['sep_distance'] = list()
    new_ticks['label'] = list()
    for i, dist in enumerate(data['ticks']['distance']):
        _dist = np.round(dist, 5)
        if _dist not in new_ticks['distance']:
            new_ticks['distance'].append(_dist)
            old_label = data['ticks']['label'][i].replace('$', '')
            if old_label == 'GAMMA':
                new_label = r'$\Gamma$'
            if old_label == 'SIGMA':
                new_label = r'$\Sigma$'
            else:
                if r'\mid' in old_label:
                    new_ticks['sep_distance'].append(_dist)
                    split = old_label.split(r'\mid')
                    sep = list()
                    for qpoint in split:
                        if qpoint == 'GAMMA':
                            qpoint = r'\Gamma'
                        if qpoint == 'SIGMA':
                            qpoint = r'\Sigma'
                        elif '_' in qpoint:
                            new_split = qpoint.split('_')
                            sep.append(r'\mathrm{' + new_split[0] + r'}_{' + new_split[1] + r'}')
                        else:
                            sep.append(r'\mathrm{' + qpoint + r'}')
                    new_label = r'$' + r'\mid'.join(sep) + r'$'
                else:
                    new_label = r'$\mathrm{' + old_label + r'}$'
            new_ticks['label'].append(new_label)
    new_ticks['label'][-1] = new_ticks['label'][-1] + r'$\mid 0$'

    font = {'weight': 'normal',
            'size': 20}
    mpl.rc('font', **font)

    fig = plt.figure(figsize=(args.width, args.height))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    gs.update(wspace=0, hspace=0)
    ax0 = plt.subplot(gs[0])
    # fig, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(12, 6))
    print('Plotting phonon bands ...')
    for i in tqdm(range(len(data['distances']))):
        for j, freq in enumerate(data['frequency'][i]):
            ax0.plot(data['distances'][i], freq, color='steelblue', zorder=0, linewidth=1)
    ax0.set_xticks(new_ticks['distance'])
    ax0.set_xticklabels(new_ticks['label'], fontsize=args.fs)
    ax0.autoscale(enable=True, axis='x', tight=True)
    ylim = np.max(data['frequency']) * 1.05
    ax0.vlines(new_ticks['distance'], 0, ylim, color='k', linestyle='--', zorder=1, linewidth=1)
    ax0.vlines(new_ticks['sep_distance'], 0, ylim, color='k', zorder=2, linewidth=2)
    ax0.set_ylim(0, ylim)
    ax0.set_ylabel('Frequency, THz')

    ax1 = plt.subplot(gs[1], sharey=ax0)
    # plt.subplots_adjust(wspace=0, hspace=0)
    # plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)

    ax1.fill(dos[1], dos[0], dos[1], np.zeros(dos[1].shape), facecolor='lightsteelblue', edgecolor='steelblue')
    ticks = np.array([np.round(np.max(dos[1]) / 2), 2 * np.round(np.max(dos[1]) / 2)])
    ax1.set_xticks(ticks)
    ax1.set_xlim(0)
    ax1.set_xlabel('Phonon DOS')
    fig.tight_layout()
    fig.savefig(os.path.join(args.path, 'band.pdf'), bbox_inches='tight')


if __name__ == '__main__':
    main()
