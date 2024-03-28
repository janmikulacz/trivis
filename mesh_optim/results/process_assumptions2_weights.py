#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import argparse
import numpy as np
import datatable as dt
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

from pathlib import Path
from datatable import f

from scipy.interpolate import interp1d

args_parser = argparse.ArgumentParser()
args_parser.add_argument('-i', '--input_data_dir', type=str,
                         default='../experiments/outputs/assumptions',
                         help='Input directory.')


def main(args):
    expansions_files = [str(path.resolve()) for path in Path(args.input_data_dir).rglob('avg_expansions_tea_d-inf.txt')]
    for expansion_file in expansions_files:
        estimates_file = f'{os.path.dirname(expansion_file)}/exp_estimates.txt'
        expansions = np.loadtxt(expansion_file)
        estimates = np.loadtxt(estimates_file)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(estimates[3:, 0], estimates[3:, 1], '.-', color='lightskyblue', linewidth=2)
        ax.legend(['$\eta_{\mathcal{T}}^h$  as a function of $d_{\mathit{samp}}$ (bottom axis)'], loc='lower center')
        ax.set_xscale('log')
        ax.set_xlabel('Maximal sample distance $d_{\mathit{samp}}$ [m]')
        ax.set_ylabel('Mean expansions count [-]')
        ax.invert_xaxis()
        ax2 = ax.twiny()
        ax2.plot(list(range(1, len(expansions) + 1)), expansions, color='black', linewidth=1)
        ax2.legend(['$\eta_{\mathcal{T}}$  as a function of $m$ (top axis)'], loc='upper center')
        ax2.set_xscale('log')
        ax2.set_xlabel('Number of queries $m$ [-]', labelpad=8)
        bottom, top = ax2.get_ylim()
        ax.set_ylim([0.95 * bottom, top])
        fig_file_str = f'figures/assumption_weights_{"simple" if "simple" in estimates_file else "with_holes"}.pdf'
        plt.savefig(fig_file_str, bbox_inches='tight')
        print(100 * (estimates[-1, 1] - expansions[-1]) / expansions[-1])




if __name__ == '__main__':
    main(args_parser.parse_args())
