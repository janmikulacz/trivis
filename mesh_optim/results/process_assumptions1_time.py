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
    expansions_files = [str(path.resolve()) for path in Path(args.input_data_dir).rglob('all_expansions_tea_d-inf.txt')]
    for expansion_file in expansions_files:
        csv_file = f'{expansion_file}.csv'
        if False and Path(csv_file).is_file():
            tab = dt.fread(csv_file)
        else:
            expansions = np.loadtxt(expansion_file)
            tab = dt.Frame(expansions)
            tab.names = ['exp', 'tb', 'tt', 'ti', 'tr']
            tab['exp'] = dt.int32
            tab.to_csv(csv_file)
        tab = tab[:, :, dt.sort('exp')]
        tab['t1'] = f.tb
        tab['t2'] = f.tb + f.tt
        tab['t3'] = f.tb + f.tt + f.ti
        tab['t4'] = f.tb + f.tt + f.ti + f.tr
        tab['t'] = f.tb + f.tt + f.ti + f.tr
        tab_mean = tab[:, {
                              'cnt': dt.count(f.exp),  # 1
                              'tt_mean': dt.mean(1e6 * f.tt),  # 2
                              'tt_std': dt.sd(1e6 * f.tt),  # 3
                              't_mean': dt.mean(1e6 * f.t),  # 4
                              't_std': dt.sd(1e6 * f.t),  # 5
                              't1_mean': dt.mean(1e6 * f.t1),  # 6
                              't1_std': dt.sd(1e6 * f.t1),  # 7
                              't2_mean': dt.mean(1e6 * f.t2),  # 8
                              't2_std': dt.sd(1e6 * f.t2),  # 9
                              't3_mean': dt.mean(1e6 * f.t3),  # 10
                              't3_std': dt.sd(1e6 * f.t3),  # 11
                              't4_mean': dt.mean(1e6 * f.t4),  # 12
                              't4_std': dt.sd(1e6 * f.t4),  # 13
                          }, dt.by('exp')]
        tab_np = tab_mean.to_numpy()
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.bar(tab_np[:, 0], tab_np[:, 1], color='gold', width=1.0)
        ax.set_xlabel('Number of expansions $\eta_q$ [-]')
        ax.set_ylabel('Query count (total: 1 million) [-]')
        ax.legend(['Query count for each value of $\eta_q$'], loc='lower right')
        ax.set_ylim(top=ax.get_ylim()[1] + 0.4 * (ax.get_ylim()[1] - ax.get_ylim()[0]))
        ax2 = ax.twinx()
        ax2.plot(tab_np[:, 0], tab_np[:, 10], color='black', linewidth=2.0)
        ax2.fill_between(tab_np[:, 0], tab_np[:, 10] - tab_np[:, 11], tab_np[:, 10] + tab_np[:, 11], alpha=0.25,
                         color='black')
        ax2.plot(tab_np[:, 0], tab_np[:, 8], color='#4a718a', linewidth=2.0)
        ax2.fill_between(tab_np[:, 0], tab_np[:, 8] - tab_np[:, 9], tab_np[:, 8] + tab_np[:, 9], alpha=0.25,
                         color='#4a718a')
        ax2.plot(tab_np[:, 0], tab_np[:, 6], color='lightskyblue', linewidth=2.0)
        ax2.fill_between(tab_np[:, 0], tab_np[:, 6] - tab_np[:, 7], tab_np[:, 6] + tab_np[:, 7], alpha=0.25,
                         color='lightskyblue')
        ax2.legend([
            '$t_\Delta + t_{\mathrm{TEA}} + t_\\times = t_q \;\;\; (\mathit{mean}\pm\mathit{std})$',
            '$t_\Delta + t_{\mathrm{TEA}} \;\;\; (\mathit{mean}\pm\mathit{std})$',
            '$t_\Delta \;\;\; (\mathit{mean}\pm\mathit{std})$',
        ], loc='best')
        ax2.set_ylabel('Computational time $t_q$ [μs]')
        ax2.set_ylim(bottom=ax2.get_ylim()[0] - 0.05 * (ax2.get_ylim()[1] - ax2.get_ylim()[0]))
        fig_file_str = f'figures/assumption_time_{"simple" if "simple" in expansion_file else "with_holes"}.pdf'
        plt.savefig(fig_file_str, bbox_inches='tight')


if __name__ == '__main__':
    main(args_parser.parse_args())
