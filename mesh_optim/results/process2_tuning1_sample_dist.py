#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import argparse
import numpy as np
import datatable as dt
import matplotlib.pyplot as plt

from pathlib import Path
from datatable import f

args_parser = argparse.ArgumentParser()
args_parser.add_argument('-i', '--input_csv_file', type=str,
                         default='csv/tuning1_sample_dist.csv',
                         help='Input table in CSV format.')


def main(args):
    # Load all results.
    tab_main = dt.Frame()
    tab_main = dt.fread(args.input_csv_file)

    # Cosmetics.
    tab_main.names = {'map_name': 'map'}
    tab_main.names = {'method_name': 'method'}
    tab_main.names = {'method_random_seed': 'rand_seed'}
    tab_main.names = {'mwt_simple_polygon_max_size': 'mwt_poly_n'}
    tab_main.names = {'mwt_iteration': 'mwt_it'}
    tab_main.names = {'mwt_result_weight': 'mwt_w'}
    tab_main.names = {'time_construction': 'mwt_t'}
    tab_main.names = {'time_weights': 'w_t'}
    tab_main.names = {'avg_expansions': 'expan'}
    tab_main.names = {'max_sample_dist': 'samp_dist'}
    tab_main.names = {'weights_shortest_edges_node_max': 'w_node_max'}
    tab_main = tab_main[:, :, dt.sort(f.samp_dist)]
    tab_main = tab_main[:, :, dt.sort(f.map)]

    tab_main['tc'] = f.w_t  # + f.mwt_t
    tab_main['expan_ref'] = 0.0
    map_names = dt.unique(tab_main['map']).to_list()[0]
    for map_name in map_names:
        expan_ref = tab_main[(f.map == map_name) & (f.method == 'CDT'), 'expan'].to_list()[0][0]
        tab_main[f.map == map_name, 'expan_ref'] = expan_ref
    tab_main['expan_rel'] = 100.0 * (f.expan - f.expan_ref) / f.expan_ref
    tab_mean = tab_main[:, {
                               'expan_rel_mean': dt.mean(f.expan_rel),
                               'expan_rel_std': dt.sd(f.expan_rel),
                               'tc_mean': dt.mean(f.tc),
                               'tc_std': dt.sd(f.tc),
                           }, dt.by('method', 'samp_dist')]
    tab_np = tab_mean[
        f.method == 'MinVT', ['samp_dist', 'expan_rel_mean', 'expan_rel_std', 'tc_mean', 'tc_std']].to_numpy()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ln1 = ax.plot(tab_np[:, 0], tab_np[:, 3], '.-', color='lightskyblue', linewidth=2)
    ax.fill_between(tab_np[:, 0], tab_np[:, 3] - tab_np[:, 4], tab_np[:, 3] + tab_np[:, 4], alpha=0.25, color='lightskyblue')
    # ax.grid()
    ax.set_ylabel('Time to compute the weights [s]')
    ax.set_xlabel('Maximal sample distance $d_{\mathit{samp}}$ [m]')
    ax.legend(['Time as a function of $d_{\mathit{samp}} \;\;\; (\mathit{mean}\pm\mathit{std})$'], loc='upper left')
    ax.set_ylim(bottom=-150, top=1150)
    # ax.legend(['Number of queries with this number of expansions'], loc='lower center')
    ax2 = ax.twinx()
    min_vt_exp_mean = np.array([tab_mean[f.method == 'MinLT', :].to_list()[2][0]] * len(tab_np[:, 0]))
    min_vt_exp_std = np.array([tab_mean[f.method == 'MinLT', :].to_list()[3][0]] * len(tab_np[:, 0]))
    ln2 = ax2.plot(tab_np[:, 0], tab_np[:, 1], '.-', color='black', linewidth=2)
    ax2.fill_between(tab_np[:, 0], tab_np[:, 1] - tab_np[:, 2], tab_np[:, 1] + tab_np[:, 2], alpha=0.25, color='black')
    #ax2.set_yscale('log')
    ax2.set_ylabel('Expansion gap $\%\eta_{\mathcal{T}}$  from the CDT [%]')
    ax2.legend(['$\%\eta_{\mathcal{T}}$  as a function of $d_{\mathit{samp}} \;\;\; (\mathit{mean}\pm\mathit{std})$'], loc='lower left')
    ax2.set_ylim(bottom=-20, top=0)
    plt.savefig('figures/sample_dist.pdf', bbox_inches='tight')



if __name__ == '__main__':
    main(args_parser.parse_args())
