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
args_parser.add_argument('-i', '--input_csv_file', type=str,
                         default='csv/tuning4_weights.csv',
                         help='Input table in CSV format.')


def main(args):
    # Load all results.
    tab_main = dt.Frame()
    if '*' in args.input_csv_file:
        for file in [str(path.resolve()) for path in Path('').glob(args.input_csv_file)]:
            print(f'Reading {file}.')
            tab_main = dt.rbind(tab_main, dt.fread(file))
    else:
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
    tab_main.names = {'weights_long_edge_threshold': 'w_threshold'}

    print(tab_main.names)

    map_names = dt.unique(tab_main['map']).to_list()[0]
    n_m = len(map_names)

    w_thresholds = dt.unique(tab_main['w_threshold']).to_list()[0]
    w_thresholds = np.array(w_thresholds)[np.array(w_thresholds) >= 0.0].tolist()

    expansions_rel_all = []
    times_rel_all = []
    for m, map_name in enumerate(map_names):
        print(f'Processing map {map_name}.')
        tab_map = tab_main[f.map == map_name, :]
        tab_map = tab_map[:, :, dt.sort(f.w_threshold)]
        expan_min_lt = tab_map[f.method == 'CDT', 'expan'].to_list()[0][0]
        time_min_vt = tab_map[(f.method == 'MinVT') & (f.w_threshold == -1), 'w_t'].to_list()[0][0]
        # print(time_min_vt)
        expansions_rel, times_rel = [], []
        for w, w_threshold in enumerate(w_thresholds):
            tab_threshold = tab_map[f.w_threshold == w_threshold, :]
            expansions_rel.append(100.0 * (tab_threshold['expan'].to_list()[0][0] - expan_min_lt) / expan_min_lt)
            times_rel.append(tab_threshold['w_t'].to_list()[0][0])
            # times_rel.append((tab_threshold['w_t'].to_list()[0][0] - time_min_vt) / time_min_vt)
        expansions_rel_all.append(expansions_rel)
        times_rel_all.append(times_rel)
    expansions_rel_means = np.mean(expansions_rel_all, axis=0)
    expansions_rel_stds = np.std(expansions_rel_all, axis=0)
    times_rel_means = np.mean(times_rel_all, axis=0)
    times_rel_stds = np.std(times_rel_all, axis=0)

    # print(np.array(w_thresholds)[expansions_rel_means - (-0.05983931937970009) < 0.005])
    # print(times_rel_means[w_thresholds.index(0.6)])
    # exit(0)

    w_thresholds = 100 * (1.0 - np.array(w_thresholds))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.fill_between(w_thresholds, times_rel_means - times_rel_stds, times_rel_means + times_rel_stds,
                    color='lightskyblue', alpha=0.2)
    ax.plot(w_thresholds, times_rel_means, '.-', color='lightskyblue', linewidth=2)
    ax.set_xlabel('How many edges are penalized $r_{\mathit{pen}}$ [%]')
    ax.set_ylabel('Time to compute the weights [s]')
    ax.legend(['Time as a function of $r_{\mathit{pen}} \;\;\; (\mathit{mean}\pm\mathit{std})$'], loc='upper left')
    ax.set_ylim(bottom=-150, top=1150)
    # ax.set_yscale('log')
    ax2 = ax.twinx()
    ax2.fill_between(w_thresholds, expansions_rel_means - expansions_rel_stds,
                     expansions_rel_means + expansions_rel_stds, color='black', alpha=0.2)
    ax2.plot(w_thresholds, expansions_rel_means, '.-', color='black', linewidth=2)
    ax2.set_ylabel('Expansion gap $\%\eta_{\mathcal{T}}$  from the CDT [%]')
    ax2.legend(['$\%\eta_{\mathcal{T}}$  as a function of $r_{\mathit{pen}} \;\;\; (\mathit{mean}\pm\mathit{std})$'], loc='lower left')
    ax2.set_ylim(bottom=-20, top=0)
    plt.savefig(f'figures/weights_tuning.pdf', bbox_inches='tight')


if __name__ == '__main__':
    main(args_parser.parse_args())
