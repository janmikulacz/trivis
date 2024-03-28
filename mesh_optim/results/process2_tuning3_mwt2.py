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
                         default='csv/tuning3_mwt2*.csv',
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
    tab_main.names = {'method_random_seed': 'rand_seed'}
    tab_main.names = {'mwt_simple_polygon_max_size': 'mwt_poly_n'}
    tab_main.names = {'mwt_iteration': 'mwt_it'}
    tab_main.names = {'mwt_result_weight': 'mwt_w'}
    tab_main.names = {'time_construction': 'mwt_t'}

    map_names = dt.unique(tab_main['map']).to_list()[0]
    n_m = len(map_names)
    mwt_poly_ns = dt.unique(tab_main['mwt_poly_n']).to_list()[0]
    if -1 in mwt_poly_ns:
        mwt_poly_ns.remove(-1)
        mwt_poly_ns.append(-1)
    n_p = len(mwt_poly_ns)
    rand_seeds = dt.unique(tab_main['rand_seed']).to_list()[0]
    n_r = len(rand_seeds)
    all_times = np.arange(0, 50.01, 0.01)
    n_t = len(all_times)
    mwt_w_rel_diffs = np.zeros((n_m, n_p, n_r, n_t))
    mwt_w_rel_diffs[:] = np.nan
    mwt_w_it = np.zeros((n_m, n_p, n_r, n_t))
    mwt_w_it[:] = np.nan
    for m, map_name in enumerate(map_names):
        print(f'Processing map {map_name}.')
        tab_map = tab_main[f.map == map_name, :]
        min_mwt_w = tab_map[f.mwt_w == dt.min(f.mwt_w), 'mwt_w'].to_list()[0][0]
        max_mwt_w = tab_map[f.mwt_w == dt.max(f.mwt_w), 'mwt_w'].to_list()[0][0]
        for p, mwt_poly_n in enumerate(mwt_poly_ns):
            tab_mwt_poly_n = tab_map[f.mwt_poly_n == mwt_poly_n, ('mwt_it', 'mwt_w', 'mwt_t', 'rand_seed')]
            for r, rand_seed in enumerate(rand_seeds):
                tab_rand_seed = tab_mwt_poly_n[f.rand_seed == rand_seed, :]
                tab_rand_seed = tab_rand_seed[:, :, dt.sort(f.mwt_it)]
                tab_rand_seed_numpy = tab_rand_seed.to_numpy()
                mwt_w_rel_diff = (tab_rand_seed_numpy[:, 1] - min_mwt_w) / (max_mwt_w - min_mwt_w)
                times = tab_rand_seed_numpy[:, 2]
                iterations = tab_rand_seed_numpy[:, 0]
                interp_mwt_rel_diff = interp1d(times, mwt_w_rel_diff)
                interp_mwt_it = interp1d(times, iterations)
                k_start = np.argwhere(all_times > times[0])
                k_start = 0 if len(k_start) == 0 else k_start[0, 0]
                k_end = np.argwhere(all_times > times[-1])
                k_end = len(all_times) if len(k_end) == 0 else k_end[0, 0]
                mwt_w_rel_diffs[m, p, r, k_start:k_end] = interp_mwt_rel_diff(all_times[k_start:k_end])
                mwt_w_it[m, p, r, k_start:k_end] = interp_mwt_it(all_times[k_start:k_end])
    mwt_w_rel_diffs_mean_map_random_seed = np.nanmean(mwt_w_rel_diffs, axis=(0, 2))
    mwt_w_rel_diffs_std_map_random_seed = np.nanstd(mwt_w_rel_diffs, axis=(0, 2))
    mwt_w_it_max_map_random_seed = np.nanmax(mwt_w_it, axis=(0, 2))

    colors = pl.cm.jet(np.linspace(0, 1, n_p))
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    legend = []
    threshold = 1e-2
    below_threshold_times = np.zeros(n_p)
    below_threshold_times[:] = np.nan
    for p, mwt_poly_n in enumerate(mwt_poly_ns):
        means = mwt_w_rel_diffs_mean_map_random_seed[p, :]
        stds = mwt_w_rel_diffs_std_map_random_seed[p, :]
        # ax.fill_between(all_times, means - stds, means + stds, color=colors[p], alpha=0.2)
        ax.plot(all_times, means, color=colors[p], linewidth=2.0)
        legend.append('$n_\mathcal{P}=%s$' % ('\infty' if mwt_poly_n == -1 else str(mwt_poly_n)))
        times = all_times[means < 1e-2]
        if len(times) > 0:
            below_threshold_times[p] = times[0]
    ax.set_xlabel('Time to compute the MWT [s]')
    ax.set_ylabel('Mesh weight [-]')
    ax.legend(legend, loc='best')
    ax.autoscale(enable=True, axis='both')
    plt.savefig(f'figures/mwt_tuning-2.pdf', bbox_inches='tight')
    ax.set_xlim((0.5, 6))
    ax.set_ylim((0, 0.3))
    ax.get_legend().remove()
    plt.savefig(f'figures/mwt_tuning-2-detail.pdf', bbox_inches='tight')
    # plt.show()

    time_optimal = np.nanmin(below_threshold_times)
    p_optimal = np.nanargmin(below_threshold_times)
    mwt_poly_n_optimal = mwt_poly_ns[p_optimal]
    print(f'Optimal mwt_poly_n value: {mwt_poly_n_optimal}')
    print(f'Optimal iterations limit: {int(mwt_w_it_max_map_random_seed[p_optimal, all_times.tolist().index(time_optimal)])}')
    # print(mwt_w_it_max_map_random_seed[p_optimal, :])

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    opt_means = mwt_w_rel_diffs_mean_map_random_seed[p_optimal, :]
    opt_stds = mwt_w_rel_diffs_std_map_random_seed[p_optimal, :]
    ax.fill_between(all_times, opt_means - opt_stds, opt_means + opt_stds, color=colors[p_optimal], alpha=0.2)
    ax.plot(all_times, opt_means, color=colors[p_optimal])
    ax.autoscale(enable=True, axis='both')
    # plt.show()


if __name__ == '__main__':
    main(args_parser.parse_args())
