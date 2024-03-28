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

args_parser = argparse.ArgumentParser()
args_parser.add_argument('-i', '--input_csv_file', type=str,
                         default='csv/tuning2_mwt.csv',
                         help='Input table in CSV format.')


def main(args):
    # Load all results.
    tab_main = dt.Frame()
    tab_main = dt.fread(args.input_csv_file)

    # Cosmetics.
    tab_main.names = {'map_name': 'map'}
    tab_main.names = {'mwt_simple_polygon_max_size': 'mwt_poly_n'}
    tab_main.names = {'mwt_iteration': 'mwt_it'}
    tab_main.names = {'mwt_result_weight': 'mwt_w'}
    tab_main.names = {'time_construction': 'mwt_t'}

    tab_main = tab_main[:, :, dt.sort(f.mwt_poly_n)]
    tab_main = tab_main[:, :, dt.sort(f.map)]

    best_values_iter_limit = {}
    best_values_time_limit = {}
    best_values = set()
    map_names = dt.unique(tab_main['map']).to_list()[0]
    for map_name in map_names:
        print(f'Processing map {map_name}.')
        tab_map = tab_main[f.map == map_name, :]
        min_mwt_w = tab_map[f.mwt_w == dt.min(f.mwt_w), 'mwt_w'].to_list()[0][0]
        max_mwt_w = tab_map[f.mwt_w == dt.max(f.mwt_w), 'mwt_w'].to_list()[0][0]
        min_mwt_t = tab_map[f.mwt_t == dt.min(f.mwt_t), 'mwt_t'].to_list()[0][0]
        max_mwt_t = tab_map[f.mwt_t == dt.max(f.mwt_t), 'mwt_t'].to_list()[0][0]
        mwt_poly_ns = dt.unique(tab_main['mwt_poly_n']).to_list()[0]
        n = len(mwt_poly_ns)
        if -1 in mwt_poly_ns:
            mwt_poly_ns.remove(-1)
            mwt_poly_ns.append(-1)
        list_map = []
        best_res = []
        for mwt_poly_n in mwt_poly_ns:
            tab_mwt_poly_n = tab_map[f.mwt_poly_n == mwt_poly_n, ('mwt_it', 'mwt_w', 'mwt_t')]
            tab_mwt_poly_n = tab_mwt_poly_n[:, :, dt.sort(f.mwt_it)]
            tab_mwt_poly_n_numpy = tab_mwt_poly_n.to_numpy()
            mwt_w_rel_diff = tab_mwt_poly_n_numpy[:,
                             1]  # (tab_mwt_poly_n_numpy[:, 1] - min_mwt_w) / (max_mwt_w - min_mwt_w)
            max_i = np.argmax(mwt_w_rel_diff < 1e-10)
            max_i = max_i + 1 if max_i != 0 else len(mwt_w_rel_diff)
            dict_mwt_poly_n = {
                'mwt_it': tab_mwt_poly_n_numpy[:max_i, 0],
                'mwt_w': tab_mwt_poly_n_numpy[:max_i, 1],
                'mwt_w_rel_diff': mwt_w_rel_diff[:max_i],
                'mwt_t': tab_mwt_poly_n_numpy[:max_i, 2],
                'mwt_t_rel': tab_mwt_poly_n_numpy[:max_i, 2] / 360
            }
            best_res.append((mwt_poly_n, dict_mwt_poly_n['mwt_t'][max_i - 1], dict_mwt_poly_n['mwt_it'][max_i - 1]))
            list_map.append(dict_mwt_poly_n)
        best_res = sorted(best_res, key=lambda val: val[1])
        for i in range(2):
            val = best_res[i][0]
            val_t = best_res[i][1]
            val_it = best_res[i][2]
            best_values.add(val)
            if val not in best_values_time_limit:
                best_values_time_limit[val] = [val_t]
            else:
                best_values_time_limit[val].append(val_t)
            if val not in best_values_iter_limit:
                best_values_iter_limit[val] = [val_it]
            else:
                best_values_iter_limit[val].append(val_it)
        colors = pl.cm.jet(np.linspace(0, 1, n))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        legend = []
        for i, mwt_poly_n in enumerate(mwt_poly_ns):
            ax.plot(list_map[i]['mwt_t'], list_map[i]['mwt_w_rel_diff'],
                    '.-',
                    color=colors[i],
                    linewidth=2,
                    )
            legend.append('$n_\mathcal{P}=%s$' % ('\infty' if mwt_poly_n == -1 else str(mwt_poly_n)))
        ax.set_xlabel('Time to compute the MWT [s]')
        ax.set_ylabel('Mesh weight [-]')
        ax.legend(legend, loc='upper right')
        ax.autoscale(enable=True, axis='both')
        plt.savefig(f'figures/mwt_tuning-{map_name}.pdf', bbox_inches='tight')
        ax.set_xlim((-3.5, 100.0))
        ax.set_ylim((187.8, 198.5))
        ax.get_legend().remove()
        plt.savefig(f'figures/mwt_tuning-{map_name}-detail.pdf', bbox_inches='tight')
    best_values = sorted(list(best_values))
    print(best_values)
    time_limit = 0
    for _, val_t_list in best_values_time_limit.items():
        time_limit = max(time_limit, max(val_t_list))
    print(f'Time limit: {time_limit}')
    iter_limit = 0
    for _, val_it_list in best_values_iter_limit.items():
        iter_limit = max(iter_limit, max(val_it_list))
    print(f'Iter limit: {iter_limit}')


if __name__ == '__main__':
    main(args_parser.parse_args())
