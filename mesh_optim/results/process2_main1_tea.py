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
                         default='csv/main1_tea.csv',
                         help='Input table in CSV format.')


def main(args):
    # Load all results.
    tab_main = dt.Frame()
    if '*' in args.input_csv_file:
        for file in [str(path.resolve()) for path in Path('').glob(args.input_csv_file)]:
            # print(f'Reading {file}.')
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
    # tab_main.names = {'avg_time_vis_reg_tea': 'q_t'}
    tab_main.names = {'weights_shortest_edges_node_max': 'w_node_max'}
    tab_main.names = {'max_sample_dist': 'samp_dist'}
    tab_main.names = {'weights_long_edge_threshold': 'w_thresh'}
    tab_main[f.map == 'scene_sp_endmaps', 'map'] = 'scene_sp_em'

    tab_main['q_t'] = f.avg_time_vis_reg_bucketing + f.avg_time_vis_reg_tea + f.avg_time_vis_reg_intersections

    tab_main = tab_main[:, :, dt.sort(f.map)]

    method_tabs = []
    method_tabs.append(('CDT', tab_main[(f.method == 'CDT'), :]))
    method_tabs.append(('MinLT', tab_main[(f.method == 'MinLT'), :]))
    method_tabs.append(('MinVT-2', tab_main[(f.method == 'MinVT') & (f.samp_dist == 2.0) & (f.w_thresh == -1.0), :]))
    method_tabs.append(('MinVT-4', tab_main[(f.method == 'MinVT') & (f.samp_dist == 4.0) & (f.w_thresh == -1.0), :]))
    method_tabs.append(('MinVT-8', tab_main[(f.method == 'MinVT') & (f.samp_dist == 8.0) & (f.w_thresh == -1.0), :]))
    method_tabs.append(('MinVTO-2', tab_main[(f.method == 'MinVT') & (f.samp_dist == 2.0) & (f.w_thresh == 0.5), :]))
    method_tabs.append(('MinVTO-4', tab_main[(f.method == 'MinVT') & (f.samp_dist == 4.0) & (f.w_thresh == 0.5), :]))
    method_tabs.append(('MinVTO-8', tab_main[(f.method == 'MinVT') & (f.samp_dist == 8.0) & (f.w_thresh == 0.5), :]))
    method_tabs.append(('MaxLT', tab_main[(f.method == 'MaxLT'), :]))
    method_tabs.append(('MaxVT-2', tab_main[(f.method == 'MaxVT') & (f.samp_dist == 2.0) & (f.w_thresh == -1.0), :]))

    map_names = method_tabs[0][1][:, {'expan': dt.mean(f.expan)}, dt.by('map')][:, :, dt.sort('expan')].to_list()[0]
    # map_names = dt.unique(tab_main['map']).to_list()[0]
    n_m = len(map_names)

    tex_tab_head = \
        '''\\begin{sidewaystable}
\t\\sidewaystablefn%
\t\\begin{center}
\t\t\\begin{minipage}{\\textheight}
\t\t\t\\setlength{\\tabcolsep}{3pt}
\t\t\t\\caption{Performance of the TEA with various triangular meshes.}
\t\t\t\\label{tab:res-tea}
\t\t\t\\begin{tiny}
\t\t\t\t\\begin{tabular*}{\\textheight}{@{\\extracolsep{\\fill}}l*{29}{r}@{\\extracolsep{\\fill}}}
\t\t\t\t\t\\toprule
\t\t\t\t\t\t& \\multicolumn{2}{r}{CDT} & \\multicolumn{3}{r}{MinLT} & \\multicolumn{3}{r}{MinVT-2} & \\multicolumn{3}{r}{MinVT-4} & \\multicolumn{3}{r}{MinVT-8} & \\multicolumn{3}{r}{MinVTO-2} & \\multicolumn{3}{r}{MinVTO-4} & \\multicolumn{3}{r}{MinVTO-8} & \\multicolumn{3}{r}{MaxLT} & \\multicolumn{3}{r}{MaxVT-2} \\\\
\t\t\t\t\t\\cmidrule(l){2-3} \\cmidrule(l){4-6} \\cmidrule(l){7-9} \\cmidrule(l){10-12} \\cmidrule(l){13-15} \\cmidrule(l){16-18} \\cmidrule(l){19-21} \\cmidrule(l){22-24} \\cmidrule(l){25-27} \\cmidrule(l){28-30}
\t\t\t\t\t\tmap & $\\eta_\\mathcal{T}$ & $\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $t_c$ & $\\%\\eta_\\mathcal{T}$ & $\\%\\bar{t}_q$ \\\\
\t\t\t\t\t\\midrule'''

    stats = {}
    tex_tab_body = ''
    for m, map_name in enumerate(map_names):
        if m != 0 and m % 5 == 0:
            tex_tab_body += f'\t\t\t\t\t\\midrule\n'
        # print(f'Processing map {map_name}.')
        ref_tab = tab_main[(f.map == map_name) & (f.method == 'CDT'), :]
        ref_expan = ref_tab['expan'].to_list()[0][0]
        ref_q_t = ref_tab['q_t'].to_list()[0][0]
        tex_line = '\t\t\t\t\t\t'
        tex_line += map_name.replace('scene_', '').replace('_', '')
        for method, method_tab in method_tabs:
            expan = method_tab[f.map == map_name, 'expan'].to_list()[0][0]
            q_t = method_tab[f.map == map_name, 'q_t'].to_list()[0][0]
            w_t = method_tab[f.map == map_name, 'w_t'].to_list()[0][0]
            mwt_t = method_tab[f.map == map_name, 'mwt_t'].to_list()[0][0]
            c_t = w_t + mwt_t
            if method != 'CDT':
                # Compute relative change.
                expan = 100 * (expan - ref_expan) / ref_expan
                q_t = 100 * (q_t - ref_q_t) / ref_q_t
            if method not in stats:
                stats[method] = []
            stats[method].append((c_t, q_t, expan))
            if method == 'CDT':
                tex_line += f' & {expan:.1f}'
                tex_line += f' & {1e6 * q_t:.2f}'
            else:
                tex_line += f' & {c_t:.0f}'
                tex_line += f' & {expan:.0f}'
                tex_line += f' & {q_t:.0f}'
        tex_line += f' \\\\\n'
        tex_tab_body += tex_line
    tex_tab_body += f'\t\t\t\t\t\\midrule\n'
    n_col = 0
    tex_line = '\t\t\t\t\t\t'
    for method, method_stats in stats.items():
        method_stats = np.array(method_stats)
        c_t, q_t, expan = np.mean(method_stats, axis=0)
        method_vals = ''
        if method == 'CDT':
            method_vals += f' & {expan:.1f}'
            method_vals += f' & {1e6 * q_t:.2f}'
        else:
            method_vals += f' & {c_t:.0f}'
            method_vals += f' & {expan:.0f}'
            method_vals += f' & {q_t:.0f}'
        tex_line += method_vals
        n_col += 3
        # print(method, method_vals)
    tex_line += f' \\\\\n'
    tex_tab_body += tex_line
    tex_tab_tail = \
        '''\t\t\t\t\t\\cmidrule(l){2-3} \\cmidrule(l){4-6} \\cmidrule(l){7-9} \\cmidrule(l){10-12} \\cmidrule(l){13-15} \\cmidrule(l){16-18} \\cmidrule(l){19-21} \\cmidrule(l){22-24} \\cmidrule(l){25-27} \\cmidrule(l){28-30}
\t\t\t\t\t\t& \\multicolumn{2}{r}{CDT} & \\multicolumn{3}{r}{MinLT} & \\multicolumn{3}{r}{MinVT-2} & \\multicolumn{3}{r}{MinVT-4} & \\multicolumn{3}{r}{MinVT-8} & \\multicolumn{3}{r}{MinVTO-2} & \\multicolumn{3}{r}{MinVTO-4} & \\multicolumn{3}{r}{MinVTO-8} & \\multicolumn{3}{r}{MaxLT} & \\multicolumn{3}{r}{MaxVT-2} \\\\
\t\t\t\t\t\\botrule
\t\t\t\t\\end{tabular*}
\t\t\t\\end{tiny}
\t\t\t\\footnotetext{
\t\t\t\t\\emph{Legend:} $\\eta_\\mathcal{T}$ is the mean number of edge expansions; $\\bar{t}_q$ is the mean query computational time; $t_c$ is the construction time of the mesh; $\\%\\mathit{val}$ means a relative percentage gap of $\\mathit{val}$ from the CDT value ($\\mathit{val}_\\mathrm{CDT}$).
\t\t\t\tEvery value of $\\eta_\\mathcal{T}$ and $\\bar{t}_q$ is based on one million queries; the values in the \\textbf{avg} row are averaged over all test maps.
\t\t\t}
\t\t\t\\footnotetext{
\t\t\t\t\\emph{Units:} $\\eta_\\mathcal{T}$ is a count (the unit is 1); $\\bar{t}_q$ is in microseconds, $t_c$ is in seconds; $\\%\\eta_\\mathcal{T}$, and $\\%\\bar{t}_q$ are in percents.
\t\t\t}
\t\t\\end{minipage}
\t\\end{center}
\\end{sidewaystable}'''

    with open('tables/tea.tex', 'w') as file:
        file.write(tex_tab_head + '\n' + tex_tab_body + '\n' + tex_tab_tail)


if __name__ == '__main__':
    main(args_parser.parse_args())
