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
                         default='csv/main2_dtea.csv',
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
    tab_main.names = {'method_vis_radius': 'd'}
    tab_main.names = {'vis_radius': 'd_tea'}
    tab_main.names = {'weights_shortest_edges_node_max': 'w_node_max'}
    tab_main.names = {'max_sample_dist': 'samp_dist'}
    tab_main.names = {'weights_long_edge_threshold': 'w_thresh'}
    tab_main[f.map == 'scene_sp_endmaps', 'map'] = 'scene_sp_em'

    tab_main['q_t'] = f.avg_time_vis_reg_bucketing + f.avg_time_vis_reg_tea + f.avg_time_vis_reg_intersections + f.avg_time_vis_reg_radial

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
    # method_tabs.append(('MaxLT', tab_main[(f.method == 'MaxLT'), :]))
    # method_tabs.append(('MaxVT-2', tab_main[(f.method == 'MaxVT') & (f.samp_dist == 2.0) & (f.w_thresh == -1.0), :]))
    method_tabs.append(
        ('d-MinVT-2', tab_main[(f.method == 'd-MinVT') & (f.samp_dist == 2.0) & (f.w_thresh == -1.0), :]))
    method_tabs.append(
        ('d-MinVT-4', tab_main[(f.method == 'd-MinVT') & (f.samp_dist == 4.0) & (f.w_thresh == -1.0), :]))
    method_tabs.append(
        ('d-MinVT-8', tab_main[(f.method == 'd-MinVT') & (f.samp_dist == 8.0) & (f.w_thresh == -1.0), :]))
    method_tabs.append(
        ('d-MinVTO-2', tab_main[(f.method == 'd-MinVT') & (f.samp_dist == 2.0) & (f.w_thresh == 0.5), :]))
    method_tabs.append(
        ('d-MinVTO-4', tab_main[(f.method == 'd-MinVT') & (f.samp_dist == 4.0) & (f.w_thresh == 0.5), :]))
    method_tabs.append(
        ('d-MinVTO-8', tab_main[(f.method == 'd-MinVT') & (f.samp_dist == 8.0) & (f.w_thresh == 0.5), :]))

    vis_radii = dt.unique(tab_main['d_tea']).to_list()[0]
    if -1 in vis_radii:
        vis_radii = vis_radii[1:]
        vis_radii.append(-1)
    vis_radii.reverse()
    n_d = len(vis_radii)

    map_names = dt.unique(tab_main['map']).to_list()[0]
    n_m = len(map_names)
    tex_tab_head = \
        '''\\begin{sidewaystable}
\t\\sidewaystablefn%
\t\\begin{center}
\t\t\\begin{minipage}{\\textheight}
\t\t\t\\setlength{\\tabcolsep}{3pt}
\t\t\t\\caption{Performance of the d-TEA with various triangular meshes for various values of $d$ (summary).}
\t\t\t\\label{tab:res-dtea}
\t\t\t\\begin{tiny}
\t\t\t\t\\begin{tabular*}{\\textheight}{@{\\extracolsep{\\fill}}l*{24}{r}@{\\extracolsep{\\fill}}}
\t\t\t\t\t\\toprule
\t\t\t\t\t\t& \\multicolumn{3}{r}{$d\\,{=}\\,\\infty$} & \\multicolumn{3}{r}{$d\\,{=}\\,128$} & \\multicolumn{3}{r}{$d\\,{=}\\,64$} & \\multicolumn{3}{r}{$d\\,{=}\\,32$} & \\multicolumn{3}{r}{$d\\,{=}\\,16$} & \\multicolumn{3}{r}{$d\\,{=}\\,8$} & \\multicolumn{3}{r}{$d\\,{=}\\,4$} & \\multicolumn{3}{r}{$d\\,{=}\\,2$} \\\\
\t\t\t\t\t\\cmidrule(l){2-4} \\cmidrule(l){5-7} \\cmidrule(l){8-10} \\cmidrule(l){11-13} \\cmidrule(l){14-16} \\cmidrule(l){17-19} \\cmidrule(l){20-22} \\cmidrule(l){23-25} 
\t\t\t\t\t\tmesh type & $\\bar{t}_c$ & $\\%\\bar{\\eta}_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $\\bar{t}_c$ & $\\%\\bar{\\eta}_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $\\bar{t}_c$ & $\\%\\bar{\\eta}_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $\\bar{t}_c$ & $\\%\\bar{\\eta}_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $\\bar{t}_c$ & $\\%\\bar{\\eta}_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $\\bar{t}_c$ & $\\%\\bar{\\eta}_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $\\bar{t}_c$ & $\\%\\bar{\\eta}_\\mathcal{T}$ & $\\%\\bar{t}_q$ & $\\bar{t}_c$ & $\\%\\bar{\\eta}_\\mathcal{T}$ & $\\%\\bar{t}_q$ \\\\
\t\t\t\t\t\\midrule'''
    all_stats = {}
    for d_i, d in enumerate(vis_radii):
        tab_d = tab_main[f.d_tea == d, :]
        stats = {}
        tex_tab_body = ''
        tex_tab_body += f'\\midrule\n'
        for m, map_name in enumerate(map_names):
            if m != 0 and m % 5 == 0:
                tex_tab_body += f'\\midrule\n'
            # print(f'Processing map {map_name}.')
            ref_tab = tab_d[(f.map == map_name) & (f.method == 'CDT'), :]
            ref_expan = ref_tab['expan'].to_list()[0][0]
            ref_q_t = ref_tab['q_t'].to_list()[0][0]
            tex_line = ''
            tex_line += map_name.replace('scene_', '').replace('_', '')
            for method, method_tab in method_tabs:
                if method_tab[f.d_tea == d, :].nrows == 0:
                    continue
                expan = method_tab[(f.map == map_name) & (f.d_tea == d), 'expan'].to_list()[0][0]
                q_t = method_tab[(f.map == map_name) & (f.d_tea == d), 'q_t'].to_list()[0][0]
                w_t = method_tab[(f.map == map_name) & (f.d_tea == d), 'w_t'].to_list()[0][0]
                mwt_t = method_tab[(f.map == map_name) & (f.d_tea == d), 'mwt_t'].to_list()[0][0]
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
        tex_tab_body += f'\\midrule\n'
        n_col = 0
        tex_line = 'avg'
        all_stats[d] = {}
        for method, method_stats in stats.items():
            method_stats = np.array(method_stats)
            c_t, q_t, expan = np.mean(method_stats, axis=0)
            all_stats[d][method] = (c_t, q_t, expan)
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
        tex_tab_body += f'\\bottomrule\n'
        # print(f'\n\n\nd = {d}\n\n')
        # print(tex_tab_body)

    tex_tab_body = ''
    for method, _ in method_tabs:
        tex_line = '\t\t\t\t\t\t'
        tex_line += f'{method}'
        for d in vis_radii:
            if method in all_stats[d]:
                c_t, q_t, expan = all_stats[d][method]
                method_vals = ''
                if method == 'CDT':
                    method_vals += f' & {c_t:.0f}'
                    method_vals += f' & {expan:.1f}'
                    method_vals += f' & {1e6 * q_t:.2f}'
                else:
                    method_vals += f' & {c_t:.0f}'
                    method_vals += f' & {expan:.1f}'
                    method_vals += f' & {q_t:.1f}'
                tex_line += method_vals
            else:
                tex_line += ' & - & - & - '
        tex_line += f' \\\\\n'
        tex_tab_body += tex_line
    tex_tab_tail = \
        '''\t\t\t\t\t\\cmidrule(l){2-4} \\cmidrule(l){5-7} \\cmidrule(l){8-10} \\cmidrule(l){11-13} \\cmidrule(l){14-16} \\cmidrule(l){17-19} \\cmidrule(l){20-22} \\cmidrule(l){23-25} 
\t\t\t\t\t\t& \\multicolumn{3}{r}{$d\\,{=}\\,\\infty$} & \\multicolumn{3}{r}{$d\\,{=}\\,128$} & \\multicolumn{3}{r}{$d\\,{=}\\,64$} & \\multicolumn{3}{r}{$d\\,{=}\\,32$} & \\multicolumn{3}{r}{$d\\,{=}\\,16$} & \\multicolumn{3}{r}{$d\\,{=}\\,8$} & \\multicolumn{3}{r}{$d\\,{=}\\,4$} & \\multicolumn{3}{r}{$d\\,{=}\\,2$} \\\\
\t\t\t\t\t\\botrule
\t\t\t\t\\end{tabular*}
\t\t\t\\end{tiny}
\t\t\t\\footnotetext{
\t\t\t\t\\emph{Legend:} $d$ is the visibility range. Otherwise, the legend is the same as in Tab.~\\ref{tab:res-tea}, but all the values are averaged over all test maps.
\t\t\t}
\t\t\t\\footnotetext{
\t\t\t\t\\emph{Units:} $d$ is in meters; $\\bar{t}_c$ is in seconds; $\\%\\bar{\\eta}_\\mathcal{T}$, and $\\%\\bar{t}_q$ are in percents except for CDT (see footnote\\footnotemark[1]).
\t\t\t}
\t\t\t\\footnotetext{
\t\t\t\t\\emph{Note:} The values marked with \\textsuperscript{$\\star$} also appear at the bottom of Tab.~\\ref{tab:res-tea} (\\textbf{avg} row). 
\t\t\t}
\t\t\t\\footnotetext[1]{
\t\t\t\tThe values in $\\%\\bar{\\eta}_\\mathcal{T}$ and $\\%\\bar{t}_q$ columns are absolute for the CDT: $\\bar{\\eta}_\\mathcal{T}$ is a mean count, $\\bar{t}_q$ is in microseconds.
\t\t\t\tThe same values are relative percentage gaps w.r.t. CDT for the other mesh types.
\t\t\t}
\t\t\\end{minipage}
\t\\end{center}
\\end{sidewaystable}'''

    with open('tables/dtea.tex', 'w') as file:
        file.write(tex_tab_head + '\n' + tex_tab_body + tex_tab_tail)


if __name__ == '__main__':
    main(args_parser.parse_args())
