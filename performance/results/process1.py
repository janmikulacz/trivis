#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import datatable as dt

from pathlib import Path
from datatable import f

args_parser = argparse.ArgumentParser()
args_parser.add_argument('-i', '--input_csv_file', type=str, default='csv/*.csv', help='Input table in CSV format.')


def multiply_and_round(tab, column, factor, decimals):
    tab[column] = dt.float64(dt.int64(factor * (10 ** decimals) * f[column])) / (10 ** decimals)


def compute_summary(tab_main, group_by):
    summary = dt.Frame()
    tab_time_prep = tab_main[:, {
                                    # Time to initialize the solver (e.g., construct map representation):
                                    'time_init': dt.min(f.time_init),
                                    # Time to compute the preprocessing (e.g., compute triangulation):
                                    'time_preprocess': dt.min(f.time_preprocess),
                                }, dt.by(['map_name', 'points', 'algorithm'])]
    tab_time_prep = tab_time_prep[:, {
                                         # Number of initializations:
                                         'num_init': dt.count(),
                                         # Mean and standard deviation of time to initialize the solver (e.g., construct map representation):
                                         'time_init_mean': dt.mean(f.time_init),
                                         'time_init_std': dt.sd(f.time_init),
                                         # Mean and standard deviation of time to compute the preprocessing (e.g., compute triangulation):
                                         'time_prep_mean': dt.mean(f.time_preprocess),
                                         'time_prep_std': dt.sd(f.time_preprocess),
                                     }, dt.by(group_by)]
    summary = dt.rbind(summary, tab_time_prep)
    del tab_time_prep
    tab_result = tab_main[:, {
                                 # Number of queries:
                                 'num_query': dt.count(),
                                 # Number and percentage of crashes:
                                 'num_crashed': dt.sum(f.crashed),
                                 '%crashed': 100.0 * dt.sum(f.crashed) / dt.count(),
                                 # Number and percentage of timeouts / infinite loops:
                                 'num_timeout': dt.sum(f.timeout),
                                 '%timeout': 100.0 * dt.sum(f.timeout) / dt.count(),
                                 # Number and percentage of 'found' vis. regions:
                                 'num_found': dt.sum(f.found),
                                 '%found': 100.0 * dt.sum(f.found) / dt.count(),
                                 # Number and percentage of 'found' vis. regions with unavailable GT:
                                 'num_f_gt_unavailable': dt.sum(f.f_gt_unavailable),
                                 '%f_gt_unavailable': 100.0 * dt.sum(f.f_gt_unavailable) / dt.count(),
                                 # Number and percentage of 'found' vis. regions with GT reporting 'not-found':
                                 'num_f_gt_not_found': dt.sum(f.f_gt_not_found),
                                 '%f_gt_not_found': 100.0 * dt.sum(f.f_gt_not_found) / dt.count(),
                                 # Number and percentage of 'found' vis. regions that match the GT:
                                 'num_f_gt_ok': dt.sum(f.f_gt_ok),
                                 '%f_gt_ok': 100.0 * dt.sum(f.f_gt_ok) / dt.count(),
                                 # Number and percentage of 'found' vis. regions not matching the GT that were also computed from weakly simple nodes:
                                 'num_f_gt_err_weakly_simple': dt.sum(f.f_gt_err_weakly_simple),
                                 '%f_gt_err_weakly_simple': 100.0 * dt.sum(f.f_gt_err_weakly_simple) / dt.count(),
                                 # Number and percentage of 'found' vis. regions not matching the GT that were also computed from snapped nodes:
                                 'num_f_gt_err_snap': dt.sum(f.f_gt_err_snap),
                                 '%f_gt_err_snap': 100.0 * dt.sum(f.f_gt_err_snap) / dt.count(),
                                 # Number and percentage of 'found' vis. regions not matching the GT that do not fall into any of the previous categories:
                                 'num_f_gt_err_other': dt.sum(f.f_gt_err_other),
                                 '%f_gt_err_other': 100.0 * dt.sum(f.f_gt_err_other) / dt.count(),
                                 # Number and percentage of 'not-found' vis. regions:
                                 'num_not_found': dt.sum(f.not_found),
                                 '%not_found': 100.0 * dt.sum(f.not_found) / dt.count(),
                                 # Number and percentage of 'not-found' vis. regions with unavailable GT:
                                 'num_n_gt_unavailable': dt.sum(f.n_gt_unavailable),
                                 '%n_gt_unavailable': 100.0 * dt.sum(f.n_gt_unavailable) / dt.count(),
                                 # Number and percentage of 'not-found' vis. regions with GT reporting 'found':
                                 'num_n_gt_found': dt.sum(f.n_gt_found),
                                 '%n_gt_found': 100.0 * dt.sum(f.n_gt_found) / dt.count(),
                                 # Number and percentage of 'not-found' vis. regions with GT also reporting 'not-found':
                                 'num_n_gt_ok': dt.sum(f.n_gt_ok),
                                 '%n_gt_ok': 100.0 * dt.sum(f.n_gt_ok) / dt.count(),
                             }, dt.by(group_by)]
    # Number and percentage of queries with a crash or timeout / infinite loop:
    tab_result['num_crashed_timeout'] = f.num_crashed + f.num_timeout
    tab_result['%crashed_timeout'] = 100.0 * f.num_crashed_timeout / f.num_query
    # Number and percentage of queries with unavailable GT:
    tab_result['num_gt_unavailable'] = f.num_f_gt_unavailable + f.num_n_gt_unavailable
    tab_result['%gt_unavailable'] = 100.0 * f.num_gt_unavailable / f.num_query
    # Number and percentage of queries that don't match the GT 'found' or 'not-found':
    tab_result['num_gt_found_missmatch'] = f.num_f_gt_not_found + f.num_n_gt_found
    tab_result['%gt_found_missmatch'] = 100.0 * f.num_gt_found_missmatch / f.num_query
    # Number and percentage of queries that match the GT:
    tab_result['num_gt_ok'] = f.num_f_gt_ok + f.num_n_gt_ok
    tab_result['%gt_ok'] = 100.0 * f.num_gt_ok / f.num_query
    # Number and percentage of queries that don't match the GT but have feature-related output:
    tab_result['num_gt_err_feature'] = f.num_f_gt_err_weakly_simple + f.num_f_gt_err_snap
    tab_result['%gt_err_feature'] = 100.0 * f.num_gt_err_feature / f.num_query
    # Number and percentage of queries that don't match the GT and don't have feature-related output:
    tab_result['num_gt_err_other'] = f.num_f_gt_err_other
    tab_result['%gt_err_other'] = 100.0 * f.num_gt_err_other / f.num_query
    del tab_result[:, group_by]
    summary = dt.cbind(summary, tab_result)
    del tab_result
    tab_time_query_found = tab_main[f.found, :]
    tab_time_query_found = tab_time_query_found[:, {
                                                       # Total time to locate all queries (in seconds):
                                                       'time_locate': dt.sum(f.time_locate),
                                                       # Total time to compute all vis. regions (in seconds):
                                                       'time_query': dt.sum(f.time_query),
                                                       # Total time to compute all intersections in the vis. regions (in seconds):
                                                       'time_intersect': dt.sum(f.time_intersect),
                                                       # Total time to compute all vis. regions including the intersections (in seconds):
                                                       'time_query_intersect': dt.sum(f.time_query_intersect),
                                                       # Total time to answer all queries (in seconds):
                                                       'time_total': dt.sum(f.time_total),
                                                       # Mean and standard deviation of time to locate a single query (in seconds):
                                                       'time_found_mean': dt.mean(f.time_total),
                                                       'time_found_std': dt.sd(f.time_total),
                                                   }, dt.by(group_by)]
    # Percentage of time spend locating the queries:
    tab_time_query_found['%time_locate'] = 100.0 * f.time_locate / f.time_total
    # Percentage of time spend computing the vis. regions:
    tab_time_query_found['%time_query'] = 100.0 * f.time_query / f.time_total
    # Percentage of time spend computing the intersections in the vis. regions:
    tab_time_query_found['%time_intersect'] = 100.0 * f.time_intersect / f.time_total
    # Percentage of time spend computing the vis. regions including the intersections:
    tab_time_query_found['%time_query_intersect'] = 100.0 * f.time_query_intersect / f.time_total
    del tab_time_query_found[:, group_by]
    summary = dt.cbind(summary, tab_time_query_found)
    del tab_time_query_found
    return summary


def main(args):
    # Load all results.
    tab_main = dt.Frame()
    if '*' in args.input_csv_file:
        for file in [str(path.resolve()) for path in Path('').glob(args.input_csv_file)]:
            print(f'Reading {file}.')
            tab_temp = dt.fread(file)
            tab_main = dt.rbind(tab_main, tab_temp)
    else:
        print(f'Reading {args.input_csv_file}.')
        tab_main = dt.fread(args.input_csv_file)
    point_sets_order = [
        'In',
        'BB',
        'Ver',
        'NearV',
        'Mid',
        'NearM',
    ]
    algorithms_renames = {
        'cgal_triangular_expansion': 'CGAL-TEA-CE',
        'cgal_triangular_expansion_inexact_constructions': 'CGAL-TEA-CI',
        'cgal_rotational_sweep': 'CGAL-RSA-CE',
        'cgal_rotational_sweep_inexact_constructions': 'CGAL-RSA-CI',
        'visilibity': 'VisiLibity1',
        'trivis': 'TriVis',
    }
    for algorithm, algorithm_new in algorithms_renames.items():
        tab_main[f.algorithm == algorithm, 'algorithm'] = algorithm_new

    # Remove the invalid results: results are known to be invalid for map scene_sp_cha_02.
    tab_main = tab_main[f.map_name != 'scene_sp_cha_02', :]

    # Compute the GT-validity metrics.
    tab_main['gt_xor_area'] = f.gt_diff_area + f.diff_gt_area
    tab_main['gt_xor_area_ratio'] = f.gt_xor_area / f.map_area
    # Categorize the results.
    tab_main['crashed'] = f.crashed
    tab_main['timeout'] = (~f.crashed) & f.timeout
    tab_main['found'] = (~f.crashed) & (~f.timeout) & f.found
    tab_main['f_gt_unavailable'] = f.found & (~f.gt_available)
    tab_main['f_gt_not_found'] = f.found & (~f.f_gt_unavailable) & (~f.gt_found)
    tab_main['f_gt_ok'] = f.found & (~f.f_gt_unavailable) & (~f.f_gt_not_found) & (f.gt_xor_area_ratio < 1e-9)
    tab_main['f_gt_err_weakly_simple'] = f.found & (~f.f_gt_unavailable) & (~f.f_gt_not_found) & (~f.f_gt_ok) & f.weakly_simple
    tab_main['f_gt_err_snap'] = f.found & (~f.f_gt_unavailable) & (~f.f_gt_not_found) & (~f.f_gt_ok) & (~f.f_gt_err_weakly_simple) & f.snapped
    tab_main['f_gt_err_other'] = f.found & (~f.f_gt_unavailable) & (~f.f_gt_not_found) & (~f.f_gt_ok) & (~f.f_gt_err_weakly_simple) & (~f.f_gt_err_snap)
    tab_main['not_found'] = (~f.crashed) & (~f.timeout) & (~f.found)
    tab_main['n_gt_unavailable'] = f.not_found & (~f.gt_available)
    tab_main['n_gt_found'] = f.not_found & (~f.n_gt_unavailable) & f.gt_found
    tab_main['n_gt_ok'] = f.not_found & (~f.n_gt_unavailable) & (~f.n_gt_found)
    # Check that the error categorization is correct.
    tab_main['check1'] = f.crashed + f.timeout + f.found + f.not_found
    if tab_main[f.check1 != 1, :].nrows > 0:
        print(f'Result categorization is incorrect based on check1. Check the following rows:')
        print(tab_main[f.check1 != 1, :])
        exit(1)
    tab_main['check2'] = (f.crashed + f.timeout  # Crashed or timeout.
                          + (f.f_gt_unavailable + f.f_gt_not_found + f.f_gt_ok + f.f_gt_err_weakly_simple + f.f_gt_err_snap + f.f_gt_err_other)  # Found categories.
                          + (f.n_gt_unavailable + f.n_gt_found + f.n_gt_ok))  # Not found categories.
    if tab_main[f.check2 != 1, :].nrows > 0:
        print(f'Result categorization is incorrect based on check2. Check the following rows:')
        print(tab_main[f.check2 != 1, :])
        exit(1)

    # Compute summary table per points+algorithm and per algorithm and merge the results.
    tab_results_per_points = compute_summary(tab_main, ['points', 'algorithm'])
    tab_results_per_algorithm = compute_summary(tab_main, ['algorithm'])
    tab_results = dt.rbind(tab_results_per_points, tab_results_per_algorithm, force=True)
    tab_results[f.points == None, 'points'] = 'Overall'

    # Table cosmetics.
    algorithms_order = [algorithm for algorithm in algorithms_renames.values()]
    tab_temp = dt.Frame()
    for points in point_sets_order + ['Overall']:
        for algorithm in algorithms_order:
            tab_temp = dt.rbind(tab_temp, tab_results[(f.points == points) & (f.algorithm == algorithm), :])
    tab_results = tab_temp

    # Save the results.
    tab_results.to_csv('results.csv')


if __name__ == '__main__':
    main(args_parser.parse_args())
