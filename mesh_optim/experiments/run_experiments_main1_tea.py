#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import argparse
from termcolor import colored
from pathlib import Path

args_parser = argparse.ArgumentParser()
args_parser.add_argument('-G', '--generate_meshes',
                         default=False, action='store_true',
                         help="The program will generate the meshes.")
args_parser.add_argument('-E', '--evaluate_meshes',
                         default=False, action='store_true',
                         help="The program will evaluate the existing meshes.")
args_parser.add_argument('-g', '--executable_gen', type=str,
                         default='../../build-Release/mesh_optim/generate_meshes',
                         help='Path to the generate_meshes executable.')
args_parser.add_argument('-e', '--executable_eval', type=str,
                         default='../../build-Release/mesh_optim/evaluate_mesh',
                         help='Path to the evaluate_mesh executable.')
args_parser.add_argument('-o', '--output_dir', type=str,
                         default='./outputs/main1_tea',
                         help='Path to the output directory.')

maps = [
    ('scene_mp_2p_01', 1.0),
    ('scene_mp_2p_02', 1.0),
    ('scene_mp_2p_04', 1.0),
    ('scene_mp_4p_01', 1.0),
    ('scene_mp_4p_03', 1.0),
    ('scene_mp_6p_01', 1.0),
    ('scene_mp_6p_03', 1.0),
    ('scene_sp_cha_02', 1.0),
    ('scene_sp_cha_03', 1.0),
    ('scene_sp_cha_04', 1.0),
    ('scene_sp_endmaps', 1.0),
    ('scene_sp_pol_01', 1.0),
    ('scene_sp_pol_02', 1.0),
    ('scene_sp_pol_03', 1.0),
    ('scene_sp_pol_05', 1.0),
    ('scene_sp_rus_01', 1.0),
    ('scene_sp_rus_02', 1.0),
    ('scene_sp_rus_04', 1.0),
    ('scene_sp_rus_06', 1.0),
    ('scene_sp_rus_07', 1.0),
    ('scene_sp_sax_01', 1.0),
    ('scene_sp_sax_03', 1.0),
    ('scene_sp_sax_04', 1.0),
    ('scene_sp_sax_05', 1.0),
    ('scene_sp_sax_06', 1.0),
]

bucket_size = 1.0
mean_bucket_triangle_count_max = -1.0  # invalid

random_seed_method = 420
eps_dist_diff_collinear = 1e-7
max_sample_dists = [
    2,
    4,
    8
]
w_long_edge_penalty = 1000.0
w_long_edge_thresholds = [
    -1.0,  # not optimized
    0.5
]
mwt_iter_limit = 200
mwt_time_limit = 10
mwt_simple_poly_max_size = 450
method_groups = [
    'CDT MinLT MaxLT',
    'MinVT MaxVT',
    # 'd-MinVT d-MaxVT'
]
vis_radii = '2 4 8 16 32'


def main(args):
    args_common = ''
    args_common = f'{args_common} --bucket-size {bucket_size}'
    args_common = f'{args_common} --mean-bucket-triangle-count-max {mean_bucket_triangle_count_max}'

    if args.generate_meshes:

        args_gen_queue = []

        args_gen_common = f'{args_common}'
        args_gen_common = f'{args_gen_common} --random-seed-method {random_seed_method}'
        args_gen_common = f'{args_gen_common} --eps-dist-diff-collinear {eps_dist_diff_collinear}'
        args_gen_common = f'{args_gen_common} --save-weights'
        args_gen_common = f'{args_gen_common} --load-weights'
        # args_gen_common = f'{args_gen_common} --precompute-weights'
        args_gen_common = f'{args_gen_common} --save-only-last-iteration'
        args_gen_common = f'{args_gen_common} --mwt-iter-limit {mwt_iter_limit}'
        args_gen_common = f'{args_gen_common} --mwt-time-limit {mwt_time_limit}'
        args_gen_common = f'{args_gen_common} --mwt-simple-polygon-max-size {mwt_simple_poly_max_size}'
        for map_name, map_scale in maps:
            args_gen_map = f'{args_gen_common}'
            args_gen_map = f'{args_gen_map} --map-name {map_name}'
            args_gen_map = f'{args_gen_map} --map-scale {map_scale}'
            for methods in method_groups:
                args_gen_methods = f'{args_gen_map}'
                args_gen_methods = f'{args_gen_methods} --methods {methods}'
                if 'd-' in methods:
                    args_gen_methods = f'{args_gen_methods} --vis-radii-mesh {vis_radii}'
                if 'VT' not in methods:
                    out_dir = f'{args.output_dir}/others'
                    args_gen_methods = f'{args_gen_methods} --out-dir {out_dir}'
                    args_gen_queue.append(args_gen_methods)
                else:
                    for samp_dist in max_sample_dists:
                        args_gen_max_sample_dist = f'{args_gen_methods}'
                        args_gen_max_sample_dist = f'{args_gen_max_sample_dist} --max-sample-dist {samp_dist}'
                        for threshold in w_long_edge_thresholds:
                            args_gen_w = f'{args_gen_max_sample_dist}'
                            if threshold < 0.0:
                                out_dir = f'{args.output_dir}/sample-dist-{samp_dist}/standard'
                                args_gen_w = f'{args_gen_w} --out-dir {out_dir}'
                            else:
                                args_gen_w = f'{args_gen_w} --weights-long-edge-penalty {w_long_edge_penalty}'
                                args_gen_w = f'{args_gen_w} --weights-long-edge-threshold {threshold}'
                                out_dir = f'{args.output_dir}/sample-dist-{samp_dist}/long-edge-thresh-{threshold}'
                                args_gen_w = f'{args_gen_w} --out-dir {out_dir}'
                            args_gen_queue.append(args_gen_w)

        print(colored('=== PRINTING ARGUMENTS QUEUE ===', 'green'))
        for i, args_gen in enumerate(args_gen_queue):
            print(f'\nArguments {i + 1}:\n{args_gen}')
        print(colored('================================', 'green'))

        print(f'Queue size: {len(args_gen_queue)}.')

        time.sleep(5)

        for i, args_gen in enumerate(args_gen_queue):
            cmd = f'{args.executable_gen} {args_gen}'
            print(colored(f'{i + 1}/{len(args_gen_queue)}: Running "{cmd}"', 'cyan'))
            os.system(cmd)
            print(colored(f'{i + 1}/{len(args_gen_queue)}: Finished "{cmd}"', 'magenta'))

    if args.evaluate_meshes:

        args_eval_queue = []

        args_eval_common = f'{args_common}'
        args_eval_common = f'{args_eval_common} --n-queries 1000000'
        args_eval_common = f'{args_eval_common} --random-seed-queries 42'
        args_eval_common = f'{args_eval_common} --evaluate-tea'
        args_eval_common = f'{args_eval_common} --vis-radii-test -1'

        for mesh_path in Path(args.output_dir).rglob('mesh.txt'):
            mesh_dir = os.path.dirname(mesh_path)
            if not os.path.isfile(f'{mesh_dir}/mesh_info.json'):
                print(colored(f'No mesh_info.json in {mesh_dir}!', 'red'))
                continue
            args_eval_mesh = f'{args_eval_common}'
            args_eval_mesh = f'{args_eval_mesh} --mesh-file {mesh_path}'
            args_eval_mesh = f'{args_eval_mesh} --mesh-info-file {mesh_dir}/mesh_info.json'
            args_eval_mesh = f'{args_eval_mesh} --out-dir {mesh_dir}'
            args_eval_queue.append(args_eval_mesh)

        print(colored('=== PRINTING ARGUMENTS QUEUE ===', 'green'))
        for args_eval in args_eval_queue:
            print(args_eval)
        print(colored('================================', 'green'))

        print(f'Queue size: {len(args_eval_queue)}.')

        time.sleep(5)

        for i, args_eval in enumerate(args_eval_queue):
            cmd = f'{args.executable_eval} {args_eval}'
            print(colored(f'{i + 1}/{len(args_eval_queue)}: Running "{cmd}"', 'cyan'))
            os.system(cmd)
            print(colored(f'{i + 1}/{len(args_eval_queue)}: Finished "{cmd}"', 'magenta'))


if __name__ == '__main__':
    main(args_parser.parse_args())
