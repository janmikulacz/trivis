#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import argparse
from termcolor import colored
from pathlib import Path
import math

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
                         default='./outputs/tuning4_weights',
                         help='Path to the output directory.')

maps = [
    ('scene_mp_2p_03', 1.0),
    ('scene_mp_4p_02', 1.0),
    ('scene_mp_6p_02', 1.0),
    ('scene_sp_cha_01', 1.0),
    ('scene_sp_pol_04', 1.0),
    ('scene_sp_pol_06', 1.0),
    ('scene_sp_rus_03', 1.0),
    ('scene_sp_rus_05', 1.0),
    ('scene_sp_sax_02', 1.0),
    ('scene_sp_sax_07', 1.0)
]

mwt_iter_limit = 200
mwt_time_limit = 10
mwt_simple_poly_max_size = 450

long_edge_thresholds = [
    math.nan,  # CDT
    -1.0,  # MinVT without optimization
    1.0,
    0.9,
    0.8,
    0.7,
    0.6,
    0.5,
    0.4,
    0.3,
    0.2,
    0.1
]


def main(args):
    args_common = ''
    args_common = f'{args_common} --bucket-size 1.0'
    args_common = f'{args_common} --mean-bucket-triangle-count-max -1.0'

    if args.generate_meshes:

        args_gen_queue = []

        args_gen_common = f'{args_common}'
        args_gen_common = f'{args_gen_common} --random-seed-method 4200'
        args_gen_common = f'{args_gen_common} --eps-dist-diff-collinear 1e-7'
        args_gen_common = f'{args_gen_common} --max-sample-dist 2'
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
            for threshold in long_edge_thresholds:
                args_gen_threshold = f'{args_gen_map}'
                if math.isnan(threshold):
                    args_gen_threshold = f'{args_gen_threshold} --methods CDT'
                    args_gen_threshold = f'{args_gen_threshold} --out-dir {args.output_dir}/standard'
                    args_gen_queue.append(args_gen_threshold)
                else:
                    args_gen_threshold = f'{args_gen_threshold} --methods MinVT'
                    if threshold >= 0.0:
                        args_gen_threshold = f'{args_gen_threshold} --weights-long-edge-penalty 1000.0'
                        args_gen_threshold = f'{args_gen_threshold} --weights-long-edge-threshold {threshold}'
                        args_gen_threshold = f'{args_gen_threshold} --out-dir {args.output_dir}/long-edge-threshold-{threshold}'
                    else:
                        args_gen_threshold = f'{args_gen_threshold} --out-dir {args.output_dir}/standard'
                    args_gen_queue.append(args_gen_threshold)

        print(colored('=== PRINTING ARGUMENTS QUEUE ===', 'green'))
        for args_gen in args_gen_queue:
            print(args_gen)
        print(colored('================================', 'green'))

        time.sleep(5)

        for i, args_gen in enumerate(args_gen_queue):
            cmd = f'{args.executable_gen} {args_gen}'
            print(colored(f'{i}/{len(args_gen_queue)}: Running "{cmd}"', 'cyan'))
            os.system(cmd)
            print(colored(f'{i}/{len(args_gen_queue)}: Finished "{cmd}"', 'magenta'))

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
        for args_gen in args_eval_queue:
            print(args_gen)
        print(colored('================================', 'green'))

        time.sleep(5)

        for i, args_eval in enumerate(args_eval_queue):
            cmd = f'{args.executable_eval} {args_eval}'
            print(colored(f'{i}/{len(args_eval_queue)}: Running "{cmd}"', 'cyan'))
            os.system(cmd)
            print(colored(f'{i}/{len(args_eval_queue)}: Finished "{cmd}"', 'magenta'))


if __name__ == '__main__':
    main(args_parser.parse_args())
