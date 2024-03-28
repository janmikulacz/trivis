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
                         default='./outputs/tuning1_sample_dist',
                         help='Path to the output directory.')

maps = [
    ('scene_mp_2p_03_simple-1', 1.0),
    ('scene_mp_4p_02_simple-1', 1.0),
    ('scene_mp_6p_02_simple-1', 1.0),
    ('scene_sp_cha_01_simple-1', 1.0),
    ('scene_sp_pol_04_simple-1', 1.0),
    ('scene_sp_pol_06_simple-1', 1.0),
    ('scene_sp_rus_03_simple-1', 1.0),
    ('scene_sp_rus_05_simple-1', 1.0),
    ('scene_sp_sax_02_simple-1', 1.0),
    ('scene_sp_sax_07_simple-1', 1.0)
]

sample_dists = [
    128,
    64,
    32,
    16,
    8,
    4,
    2,
    1,
    0.5,
    0.25,
    0.125
]


def main(args):
    args_common = ''
    args_common = f'{args_common} --bucket-size 1.0'
    args_common = f'{args_common} --mean-bucket-triangle-count-max -1.0'

    if args.generate_meshes:

        args_gen_queue = []

        args_gen_common = f'{args_common}'
        args_gen_common = f'{args_gen_common} --random-seed-method 420'
        args_gen_common = f'{args_gen_common} --mwt-iter-limit 1'
        args_gen_common = f'{args_gen_common} --eps-dist-diff-collinear 1e-7'
        args_gen_common = f'{args_gen_common} --save-weights'
        args_gen_common = f'{args_gen_common} --load-weights'
        # args_gen_common = f'{args_gen_common} --precompute-weights'
        for map_name, map_scale in maps:
            args_gen_map = f'{args_gen_common}'
            args_gen_map = f'{args_gen_map} --map-name {map_name}'
            args_gen_map = f'{args_gen_map} --map-scale {map_scale}'
            args_gen_queue.append(f'{args_gen_map} --out-dir {args.output_dir}/no-sampling --methods CDT MinLT')
            for sample_dist in sample_dists:
                args_gen_sample_dist = f'{args_gen_map}'
                args_gen_sample_dist = f'{args_gen_sample_dist} --max-sample-dist {sample_dist}'
                args_gen_sample_dist = f'{args_gen_sample_dist} --out-dir {args.output_dir}/max-sample-dist-{sample_dist}'
                args_gen_sample_dist = f'{args_gen_sample_dist} --methods MinVT'
                args_gen_queue.append(args_gen_sample_dist)

        print(colored('=== PRINTING ARGUMENTS QUEUE ===', 'green'))
        for args_gen in args_gen_queue:
            print(args_gen)
        print(colored('================================', 'green'))

        time.sleep(5)

        n = len(args_gen_queue)
        for i, args_gen in enumerate(args_gen_queue):
            cmd = f'{args.executable_gen} {args_gen}'
            print(colored(f'{i}/{n}: Running "{cmd}"', 'cyan'))
            os.system(cmd)
            print(colored(f'{i}/{n}: Finished "{cmd}"', 'magenta'))

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

        for args_eval in args_eval_queue:
            cmd = f'{args.executable_eval} {args_eval}'
            print(colored(f'Running "{cmd}"', 'cyan'))
            os.system(cmd)
            print(colored(f'Finished "{cmd}"', 'magenta'))


if __name__ == '__main__':
    main(args_parser.parse_args())
