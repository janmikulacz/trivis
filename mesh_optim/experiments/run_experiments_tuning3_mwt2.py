#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import argparse
from termcolor import colored
from pathlib import Path

args_parser = argparse.ArgumentParser()
args_parser.add_argument('-g', '--executable_gen', type=str,
                         default='../../build-Release/mesh_optim/generate_meshes',
                         help='Path to the generate_meshes executable.')
args_parser.add_argument('-o', '--output_dir', type=str,
                         default='./outputs/tuning3_mwt2',
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

simple_poly_max_sizes = list(range(50, 600 + 1, 50))
random_seeds = [420 + r for r in list(range(10))]
time_limit = 50
iter_limit = 2000


def main(args):
    args_gen_queue = []

    args_gen_common = ''
    args_gen_common = f'{args_gen_common} --bucket-size 1.0'
    args_gen_common = f'{args_gen_common} --mean-bucket-triangle-count-max -1.0'
    # args_gen_common = f'{args_gen_common} --random-seed-method 420'
    args_gen_common = f'{args_gen_common} --eps-dist-diff-collinear 1e-7'
    args_gen_common = f'{args_gen_common} --max-sample-dist 2'
    args_gen_common = f'{args_gen_common} --save-weights'
    args_gen_common = f'{args_gen_common} --load-weights'
    args_gen_common = f'{args_gen_common} --dont-save-mesh'  # only mesh_info.json will be generated
    args_gen_common = f'{args_gen_common} --methods MinVT'
    args_gen_common = f'{args_gen_common} --mwt-iter-limit {iter_limit}'
    args_gen_common = f'{args_gen_common} --mwt-time-limit {time_limit}'
    for random_seed in random_seeds:
        args_gen_random_seed = f'{args_gen_common}'
        args_gen_random_seed = f'{args_gen_random_seed} --random-seed-method {random_seed}'
        for map_name, map_scale in maps:
            args_gen_map = f'{args_gen_random_seed}'
            args_gen_map = f'{args_gen_map} --map-name {map_name}'
            args_gen_map = f'{args_gen_map} --map-scale {map_scale}'
            for simple_poly_max_size in simple_poly_max_sizes:
                args_gen_max_size = f'{args_gen_map}'
                if simple_poly_max_size <= 0:
                    args_gen_max_size = f'{args_gen_max_size} --out-dir {args.output_dir}/mwt-dp-full'
                    args_gen_queue.append(args_gen_max_size)
                else:
                    args_gen_max_size = f'{args_gen_max_size} --mwt-simple-polygon-max-size {simple_poly_max_size}'
                    args_gen_max_size = f'{args_gen_max_size} --out-dir {args.output_dir}/rand-seed-{random_seed}/mwt-dp-polysize-{simple_poly_max_size}'
                    args_gen_queue.append(args_gen_max_size)

    print(colored('=== PRINTING ARGUMENTS QUEUE ===', 'green'))
    for args_gen in args_gen_queue:
        print(args_gen)
    print(colored('================================', 'green'))

    time.sleep(5)

    for args_gen in args_gen_queue:
        cmd = f'{args.executable_gen} {args_gen}'
        print(colored(f'Running "{cmd}"', 'cyan'))
        os.system(cmd)
        print(colored(f'Finished "{cmd}"', 'magenta'))


if __name__ == '__main__':
    main(args_parser.parse_args())
