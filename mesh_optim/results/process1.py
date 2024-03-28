#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import argparse
from termcolor import colored
import datatable as dt

from pathlib import Path

args_parser = argparse.ArgumentParser()
args_parser.add_argument('-i', '--input_data_dir', type=str, default='exp_data', help='Directory with experimental data.')
args_parser.add_argument('-o', '--output_csv_file', type=str, default='exp_data/tuning4_weights.csv', help='Output table in CSV format.')
args_parser.add_argument('--append', default=False, action='store_true', help="Append to the end of the CSV table.")
args_parser.add_argument('--only_info_files', default=False, action='store_true', help="This script will process only mesh_info.json files.")


def process_results(data_dir: str, output_file: str, only_info_files: bool, append: bool):
    print(f'PROCESSING DATA: {data_dir} ...')
    all_files = []
    if not only_info_files:
        all_files = [str(path.resolve()) for path in Path(data_dir).rglob('mesh_eval_*.json')]
    no_eval_files = len(all_files) == 0
    if no_eval_files:
        all_files = [str(path.resolve()) for path in Path(data_dir).rglob('mesh_info.json')]
    n_files = len(all_files)
    print(f'N files: {n_files}')

    cnt = 0
    tab_main = dt.Frame()
    for curr_file in all_files:
        cnt = cnt + 1
        print(f'File {cnt} / {n_files}: {curr_file}')
        path_curr_file = Path(curr_file)
        if not path_curr_file.is_file():
            print(colored(f'File {curr_file} does not exist!', 'red'))
        mesh_info_file = Path(f'{os.path.dirname(path_curr_file)}/mesh_info.json')
        if not mesh_info_file.is_file():
            print(colored(f'File {mesh_info_file} does not exist!', 'red'))
            continue
        tab = dt.Frame([{'exp_data': data_dir}])
        with open(mesh_info_file, 'r') as file:
            data = json.load(file)
            for key, value in data.items():
                tab[key] = value
        if not no_eval_files:
            with open(path_curr_file, 'r') as file:
                data = json.load(file)
                for key, value in data.items():
                    tab[key] = value
        for i in range(tab.nrows):
            if tab[i, :][dt.rowany(dt.math.isna(dt.f[:])), :].nrows == 0:
                tab_main = dt.rbind([tab_main, tab[i, :]])
    tab_main.to_csv(output_file, append=append)


def main(args):
    process_results(args.input_data_dir, args.output_csv_file, args.only_info_files, args.append)


if __name__ == '__main__':
    main(args_parser.parse_args())
