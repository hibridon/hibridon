#!/usr/bin/env python3
# from typing import Optional, List, Tuple
import argparse
from pathlib import Path
# import re
from logging import error
import subprocess


def remove_time_dependent_lines(input_path: Path, output_path: Path):
    ignored_patterns = []  # patterns containing variable numbers
    ignored_patterns.append('HIBRIDON SCATTERING CODE V ')
    ignored_patterns.append('SYS: ')  # eg 'SYS: Darwin 20.5.0 x86_64'
    # darwin info patterns
    ignored_patterns.append('Model Name: ')  # eg 'Model Name: MacBook Pro'
    ignored_patterns.append('Model Identifier: ')  # eg 'Model Identifier: MacBookPro14,3'
    ignored_patterns.append('Processor Name: ')  # eg 'Processor Name: Quad-Core Intel Core i7'
    ignored_patterns.append('Number of Processors: ')  # eg 'Number of Processors: 1'
    ignored_patterns.append('Processor Speed: ')
    ignored_patterns.append('Total Number of Cores: ')  # eg 'Total Number of Cores: 4'
    ignored_patterns.append('L2 Cache ')  # eg 'L2 Cache (per Core): 256 KB'
    ignored_patterns.append('L3 Cache: ')  # eg 'L3 Cache: 8 MB'
    ignored_patterns.append('Memory: ')  # eg 'Memory: 16 GB'
    # linux info patterns
    ignored_patterns.append('Model name:')  # eg 'Model name:                         Intel(R) Core(TM) i5-8350U CPU @ 1.70GHz'
    ignored_patterns.append('CPU MHz:')  # eg 'CPU MHz:                            3176.744'
    ignored_patterns.append('Socket')  # eg 'Socket(s):                          1'
    ignored_patterns.append('Core\\(s\\) per socket:')  # eg 'Core(s) per socket:                 4'
    ignored_patterns.append('L2 cache:')  # eg 'L2 cache:                           1 MiB'
    ignored_patterns.append('L3 cache:')  # eg 'L3 cache:                           6 MiB'
    ignored_patterns.append('MemTotal:')  # eg 'MemTotal:       16255508 kB'

    # time patterns
    ignored_patterns.append('DATE')
    ignored_patterns.append('ELAPSED')
    ignored_patterns.append('WRITTEN')

    # other patterns
    ignored_patterns.append('Maximum number of threads used by OpenMP: ')  # eg 'Maximum number of threads used by OpenMP: 8'

    subprocess.run(f'cat {input_path} | grep -E -v \'({"|".join(ignored_patterns)})\' > {output_path}', shell=True, check=True)


def check_output(out_path: Path, ref_path: Path, backend_path: Path):
    if not backend_path.exists():
        error(f'failed to find the backend program {backend_path}')
    file_ext = out_path.suffix
    if file_ext != ref_path.suffix:
        error(f'{out_path} and {ref_path} are expected to have the same file extension')
    print(file_ext)
    if file_ext == '.stdout':
        preprocessed_out_path = out_path.with_suffix('.tistdout')
        preprocessed_ref_path = (out_path.parent / (out_path.stem + '_ref')).with_suffix('.tistdout')
        remove_time_dependent_lines(out_path, preprocessed_out_path)
        remove_time_dependent_lines(ref_path, preprocessed_ref_path)
    else:
        preprocessed_out_path = out_path
        preprocessed_ref_path = ref_path
    completed_process = subprocess.run([backend_path, preprocessed_ref_path, preprocessed_out_path])
    assert completed_process.returncode == 0


def main():
    parser = argparse.ArgumentParser(description='checks that an hibridon output file\'s results are close enough to the expected reference')
    parser.add_argument('--backend-path', type=Path, required=True, help='path to the check_outputs executable (the backend that performs the actual check)')
    parser.add_argument('--out-path', type=Path, required=True, help='path to the hibridon output file that is checked against the reference file')
    parser.add_argument('--ref-path', type=Path, required=True, help='path to the hibridon output file that is used as a reference file')
    parser.epilog = 'example: chkoutf.py --backend-path=./tests/check_outputs --out-path=./tests/o3ph2/o3ph2_quick.stdout --ref-path=../hibridon.git/tests/o3ph2/o3ph2_quick.stdout'

    args = parser.parse_args()
    exit_code = 0
    check_output(args.out_path, args.ref_path, args.backend_path)
    return exit_code


main()
