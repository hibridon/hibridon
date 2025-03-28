#!/usr/bin/env python3
from typing import Optional, List, Tuple
import argparse
from pathlib import Path
import re
"""
zwcheck.py analyzes 
    ! disable-warnings: 

"""


DISABLE_WARNINGS_KEYWORD = 'disable-warnings'  # use this keyword in source file to make zwcheck.py ignore some false psitive warnings caused by the compilation of the line it's on.


def remove_ansi(ansi_string: str):
    # # 7-bit C1 ANSI sequences
    # ansi_escape = re.compile(r'''
    #     \x1B  # ESC
    #     (?:   # 7-bit C1 Fe (except CSI)
    #         [@-Z\\-_]
    #     |     # or [ for CSI, followed by a control sequence
    #         \[
    #         [0-?]*  # Parameter bytes
    #         [ -/]*  # Intermediate bytes
    #         [@-~]   # Final byte
    #     )
    # ''', re.VERBOSE)
    # result = ansi_escape.sub('', sometext)

    # or, without the VERBOSE flag, in condensed form:

    # ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    # result = ansi_escape.sub('', sometext)
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    non_ansi_string = ansi_escape.sub('', ansi_string)
    return non_ansi_string


class FlagWarning():
    type_id: Optional[str]  # eg 'f951'
    severity: Optional[str]  # eg 'warning'
    description: Optional[str]  # eg 'Flag ‘-fno-automatic’ overwrites ‘-frecursive’ implied by ‘-fopenmp’'

    def __str__(self):
        return f'{self.type_id}: {self.description}'


WarningTypeId = str  # eg 'do-subscript'


class Warning():
    compile_option: Optional[str]  # the compilation option that triggered the warning, eg '-Wdo-subscript'
    description: Optional[str]  # eg 'Array reference at (1) out of bounds (0 < 1) in loop beginning at (2)'
    src_file_path: Optional[Path]
    src_line_number: Optional[int]  # the line number in the source code causing the warning
    src_col_number: Optional[int]  # the column number in the source code causing the warning
    src_line: Optional[str]  # eg '    vdif(i-1)=damp*vdif(i-1)  ! disable-warnings:do-subscript'
    is_disabled: bool
    details: List[str]  # the details lines outputted by the warning

    def __init__(self):
        self.compile_option = None
        self.description = None
        self.src_file_path = None
        self.src_line_number = None
        self.src_col_number = None
        self.src_line = None
        self.is_disabled = False
        self.details = []

    def get_type_id(self) -> WarningTypeId:
        assert self.compile_option is not None
        match = re.match(r'^-W(?P<warning_type_id>[^ ]*)$', self.compile_option)
        assert match, f'failed to parse compile option: {self.compile_option}'
        return match['warning_type_id']

    def __str__(self):
        return f'{self.get_type_id()} in {self.src_file_path}:{self.src_line_number}: {self.description} {self.is_disabled}'

    def print_details(self):
        for details_line in self.details:
            print(details_line)


def get_fortran_comments(fortran_src_line: str) -> str:
    line_parts = fortran_src_line.split('!', maxsplit=1)
    comments = ''
    if len(line_parts) > 1:
        comments = line_parts[1]
    return comments


def parse_warning_detail_line(detail_line: str, warning: Warning, disabled_warnings: List[WarningTypeId]):
    # find the culprit source code line number and the possible disable annotation in detail_line
    # print(f'detail_line: {detail_line}')
    warning.details.append(detail_line)
    match = re.match(r'^ *(?P<src_line_number>[0-9]+) *\|(?P<src_line>.*)$', detail_line)
    if match:
        # eg '  186 |     vdif(i-1)=damp*vdif(i-1)  ! disable-warnings:do-subscript'
        warning.src_line_number = match['src_line_number']
        warning.src_line = match['src_line']
        comments = get_fortran_comments(warning.src_line)
        if comments != '':
            # print(comments)
            for part in comments.split(' '):
                match = re.match(r'^' + DISABLE_WARNINGS_KEYWORD + r':(?P<disabled_warnings>.*)$', part)
                if match:
                    for disabled_warning in match['disabled_warnings'].split(','):
                        disabled_warnings.append(disabled_warning)
                        # print(f'disabled warning in {detail_line}')
                    # assert False
    else:
        # Warning: Array reference at (1) out of bounds (0 < 1) in loop beginning at (2) [-Wdo-subscript]
        match = re.match(r'^Warning: (?P<warning_description>[^\[]+)\[(?P<compile_option>[^\]]+)\]$', detail_line)
        # match = re.match(r'^Warning: (?P<warning_description>[^[]+)]\[(?P<compile_option>[^]]+)\]$', non_ansi_line)
        if match:
            # print(match['warning_description'])
            warning.compile_option = match['compile_option']
            warning.description = match['warning_description']
            if warning.get_type_id() in disabled_warnings:
                warning.is_disabled = True


def parse_warnings(make_stdout_file_path: Path) -> Tuple[List[Warning], List[FlagWarning]]:
    flag_warnings: List[FlagWarning] = []
    warnings: List[Warning] = []
    current_warning = None
    with open(make_stdout_file_path, 'r', encoding='utf8') as make_stdout:
        for line in make_stdout:
            # print(line)
            non_ansi_line = remove_ansi(line)
            if re.match(r'^\[ *[0-9]+%\]', non_ansi_line):
                continue
            # f951: Warning: Flag ‘-fno-automatic’ overwrites ‘-frecursive’ implied by ‘-fopenmp’
            match = re.match(r'^(?P<warning_id>f[0-9]+): *(?P<severity>[^ ]+): (?P<description>.*)$', non_ansi_line)
            if match:
                flag_warning = FlagWarning()
                flag_warning.type_id = match['warning_id']
                flag_warning.severity = match['severity'].lower()
                assert flag_warning.severity in ['warning']
                flag_warning.description = match['description']
                flag_warnings.append(flag_warning)
            else:
                match = re.match(r'^(?P<src_file_path>[^:]*):(?P<src_line_number>[0-9]+):(?P<src_col_number>[0-9]+):', non_ansi_line)
                if match:
                    # eg '/home/graffy/work/hibridon/hibridon/tests/arno_bound/pot_arno_ccsdt.F90:186:24:'
                    current_warning = Warning()
                    disabled_warnings: List[WarningTypeId] = []
                    current_warning.src_file_path = Path(match['src_file_path'])
                    current_warning.src_line_number = Path(match['src_line_number'])
                    # assert False
                    warnings.append(current_warning)
                else:
                    if current_warning is not None:
                        # we expect non_ansi_line to be related to current_warning
                        parse_warning_detail_line(non_ansi_line, current_warning, disabled_warnings)
    return warnings, flag_warnings


def check_warnings(make_stdout_file_path: Path, show_warning_details: bool) -> bool:
    warnings, flag_warnings = parse_warnings(make_stdout_file_path)

    non_disabled_warnings = [warning for warning in warnings if not warning.is_disabled]

    ignored_flag_warnings = {
        'f951'  # f951: Warning: Flag ‘-fno-automatic’ overwrites ‘-frecursive’ implied by ‘-fopenmp’
    }

    non_disabled_flag_warnings = [warning for warning in flag_warnings if warning.type_id not in ignored_flag_warnings]

    print(f'{len(non_disabled_warnings)}/{len(warnings)} non disabled source code warnings:')
    warning_index = 1
    for warning in non_disabled_warnings:
        print(f'  {warning_index}: {warning}')
        if show_warning_details:
            warning.print_details()
            print(f'note: if this is a false positive, you can make zwcheck.py ignore this warning by adding the comment "{DISABLE_WARNINGS_KEYWORD}:{warning.get_type_id()}" on the line {warning.src_file_path}:{warning.src_line_number}')
        warning_index += 1
    print(f'{len(non_disabled_flag_warnings)}/{len(flag_warnings)} non disabled compilation flag warnings:')
    warning_index = 1
    for warning in non_disabled_flag_warnings:
        print(f'  {warning_index}: {warning}')

    return len(non_disabled_warnings) + len(non_disabled_flag_warnings)


def main():
    parser = argparse.ArgumentParser(description='analyzes make\'s stdout to detect non disabled warnings')
    parser.add_argument('--make-stdout', type=Path, required=True)
    parser.add_argument('--show-warnings', type=str, default='true')

    args = parser.parse_args()
    show_warnings = {
        'true': True,
        'false': False
    }[args.show_warnings]
    exit_code = 0
    num_remaining_warnings = check_warnings(args.make_stdout, show_warnings)
    if num_remaining_warnings > 0:
        exit_code = 1
    exit(exit_code)


main()
