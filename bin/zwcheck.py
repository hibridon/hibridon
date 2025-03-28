#!/usr/bin/env python3
from typing import Optional, List, Tuple
import argparse
from pathlib import Path
import re
"""
zwcheck.py (zw stands for zero warning) analyzes the output of make built for a cmake project involving compilation. It is designed to detects warnings and fails if they are not explicitely flagged as false positives. Therefore, zwcheck.py can help ensure a zero warning policy.

This postprocessing is a workaround to the fact that gfortran, unlike gcc (via pragma) doesn't have any mechanism to disable a warning on a given source code line

zwcheck.py will ignore warnings explicitly flagged as disabled. To explicitely flag the warnings 'do-subscript' and 'maybe-uninitialized' as disabled on the line 42 of toto.F90, the user simply has to add the following comment to line 42 of toto.F90:

```
! disable-warnings:do-subscript,maybe-uninitialized
```
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

    # handle ansi links
    # A non-xterm extension is the hyperlink, ESC ]8;;link ST
    # eg 'ESC]8;;https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html#index-Wcharacter-truncation^G-Wcharacter-truncation'
    # ESC = 0x1B
    # BEL (^G) = 0x07
    ansi_string_wout_link = re.sub(r'\x1B\]8;;[^\x07]*\x07', '', ansi_string)

    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    non_ansi_string = ansi_escape.sub('', ansi_string_wout_link)
    return non_ansi_string


class FlagWarning():
    type_id: Optional[str]  # eg 'f951'
    severity: Optional[str]  # eg 'warning'
    description: Optional[str]  # eg 'Flag ‘-fno-automatic’ overwrites ‘-frecursive’ implied by ‘-fopenmp’'

    def __str__(self):
        return f'{self.type_id}: {self.description}'


WarningTypeId = str  # eg 'do-subscript'


class Alert():
    _alert_type: Optional[str]  # one of {'warning', 'note'}
    compile_option: Optional[str]  # the compilation option that triggered the warning, eg '-Wdo-subscript'
    description: Optional[str]  # eg 'Array reference at (1) out of bounds (0 < 1) in loop beginning at (2)'
    src_file_path: Optional[Path]
    src_line_number: Optional[int]  # the line number in the source code causing the warning
    src_col_number: Optional[int]  # the column number in the source code causing the warning
    src_line: Optional[str]  # eg '    vdif(i-1)=damp*vdif(i-1)  ! disable-warnings:do-subscript'
    is_disabled: bool
    details: List[str]  # the details lines outputted by the warning
    related_notes: List[Warning]

    def __init__(self):
        self._alert_type = None
        self.compile_option = None
        self.description = None
        self.src_file_path = None
        self.src_line_number = None
        self.src_col_number = None
        self.src_line = None
        self.is_disabled = False
        self.details = []
        self.related_notes = []

    @property
    def alert_type(self) -> str:
        return self._alert_type

    @alert_type.setter
    def alert_type(self, alert_type: str):
        assert alert_type in {'warning', 'note'}
        self._alert_type = alert_type

    def get_type_id(self) -> WarningTypeId:
        assert self.compile_option is not None, f'self = {self}'
        match = re.match(r'^-W(?P<warning_type_id>[^ ]*)$', self.compile_option)
        assert match, f'failed to parse compile option: {self.compile_option}'
        return match['warning_type_id']

    def add_related_note(self, note: Warning):
        assert note.alert_type == 'note'
        self.related_notes.append(note)

    def __str__(self):
        return f'{self.compile_option} in {self.src_file_path}:{self.src_line_number}: {self.description} {self.is_disabled}'

    def print_details(self):
        for details_line in self.details:
            print(details_line)


def get_fortran_comments(fortran_src_line: str) -> str:
    line_parts = fortran_src_line.split('!', maxsplit=1)
    comments = ''
    if len(line_parts) > 1:
        comments = line_parts[1]
    return comments


def parse_warning_detail_line(detail_line: str, alert: Alert, disabled_warnings: List[WarningTypeId]):
    # find the culprit source code line number and the possible disable annotation in detail_line
    # print(f'detail_line: {detail_line}')
    alert.details.append(detail_line)
    match = re.match(r'^ *(?P<src_line_number>[0-9]+) *\|(?P<src_line>.*)$', detail_line)
    if match:
        # eg '  186 |     vdif(i-1)=damp*vdif(i-1)  ! disable-warnings:do-subscript'
        alert.src_line_number = match['src_line_number']
        alert.src_line = match['src_line']
        comments = get_fortran_comments(alert.src_line)
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
            # print(f'compile option line : "{detail_line}"')
            # print(match['warning_description'])
            alert.alert_type = 'warning'
            alert.compile_option = match['compile_option']
            alert.description = match['warning_description']
            if alert.get_type_id() in disabled_warnings:
                alert.is_disabled = True
        else:
            # handle notes such as in
            # /opt/ipr/cluster/work.local/graffy/hibridon/issue174/hibridon.git/lib/hitensor.F90:2077:22:

            #  2077 |   if (iabsty.eq.4 .or. ibasty.eq.19) then
            #       |                      ^
            # Warning: ‘iabsty’ may be used uninitialized [-Wmaybe-uninitialized]
            # /opt/ipr/cluster/work.local/graffy/hibridon/issue174/hibridon.git/lib/hitensor.F90:2077:12:

            #  2077 |   if (iabsty.eq.4 .or. ibasty.eq.19) then
            #       |            ^
            # note: ‘iabsty’ was declared here

            match = re.match(r'^note: (?P<note_description>.*)$', detail_line)
            if match:
                alert.alert_type = 'note'
                alert.description = match['note_description']


def parse_warnings(make_stdout_file_path: Path) -> Tuple[List[Alert], List[FlagWarning]]:
    flag_warnings: List[FlagWarning] = []
    warnings: List[Alert] = []
    current_alert = None
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
                    if current_alert is not None:
                        if current_alert.alert_type == 'warning':
                            warnings.append(current_alert)
                            current_alert = None
                        elif current_alert.alert_type == 'note':
                            warnings[-1].add_related_note(current_alert)
                        else:
                            # assert False
                            pass
                    current_alert = Alert()
                    disabled_warnings: List[WarningTypeId] = []
                    current_alert.src_file_path = Path(match['src_file_path'])
                    current_alert.src_line_number = Path(match['src_line_number'])
                    # assert False
                else:
                    if current_alert is not None:
                        # we expect non_ansi_line to be related to current_warning
                        parse_warning_detail_line(non_ansi_line, current_alert, disabled_warnings)
    if current_alert.alert_type == 'warning':
        warnings.append(current_alert)
        current_alert = None
    elif current_alert.alert_type == 'note':
        warnings[-1].add_related_note(current_alert)
    else:
        # assert False
        pass
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
        print(f'  warning {warning_index}: {warning}')
        if show_warning_details:
            warning.print_details()
            note_index = 1
            for note in warning.related_notes:
                print(f'  warning {warning_index} note {note_index}: {note.src_file_path}:{note.src_line_number} :')
                note.print_details()
                note_index += 1
            print(f'note: if warning {warning_index} is a false positive, you can make zwcheck.py ignore this warning by adding the comment "{DISABLE_WARNINGS_KEYWORD}:{warning.get_type_id()}" on the line {warning.src_file_path}:{warning.src_line_number}')
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
