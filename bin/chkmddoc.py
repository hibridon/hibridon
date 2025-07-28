#!/usr/bin/env python3
# from typing import Optional, List, Tuple
import argparse
from pathlib import Path
import re
import sys
# from logging import error
# import subprocess


def check_link(target: str, link: str, file_path: Path, line_number: int) -> int:
    num_errors = 0
    if re.match(r'^[.a-zA-Z/\-0-9]+$', target):
        # print(f'{file_path}:{line_number}: {target} looks like a local file path')
        target_as_path = Path(target)
        if not target_as_path.suffix != '':
            target_as_path = target_as_path.with_suffix('.md')
        if not target_as_path.is_absolute():
            target_as_path = file_path.parent / target_as_path
        if not target_as_path.exists():
            print(f'{file_path}:{line_number}: {target_as_path} file not found')
            num_errors += 1
    else:
        # print(f'{file_path}:{line_number}: {target} doesn\'t look like a file path')
        pass
    return num_errors


def check_markdown_file(md_file_path: Path) -> int:
    num_errors = 0
    md_file_path = md_file_path.absolute()
    assert md_file_path.exists()
    with open(md_file_path, 'r', encoding='utf8') as md_file:
        line_index = 0
        for line in md_file:
            line_index += 1
            matches = re.findall(r'\[(?P<link>[^]]+)\]\((?P<target>[^)]*)\)', line)
            for m in matches:
                link = m[0]  # m.group('link')
                target = m[1]  # m.group('target'),
                num_errors += check_link(target, link, md_file_path, line_number=line_index)
    return num_errors


def check_markdown_doc(md_doc_root: Path) -> int:
    num_errors = 0
    for file in md_doc_root.rglob('**/*.md'):
        num_errors += check_markdown_file(file)
    return num_errors


def main():
    parser = argparse.ArgumentParser(description='performs checks (dead links) on a markdown documentation')
    parser.add_argument('--md-doc-root', type=Path, required=True, help='path to the root folder containing the documentation in markdown form')
    parser.epilog = 'example: chkmddoc.py --md-doc-root=./docs'

    args = parser.parse_args()
    num_errors = check_markdown_doc(args.md_doc_root)
    exit_code = 0 if num_errors == 0 else 1
    sys.exit(exit_code)


main()
