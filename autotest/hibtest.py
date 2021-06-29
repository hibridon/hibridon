#!/usr/bin/env python2

import sys

from test_exec import HibTest
from pathlib import Path


def run_tests(tests: list(str), hib_dir: Path):
    ncmp = {}
    fail_count = 0
    test_count = 0
    tests = []
    for test in tests:
        groupfile = hib_dir / "autotest" / ("testgroup_" + test)
        if groupfile.exists():
            for line in open(groupfile, "r"):
                stripped_line = line.strip()
                if stripped_line:
                    tests.append(stripped_line)
        else:
            tests.append(test)
    for test in tests:
        conf_file_path = hib_dir / "autotest" / test / "hibautotest.conf"
        if not conf_file_path.exists():
            print("** WARNING: Test configuration file ", conf_file_path, " does not exist.")
            continue
        task = HibTest(test, hib_dir)
        nc = task.execute_test()
        test_count += 1
        if nc is None:
            fail_count += 1
            continue
        for key in nc:
            ncmp[key] = nc[key]
    print
    print("%d/%d test(s) passed, %d file(s) not "
          "checked" % (test_count - fail_count, test_count, len(ncmp)))
    print
    return ncmp


def run_diff(ncmp):
    # try:
    #     os.unlink("test_diff.test")
    # except OSError:
    #     pass
    # for f1 in ncmp:
    #     f2 = ncmp[f1]
    #     cmd1 = "echo " + f1 + " " + f2 + " >> test_diff.test"
    #     cmd2 = "diff " + f1 + " " + f2 + " >> test_diff.test"
    #     os.system(cmd1)
    #     os.system(cmd2)
    for f1 in ncmp:
        print("   ", f1)
    return


def main():
    hib_dir = Path(sys.argv[1])
    if len(sys.argv) == 2:
        tests = ["regular"]
    else:
        tests = sys.argv[2:]
    ncmp = run_tests(tests, hib_dir)
    if len(ncmp) > 0:
        print("Please check manually the following file(s):")
        run_diff(ncmp)
    return


if __name__ == "__main__":
    main()
