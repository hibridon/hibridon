#!/usr/bin/env python2

def read_number_lines(f, ignore_entry):
    l = []
    for line in open(f, "rU"):
        try:
            num = []
            strs = line.split()
            for s in strs:
                if s in ignore_entry:
                    continue
                s.replace("D", "E")
                num.append(float(s))
            if num:
                l.append(num)
        except ValueError:
            pass
    return l


def compare_array(arr1, arr2, tol, eps, ignore_sign):
    if len(arr1) != len(arr2):
        return False
    for i in range(len(arr1)):
        n1 = arr1[i]
        n2 = arr2[i]
        if ignore_sign:
            n1 = abs(n1)
            n2 = abs(n2)
        if n1 != 0 and abs(n1 - n2) / n1 > tol:
            if abs(n1) > eps:
                return False
    return True


def compare_number(file1, file2, tol, eps, ignore_sign=False, ignore_entry=[]):
    l1 = read_number_lines(file1, ignore_entry)
    l2 = read_number_lines(file2, ignore_entry)
    if len(l1) != len(l2):
        return False
    for i in range(len(l1)):
        if not compare_array(l1[i], l2[i], tol, eps, ignore_sign):
            return False
    return True
