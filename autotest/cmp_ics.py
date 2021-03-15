#!/usr/bin/env python2

class HibIntc:
    """I/O of integral cross section."""

    def read_from_ics(self, FILENAME):
        """
        Reading hibridon .ics file for the full matrix.

        FILENAME: path to the .ics file.

        """
        f = open(FILENAME, 'rU')
        # Line 1: date & time
        self.date = f.readline().strip()
        # Line 2: job name
        self.job_name = f.readline().strip()
        # Line 3: pot name
        self.pot_name = f.readline().strip()
        # Line 4: energy and ?
        # TODO: the second parameter
        try:
            self.energy = float(f.readline().split()[0])
        except IndexError:
            print FILENAME, "not successfully read."
        # Line 5, 6: parameters
        # TODO: what are they?
        for i in range(2): f.readline()
        # Line 7: number of levels and available levels
        line = f.readline().split()
        self.N_ALL = int(line[0])
        self.N = int(line[1])
        # Read the i, j of the levels        #TODO
        # TODO: when error?
        self.levels_all = []
        while (len(self.levels_all) != self.N_ALL):
            line = f.readline().split()
            for i in range(len(line) // 2):
                self.levels_all.append((int(line[2 * i]),
                                        int(line[2 * i + 1])))
        # Read the energies of the levels
        # TODO: when error?
        self.level_energies = []
        while (len(self.level_energies) != self.N_ALL):
            line = f.readline().split()
            for number in line:
                self.level_energies.append(float(number) * 219474.6)
        # Get the available energy dictionary.
        count = 0
        self.levels = {}
        for i in range(self.N_ALL):
            if self.level_energies[i] <= self.energy:
                self.levels[self.levels_all[i]] = count
                count += 1
        # TODO: what if count != len?
        # Read the big matrix!
        # TODO: when error?
        self.intcrs = []
        i = 0
        while (i != self.N):
            self.intcrs.append([])
            while (len(self.intcrs[i]) != self.N):
                line = f.readline().split()
                for number in line:
                    self.intcrs[i].append(float(number))
            i += 1
        f.close()
        return

    def __init__(self, input_file):
        self.read_from_ics(input_file)
        return


def compare_ics(file1, file2, tol, eps):
    obj1 = HibIntc(file1)
    obj2 = HibIntc(file2)
    if obj1.N != obj2.N:
        return False
    for i in range(obj1.N):
        for j in range(obj2.N):
            ics1 = obj1.intcrs[i][j]
            ics2 = obj2.intcrs[i][j]
            if abs((ics2 - ics1)) / ics1 > tol and abs(ics1) > eps:
                return False
    return True
