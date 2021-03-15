#!/usr/bin/env python2

import os
import subprocess
import tempfile
import shutil
import sys
import ConfigParser

from cmp_ics import compare_ics
from cmp_number import compare_number


class HibTest:
    def __init__(self, test_name, hib_dir):
        self.hib_dir = hib_dir
        # Read test config file
        self.folder = os.path.join(hib_dir, "autotest", test_name)
        config_file = os.path.join(self.folder, "hibautotest.conf")
        self.config = ConfigParser.ConfigParser()
        self.config.read(config_file)
        return
    
    def make_executable(self, cmdlist):
        proc = subprocess.Popen(
            cmdlist, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return proc.wait()

    def copy_file(self, filelist):
        for source in filelist:
            fsrc = os.path.join(self.tmpd, source)
            fdes = os.path.join(self.folder, source + ".test")
            try:
                shutil.copyfile(fsrc, fdes)
            except IOError:
                return -1
        return 0

    def compare_file(self, filelist, tol, eps):
        for source in filelist:
            fstd = os.path.join(self.folder, source)
            fexe = fstd + ".test"
            if fstd.endswith("ics"):
                if not compare_ics(fstd, fexe, tol, eps):
                    return -1
            elif fstd.endswith("dcs") or fstd.endswith("xsc") or \
                    fstd.endswith("evl") or fstd.endswith("flx") or \
                    fstd.endswith("mcs") or fstd.endswith("tcb") or \
                    fstd.endswith("tcs") or fstd.endswith("psc") or \
                    fstd.endswith("xxsc") or fstd.endswith("pcs") or \
                    fstd.endswith("hfx") or fstd.endswith("xms") or \
                    fstd.endswith("ppb") or fstd.endswith("trn"):
                if not compare_number(fstd, fexe, tol, eps):
                    return -1
            elif fstd.endswith("psi"):
                if not compare_number(fstd, fexe, tol, eps, ignore_sign=True):
                    return -1
            elif fstd.endswith("stmix"):
                if not compare_number(fstd, fexe, tol, eps, ignore_entry=["S", "T"]):
                    return -1
            else:
                self.ncmp[fstd] = fexe
        return 0;

    def execute_section(self, section):
        print "Running test for", self.config.get(section, "title"),
        sys.stdout.flush()
        # Make hibridon executable
        cmdlist = self.config.get(section, "makecmd").split()
        if self.make_executable(cmdlist):
            print "... \033[31mFAILED\033[0m:"
            print "    Unable to make executable"
            return None
        hibbin = self.config.get(section, "exec").split()[0]
        bin_path1 = os.path.join(self.hib_dir, "bin", hibbin)
        bin_path2 = os.path.join(self.hib_dir, "bin/progs", hibbin)
        if not os.path.isfile(bin_path1):
            if os.path.isfile(bin_path2):
                shutil.copyfile(bin_path2, bin_path1)
                subprocess.call(["chmod", "+x", bin_path1])
                os.unlink(bin_path2)
            else:
                print "... \033[31mFAILED\033[0m:"
                print "    Unable to make executable"
                return None
        # Copy input files
        inpfiles = self.config.get(section, "input").split()
        for inpfile in inpfiles:
            shutil.copyfile(os.path.join(self.folder, inpfile),
                            os.path.join(self.tmpd, inpfile))
        # Copy pot data files
        potdata_dir = os.path.join(self.hib_dir, "bin/progs/potdata")
        try:
            potdatas = self.config.get(section, "potdata").split()
        except ConfigParser.NoOptionError:
            potdatas = []
        if potdatas and os.path.isdir(potdata_dir):
            potdata_dir_tmp = os.path.join(self.tmpd, "potdata")
            if not os.path.isdir(potdata_dir_tmp):
                os.makedirs(potdata_dir_tmp)
            for potdata in potdatas:
                shutil.copyfile(
                    os.path.join(potdata_dir, potdata),
                    os.path.join(self.tmpd, "potdata", potdata))
        # Execute hibridon
        os.chdir(self.tmpd)
        cmd = self.config.get(section, "exec")
        os.system(cmd)
        os.chdir(self.hib_dir)
        if os.path.isfile(bin_path1):
            os.unlink(bin_path1)
        # Fetch output files
        outputfiles = self.config.get(section, "output").split()
        if self.copy_file(outputfiles) != 0:
            print "... \033[31mFAILED\033[0m:"
            print "    Unable to fetch output files, check " + self.tmpd
            return None
        # Compare files
        self.ncmp = {}
        try:
            tolerance = float(self.config.get(section, "tolerance"))
        except ConfigParser.NoOptionError:
            tolerance = 1.0e-4
        try:
            epsilon = float(self.config.get(section, "epsilon"))
        except ConfigParser.NoOptionError:
            epsilon = 1.0e-10
        if self.compare_file(outputfiles, tolerance, epsilon) != 0:
            print "... \033[31mFAILED\033[0m:"
            print "    An output file does not match standard file"
            return None
        if len(self.ncmp) > 0:
            print "... \033[92mPASSED\033[0m (%d/%d NOT CHECKED)" % (len(self.ncmp),
                                                      len(outputfiles))
        else:
            print "... \033[92mPASSED\033[0m (ALL CHECKED)"
        return self.ncmp

    def execute_test(self):
        ncmp = {}
        self.tmpd = tempfile.mkdtemp()
        for section in self.config.sections():
            nc = self.execute_section(section)
            if nc is None:
                return None
            for key in nc:
                ncmp[key] = nc[key]
        shutil.rmtree(self.tmpd)
        return ncmp
