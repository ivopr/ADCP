#!/usr/bin/env python3
################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
##
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2019
##
################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2019
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/utils/runADFR.py,v 1.36 2017/11/15 00:41:43 sanner Exp $
#
# $Id: runADFR.py,v 1.36 2017/11/15 00:41:43 sanner Exp $
#
#

import os
import sys
import numpy
import platform
import datetime
import shutil
import random
import subprocess
import zipfile
import multiprocessing

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')


class runADCP:

    def myprint(self, s, newline=True):
        sys.stdout.write(s)
        if newline:
            sys.stdout.write('\n')

    def myexit(self):
        if self.targetFile is not None:
            print("clean up unzipped map files")
            try:
                shutil.rmtree('./tmp_%s' % self.jobName)
            except OSError:
                pass
            for element in ['C', 'A', 'SA', 'N', 'NA', 'OA', 'HD', 'd', 'e']:
                mapfile = 'rigidReceptor.%s.map' % element
                if os.path.isfile(mapfile):
                    os.remove(mapfile)
            if os.path.isfile('transpoints'):
                os.remove('transpoints')
            if os.path.isfile('translationPoints.npy'):
                os.remove('translationPoints.npy')
            if os.path.isfile('con'):
                os.remove('con')
        sys.exit(0)

    def __init__(self):
        self.ncpu = multiprocessing.cpu_count()
        system_info = platform.uname()
        _platform = system_info[0]

        if _platform == 'Windows':
            self.shell = False
            self._argv = ['./adcp -t 2']
        else:
            self.shell = True
            self._argv = ['./adcp -t 2']

        self.completedJobs = 0
        self.numberOfJobs = 0
        self.outputBaseName = None
        self.jobName = 'NoName'
        self.targetFile = None

    def __call__(self, **kw):
        seed = None
        rncpu = None
        nbRuns = kw.get('nbRuns', 50)
        numSteps = kw.get('numSteps', 2500000)
        jobName = kw.get('jobName', 'NoName')
        dryRun = kw.get('dryRun', False)
        overwriteFiles = kw.get('overwriteFiles', False)
        sequence = kw.get('sequence', None)
        input_pdb = kw.get('input', None)
        cyclic = kw.get('cyclic', False)
        cystein = kw.get('cystein', False)
        seedValue = kw.get('seedValue', -1)
        maxCores = kw.get('maxCores', None)
        targetFile_param = kw.get('target', None)

        rncpu = maxCores

        if rncpu is None:
            ncores = self.ncpu
            self.myprint('Detected %d cores, using %d cores' % (self.ncpu, ncores))
        else:
            assert rncpu > 0, "ERROR: maxCores a positive number, got %d" % rncpu
            ncores = min(self.ncpu, rncpu)
            self.myprint('Detected %d cores, request %d cores, using %d cores' % (self.ncpu, rncpu, ncores))

        if nbRuns is not None:
            self.nbRuns = nbRuns

        self.numberOfJobs = nbRuns
        self._jobStatus = [None] * nbRuns

        seed = seedValue

        if seed is None or seed == -1:
            seed = str(random.randint(1, 999999))
        else:
            seed = str(seed)

        ramaprob_path = os.path.join(DATA_DIR, 'ramaprob.data')
        if not os.path.isfile(ramaprob_path):
            ramaprob_path = 'ramaprob.data'
        if not os.path.isfile(ramaprob_path):
            print("ERROR: cannot find probability data for ramachandra plot")
            self.myexit()

        if jobName is not None:
            self.jobName = jobName

        self.targetFile = targetFile = targetFile_param
        if targetFile is None and not os.path.isfile("transpoints"):
            print("ERROR: no receptor files found")
            self.myexit()
        elif targetFile is not None:
            with zipfile.ZipFile(targetFile, 'r') as zip_ref:
                zip_ref.extractall('./tmp_%s/' % jobName)
            for element in ['C', 'A', 'SA', 'N', 'NA', 'OA', 'HD', 'd', 'e']:
                try:
                    shutil.copy(os.path.join('./tmp_%s/' % jobName, targetFile[:-4],
                                             'rigidReceptor.%s.map' % element), os.getcwd())
                except IOError:
                    print("WARNING: cannot locate map file for %s" % element)
            shutil.copy(os.path.join('./tmp_%s/' % jobName, targetFile[:-4],
                                     'translationPoints.npy'), os.getcwd())
            ttt = numpy.load('translationPoints.npy', allow_pickle=True)
            with open('transpoints', 'w+') as fff:
                fff.write('%s\n' % len(ttt))
                numpy.savetxt(fff, ttt, fmt='%7.3f')
        else:
            for element in ['C', 'A', 'SA', 'N', 'NA', 'OA', 'HD', 'd', 'e']:
                if not os.path.isfile("rigidReceptor.%s.map" % element):
                    print("WARNING: cannot locate map file rigidReceptor.%s.map" % element)

        with open('con', 'w+') as fff:
            fff.write('1\n')

        for i in range(nbRuns):
            fname = '%s_%d.pdb' % (jobName, i + 1)
            if os.path.isfile(fname):
                if not overwriteFiles:
                    print("ERROR: output file exists %s" % fname)
                    self.myexit()
                else:
                    print("Warning: overwriting output file %s" % fname)

        argv = list(self._argv)

        if sequence is None:
            if input_pdb is None or input_pdb[-3:] != 'pdb':
                print("ERROR: no input for peptide found")
                self.myexit()
            else:
                argv.append('-f')
                argv.append('%s' % input_pdb)
        else:
            argv.append('%s' % sequence)

        argv.append('-r')
        if numSteps is not None:
            pass
        argv.append('1x%s' % numSteps)

        ADCPDefaultOptions = "-p Bias=NULL,external=5,con,1.0,1.0"
        if cyclic:
            ADCPDefaultOptions += ",external2=4,con14,1.0,1.0"
        if cystein:
            ADCPDefaultOptions += ",SSbond=80,2.2,20,0.5"
        ADCPDefaultOptions += ",Opt=1,0.25,0.75,0.0"
        argv.append(ADCPDefaultOptions)

        argv.extend(['-s', '-1', '-o', jobName, ' '])

        from time import time, sleep
        t0 = time()
        runStatus = [None] * nbRuns
        runEnergies = [999.] * nbRuns
        procToRun = {}
        nbStart = 0
        nbDone = 0

        self.myprint("Performing search (%d ADCP runs with %d steps each) ..." %
                     (nbRuns, numSteps))
        print("0%   10   20   30   40   50   60   70   80   90   100%")
        print("|----|----|----|----|----|----|----|----|----|----|")

        jobNum = 1
        for jobNum in range(1, min(nbRuns, ncores) + 1):
            if seed == "-1":
                argv[-4] = str(random.randint(1, 999999))
            else:
                argv[-4] = str(int(seed) + jobNum - 1)
            argv[-2] = '%s_%d.pdb' % (jobName, jobNum)
            argv[-1] = '> %s_%d.out 2>&1' % (jobName, jobNum)
            if dryRun:
                print('/n*************** command ***************************\n')
                print(' '.join(argv))
                print()
                self.myexit()

            process = subprocess.Popen(' '.join(argv),
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       bufsize=1, shell=self.shell, cwd=os.getcwd())
            procToRun[process] = jobNum - 1
            nbStart += 1

        while nbDone < nbRuns:
            finished = []
            for proc, jnum in list(procToRun.items()):
                if proc.poll() is not None:
                    if proc.returncode != 0:
                        runStatus[jnum] = ('Error', '%s%04d' % (jobName, jnum + 1))
                        status = 'FAILED'
                    else:
                        status = 'OK'
                        error = ''
                        runStatus[jnum] = ('OKAY', '%s%04d' % (jobName, jnum + 1))
                        with open('%s_%d.out' % (jobName, jnum + 1)) as f:
                            lines = f.readlines()
                        for ln in lines:
                            if ln.startswith('best target energy'):
                                runEnergies[jnum] = float(ln.rstrip().split()[3])

                    nbDone += 1
                    del procToRun[proc]

                    self._jobStatus[jnum] = 2
                    self.completedJobs += 1
                    percent = float(self.completedJobs) / self.numberOfJobs
                    sys.stdout.write('%s\r' % ('*' * int(50 * percent)))
                    sys.stdout.flush()

                    if nbStart < nbRuns:
                        jobNum += 1
                        if seed == "-1":
                            argv[-4] = str(random.randint(1, 999999))
                        else:
                            argv[-4] = str(int(seed) + jobNum - 1)
                        argv[-2] = '%s_%d.pdb' % (jobName, jobNum)
                        argv[-1] = '> %s_%d.out 2>&1' % (jobName, jobNum)
                        try:
                            os.remove(argv[-1])
                        except OSError:
                            pass

                        process = subprocess.Popen(' '.join(argv),
                                                   stdout=subprocess.PIPE,
                                                   stderr=subprocess.PIPE,
                                                   bufsize=1, shell=self.shell, cwd=os.getcwd())
                        procToRun[process] = jobNum - 1
                        nbStart += 1
            sleep(1)

        dt = time() - t0
        h, m, s = str(datetime.timedelta(seconds=dt)).split(':')
        self.myprint('Docking performed in %.2f seconds, i.e. %s hours %s minutes %s seconds ' %
                     (dt, h, m, s))

        sort_index = numpy.argsort(runEnergies)
        for i in range(5):
            self.myprint('No. %d energy found is %3.1f kcal/mol at %s_%d.pdb ' %
                         (i + 1, runEnergies[sort_index[i]] * 0.59219, jobName, sort_index[i] + 1))
        self.myexit()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='AutoDock CrankPep',
                                     usage="usage: python %(prog)s -s GaRyMiChEL -t rec.trg -o output")

    parser.add_argument("-s", "--sequence", dest="sequence")
    parser.add_argument("-i", "--input", dest="input")
    parser.add_argument("-t", "--target", dest="target")
    parser.add_argument("-n", "--numSteps", type=int,
                        default=2500000, dest="numSteps")
    parser.add_argument("-N", "--nbRuns", type=int,
                        default=50, dest="nbRuns")
    parser.add_argument("-c", "--maxCores", type=int, dest="maxCores")
    parser.add_argument("-o", "--jobName", dest="jobName")
    parser.add_argument(
        '-y', "--dryRun", dest="dryRun", action="store_true",
        default=False,
        help="print the first adcp command line and exit")
    parser.add_argument(
        '-cyc', "--cyclic", dest="cyclic", action="store_true",
        default=False,
        help="option for cyclic peptide through backbone")
    parser.add_argument(
        '-cys', "--cystein", dest="cystein", action="store_true",
        default=False,
        help="option for cyclic peptide through CYS-S-S-CYS")
    parser.add_argument(
        '-O', "--overwriteFiles", dest="overwriteFiles",
        action="store_true", default=False,
        help="overwrite existing output files silently")
    parser.add_argument(
        '-S', "--seed", dest="seedValue", type=int, default=-1,
        help="seed for random number generator")
    parser.add_argument(
        '--version', action='version', version='%(prog)s 0.1')
    kw = vars(parser.parse_args())

    runner = runADCP()
    runner(**kw)
