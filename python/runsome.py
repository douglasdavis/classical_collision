#!/usr/bin/env python2

import subprocess

default_settings='--p-radius 1 --t-radius 2 --p-mass 0.25 --t-mass 2 -s 1'

cur = .5
for i in xrange(100):
    subprocess.call("./main "+str(default_settings)+" -l "+str(cur)+" >> data.dat",shell=True)
    cur = cur + 0.025
