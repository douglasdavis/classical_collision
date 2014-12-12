#!/usr/bin/env python2

import subprocess

cur = .5
for i in xrange(100):
    subprocess.call("./main -c 6 -l "+str(cur)+" >> data.dat",shell=True)
    cur = cur + 0.025
