#!/usr/bin/env python2

import subprocess

cur = 0.001
for i in xrange(100):
    subprocess.call("./main -e "+str(cur)+" >> data.dat",shell=True)
    cur = cur + 0.0002
