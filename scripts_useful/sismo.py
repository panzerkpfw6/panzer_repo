#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy             as np
import os
import sys, argparse, glob
import shutil
import math
import fileinput

from matplotlib.widgets import Slider, Button, RadioButtons
from skimage            import data
from skimage            import io
from array              import array
from inspect            import currentframe as cf

def plot_sismo(file_name):
  if not os.path.exists(file_name):
    raise Exception("'%s' file not found" % (file_name) +
                    " @ %s:%d" % (cf().f_code.co_filename, cf().f_lineno))  
  with open(file_name) as f:
    lines = f.readlines()
    x     =  [line.split()[0] for line in lines]
    y     =  [line.split()[1] for line in lines]
    z     =  [line.split()[2] for line in lines]
    
    y     = np.array(y)
    z     = np.array(z)
    plt.plot(x,y)
    plt.plot(x,z)
    plt.show()

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = 'sismo vizu tool.')
  parser.add_argument('filename', type=str,
                      help='the path to the wave output.')
  args = parser.parse_args(sys.argv[1:])
  try:
    plot_sismo(args.filename)
  except Exception as e:
    print('[SISMO]: %s' % e)
    exit(1)
