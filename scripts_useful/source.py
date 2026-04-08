#!/usr/bin/env python
import matplotlib.pyplot as plt
import os, sys, argparse, glob
 
def plot_source(file_name):
  if not os.path.exists(file_name):
    raise Exception("'%s' file not found" % (file_name) +
                    " @ %s:%d" % (cf().f_code.co_filename, cf().f_lineno))
  plt.plotfile(file_name, delimiter=' ', cols=(0, 1), names=('time_steps', 'amplitude'), marker='o')
  plt.show()


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = 'source function viewer.')
  parser.add_argument('filename', type=str, help='the path to the source file.')
  args = parser.parse_args(sys.argv[1:])
  try:
    plot_source(args.filename)
  except Exception as e:
    print('[SIMWAVE SOURCE VIEWER]: %s' % e)
    exit(1)
