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

def remove_keymap_conflicts(new_keys_set):
  for prop in plt.rcParams:
    if prop.startswith('keymap.'):
      keys = plt.rcParams[prop]
      remove_list = set(keys) & new_keys_set
      for key in remove_list:
        keys.remove(key)

def previous_slice(ax):
  """Go to the previous slice."""
  volume = ax.volume
  ax.index = (ax.index - 1) % volume.shape[0]  # wrap around using %
  ax.images[0].set_array(volume[ax.index])

def next_slice(ax):
  """Go to the next slice."""
  volume = ax.volume
  ax.index = (ax.index + 1) % volume.shape[0]
  ax.images[0].set_array(volume[ax.index])

def process_key(event):
  fig = event.canvas.figure
  ax  = fig.axes[0]
  if event.key == 'j':
    previous_slice(ax)
  elif event.key == 'k':
    next_slice(ax)
  fig.canvas.draw()

def slice_keyb_viewer(volume, axis=0, **kwargs):
  remove_keymap_conflicts({'j', 'k'})
  fig, ax   = plt.subplots()
  ax.volume = volume
  ax.set_title("stencil")
  ax.index  = 0 #volume.shape[0]
  ax.imshow(volume[ax.index])#, cmap='nipy_spectral')
  fig.canvas.mpl_connect('key_press_event', process_key)
  plt.show()

def slice_slid_viewer(cube, axis=0, **kwargs):
  """
  Display a 3d ndarray with a slider to move along the third dimension.
  Extra keyword arguments are passed to imshow
  """
  # check dim
  if not cube.ndim == 3:
    raise Exception("cube should be an array with dim == 3" +
                    " @ %s:%d" % (cf().f_code.co_filename, cf().f_lineno))

  # generate figure
  fig = plt.figure()
  ax  = plt.subplot(111)
  ax.set_title("stencil")
  fig.subplots_adjust(left=0.25, bottom=0.25)

  # select first image
  s  = [slice(0, 1) if i == axis else slice(None) for i in range(3)]
  im = cube[s].squeeze()

  # display image
  l  = ax.imshow(im, **kwargs)

  # define slider
  axcolor = 'lightgoldenrodyellow'
  ax      = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)

  slider  = Slider(ax, 'Axis %i index' % axis, 0, cube.shape[axis] - 1,
                   valinit=0, valfmt='%i')
  def update(val):
    ind = int(slider.val)
    s   = [slice(ind, ind + 1) if i == axis else slice(None) for i in range(3)]
    im  = cube[s].squeeze()
    l.set_data(im, **kwargs)
    fig.canvas.draw()
  slider.on_changed(update)
  plt.show()
  
def slice_auto_viewer(volume, axis=0, **kwargs):
  fig, ax   = plt.subplots()
  ax.volume = volume
  ax.set_title("stencil")
  ax.index  = 0
  l = ax.imshow(volume[0])
  fig.canvas.draw() 
  for i in range(1, volume.shape[0]-1):
    ax.index  = i
    l.set_data(volume[i])
    fig.canvas.draw()
    plt.pause(0.001)
 
def vizu(file_name, nx, ny, keyb, slid, auto):
  if not os.path.exists(file_name):
    raise Exception("'%s' file not found" % (file_name) +
                    " @ %s:%d" % (cf().f_code.co_filename, cf().f_lineno))
  tab       = np.fromfile(file_name, dtype=np.float32, count=-1)
  nb_slices = tab.size//(nx*ny)
  print ("... nb entries/slice  : %d" % (nx*ny))
  print ("... nb slices         : %d" % (nb_slices))
  a = np.reshape(tab, (nx,ny,nb_slices), order='C')
  b = a.transpose([2,1,0])
  if keyb == True:   slice_keyb_viewer(b, 0)
  elif slid == True: slice_slid_viewer(b, 0)
  else:              slice_auto_viewer(b, 0)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = 'stencil vizu tool.')
  parser.add_argument('filename', type=str,
                      help='the path to the wave output.')
  parser.add_argument('--nx', type=int, default=100,
                      help='the fastest dimension count.')
  parser.add_argument('--ny', type=int, default=100,
                      help='the second fastest dimension count.')
  parser.add_argument('--keyb', action='store_true',
                      help='keyboard controlled viewer.')
  parser.add_argument('--slid', action='store_true',
                      help='viewer with a slider controlled by mouse.')
  parser.add_argument('--auto', action='store_true',
                      help='automatic viewer.')

  args = parser.parse_args(sys.argv[1:])
  try:
    vizu(args.filename, args.nx, args.ny, args.keyb, args.slid, args.auto)
  except Exception as e:
    print('[SIMWAVE VIZU]: %s' % e)
    exit(1)
