#!/usr/bin/env python
import os,sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import math
import segyio
import numpy as np
import subprocess
### import scipy.fft as fft
import pyfftw.interfaces.numpy_fft as fft
import fnmatch
from scipy import signal
#######################################

def load_data(filename, dims):
    """
    Load a 3D seismic data array from a binary file.
    
    Parameters:
        filename (str): Path to the binary file
        dims (tuple or list): Dimensions of the array (nx, ny, nt)
    
    Returns:
        data (numpy.ndarray): 3D array of shape (nx, ny, nt)
    """
    # Unpack dimensions
    nx, ny, nt = dims
    ccnt = nx * ny * nt  # Total number of elements
    
    # Open and read binary file directly into a 1D array
    with open(filename, 'rb') as fileID:
        data = np.fromfile(fileID, dtype=np.float32, count=ccnt)
    
    # Verify the file size
    filesize = os.path.getsize(filename)
    expected_size = ccnt * 4  # 4 bytes per float32
    # if filesize != expected_size:
    #     raise ValueError(f"File size {filesize} does not match expected size {expected_size}")
    
    # Reshape and transpose to (nx, ny, nt)
    # Data in file is ordered as (nt, nx, ny), so reshape to that first
    data = data.reshape(nt, nx, ny)
    # Transpose to (nx, ny, nt)
    data = data.transpose(1, 2, 0)
    
    return data

def load_data_orig(filename, dims):
    # dims is a tuple or list: (nx, ny, nt)
    nx, ny, nt = dims
    ccnt = nx * ny * nt  # Total number of elements
    
    # Open and read binary file
    with open(filename, 'rb') as fileID:
        A = np.fromfile(fileID, dtype=np.float32, count=ccnt)
    
    # Get file size (equivalent to MATLAB's dir().bytes)
    filesize = os.path.getsize(filename)
    
    # Initialize output array with NaNs (like MATLAB's nan(dims))
    data = np.full(dims, np.nan, dtype=np.float32)
    
    # Fill the data array
    tmp = 0  # Index starts at 0 (Python convention)
    for k in range(nt):  # z-dimension (0 to nt-1)
        for i in range(nx):  # x-dimension (0 to nx-1)
            data[i, :, k] = A[tmp:(tmp + ny)]  # Slice ny elements into y-dimension
            tmp += ny
    
    return data

def save_filtered_data(filtered_data, output_fname, dims):
    """
    Save a 3D filtered data array to a binary file in float32 format, trace by trace.
    
    Parameters:
        filtered_data (numpy.ndarray): 3D array of filtered data (shape: nx, ny, nt)
        output_fname (str): Path to the output binary file
        dims (tuple or list): Dimensions of the array (nx, ny, nt)
    """
    # Unpack dimensions
    nx, ny, nt = dims
    
    # Ensure the input is a NumPy array, float32, and matches dims
    filtered_data = np.asarray(filtered_data, dtype=np.float32)
    print(f"Input data shape: {filtered_data.shape}, dims: {dims}")
    if filtered_data.shape != (nx, ny, nt):
        raise ValueError(f"Array shape {filtered_data.shape} does not match dims {dims}")
    
    # Check if output directory exists, create if it doesn't
    output_dir = os.path.dirname(output_fname)
    if output_dir and not os.path.exists(output_dir):
        print(f"Creating directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
    
    # Reorder the array to match the trace-by-trace order: (nx, ny, nt) -> (nt, nx, ny)
    # Then flatten to 1D for writing
    reordered_data = filtered_data.transpose(2, 0, 1).ravel()
    print(f"Reordered data shape: {reordered_data.shape}, total elements: {reordered_data.size}")
    
    # Write the entire array in one go
    try:
        with open(output_fname, 'wb') as fileID:
            reordered_data.tofile(fileID)
        print(f"Successfully wrote file: {output_fname}")
        print(f"File size: {os.path.getsize(output_fname)} bytes")
    except Exception as e:
        print(f"Error writing file {output_fname}: {e}")
        raise
    
    return None

def fkk_filter(data, dt, dx, dy, v_direct, mute_width=0.4, taper_width=0.15, plot_fkk=False):
    """
    Apply 3D F-K-K filtering to remove direct waves from seismic data with tapering.
    
    Parameters:
        data (ndarray): 3D seismic data (n_inline, n_crossline, n_samples)
        dt (float): Time sampling interval (s)
        dx (float): Inline spacing (m)
        dy (float): Crossline spacing (m)
        v_direct (float): Direct wave velocity (m/s)
        mute_width (float): Width of the mute zone in F-K-K domain (fraction of max wavenumber)
        taper_width (float): Width of the taper zone (fraction of max wavenumber)
        plot_fkk (bool): If True, plot the F-K-K spectrum before applying the taper
    
    Returns:
        filtered_data (ndarray): Filtered 3D seismic data
    """
    n_inline, n_crossline, n_samples = data.shape

    # Low-pass filter to avoid aliasing
    k_max = 1 / (2 * dx)  # Maximum wavenumber
    f_max_alias = v_direct * k_max  # Maximum frequency without aliasing
    print(f"Maximum frequency to avoid aliasing: {f_max_alias:.2f} Hz")
    sos = signal.butter(4, f_max_alias, fs=1/dt, btype='low', output='sos')
    data = signal.sosfiltfilt(sos, data, axis=2)

    # 3D FFT to F-K-K domain
    data_fkk = fft.rfftn(data, axes=(0, 1, 2))
    freq_shape = data_fkk.shape[-1]

    # Frequency and wavenumber axes
    freqs = fft.rfftfreq(n_samples, dt)
    kx = fft.fftfreq(n_inline, dx)
    ky = fft.fftfreq(n_crossline, dy)

    # Broadcasting for 3D grid
    kx = kx.reshape(-1, 1, 1)
    ky = ky.reshape(1, -1, 1)
    freqs = freqs.reshape(1, 1, -1)

    # Direct wave in F-K-K domain: |f| = v * sqrt(kx^2 + ky^2)
    k_magnitude = np.sqrt(kx**2 + ky**2)
    max_k = np.max(k_magnitude)
    freq_diff = np.abs(freqs - v_direct * k_magnitude)

    # Plot the F-K-K spectrum if requested
    if plot_fkk:
        amp_spectrum = np.abs(data_fkk)
        amp_spectrum = np.log10(amp_spectrum + 1e-10)

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

        # kx-ky plane at a fixed frequency
        freq_idx = freq_shape // 2
        slice_kx_ky = amp_spectrum[:, :, freq_idx]
        im1 = ax1.imshow(slice_kx_ky, cmap='jet', extent=[ky.min(), ky.max(), kx.max(), kx.min()])
        ax1.set_title(f'kx-ky at f={freqs[0, 0, freq_idx]:.2f} Hz')
        ax1.set_xlabel('ky (1/m)')
        ax1.set_ylabel('kx (1/m)')
        fig.colorbar(im1, ax=ax1, label='Log10(Amplitude)')

        # kx-f plane at a fixed ky
        ky_idx = n_crossline // 2
        slice_kx_f = amp_spectrum[:, ky_idx, :]
        im2 = ax2.imshow(slice_kx_f, cmap='jet', extent=[freqs.min(), freqs.max(), kx.max(), kx.min()])
        ax2.set_title(f'kx-f at ky={ky[0, ky_idx, 0]:.2e} 1/m')
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('kx (1/m)')
        fig.colorbar(im2, ax=ax2, label='Log10(Amplitude)')
        ax2.contour(freq_diff[:, ky_idx, :] < mute_threshold, levels=[0], colors='white', linestyles='--')

        # ky-f plane at a fixed kx
        kx_idx = n_inline // 2
        slice_ky_f = amp_spectrum[kx_idx, :, :]
        im3 = ax3.imshow(slice_ky_f, cmap='jet', extent=[freqs.min(), freqs.max(), ky.max(), ky.min()])
        ax3.set_title(f'ky-f at kx={kx[kx_idx, 0, 0]:.2e} 1/m')
        ax3.set_xlabel('Frequency (Hz)')
        ax3.set_ylabel('ky (1/m)')
        fig.colorbar(im3, ax=ax3, label='Log10(Amplitude)')
        ax3.contour(freq_diff[kx_idx, :, :] < mute_threshold, levels=[0], colors='white', linestyles='--')

        plt.tight_layout()
        plt.show()

    # Define mute and taper thresholds
    mute_threshold = mute_width * v_direct * max_k
    taper_threshold = (mute_width + taper_width) * v_direct * max_k

    # Create a smooth taper
    taper = np.ones_like(data_fkk, dtype=np.float32)
    mute_mask = freq_diff < mute_threshold
    taper_mask = (freq_diff >= mute_threshold) & (freq_diff < taper_threshold)

    taper[mute_mask] = 0.0
    taper[taper_mask] = (freq_diff[taper_mask] - mute_threshold) / (taper_threshold - mute_threshold)

    # Apply the taper
    data_fkk *= taper

    # Inverse 3D FFT
    filtered_data = fft.irfftn(data_fkk, s=(n_inline, n_crossline, n_samples), axes=(0, 1, 2))

    return filtered_data

# Main processing
def main():
    # Parameters
    nx=676;ny=676;nt=800;
    dims=[nx,ny,nt]
    print(f"data: {nx} inlines, {ny} crosslines, {nt} samples")
    v_direct = 1500.0  # Direct wave velocity (m/s, e.g., water)
    dx = 25.0          # Inline spacing (m, adjust based on your data)
    dy = 25.0          # Crossline spacing (m, adjust based on your data)
    dt=0.001
    x=np.arange(nx)*dx
    y=np.arange(ny)*dy
    t=np.arange(nt)*dt*1000

    print(os.getcwd())
    data_dir='../../data/orig_data'

    # save_data_dir='../../data'
    ###
    result = subprocess.run("mkdir ../../data/filtered_real_data", shell=True, capture_output=True, text=True)
    save_data_dir='../../data/filtered_real_data'

    file_list=fnmatch.filter(os.listdir(data_dir), '*sismos*')
    file_list=sorted(file_list)
    file_list=['sismos_20470.raw']
    
    for file in file_list:
        input_file=os.path.join(data_dir,file)
        output_file=os.path.join(save_data_dir,file)
        print("processing ",file)
        # Load data
        data = load_data(input_file,dims)

        # Apply F-K-K filtering
        filtered_data = fkk_filter(data, dt, dx, dy, v_direct)
        # filtered_data=data
        print("F-K-K filtering completed")

        # Save filtered data
        # save_filtered_data(filtered_data, output_file, dims)

        # filtered_data2 = load_data(output_file,dims)

        # plot
        data1=data[310,:,:]
        data2=filtered_data[310,:,:]

        val=1e-3
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns, figure size (width, height)
        # First subplot (like imagesc)
        im1 = ax1.imshow(data1.T, cmap='viridis', aspect='auto',vmin=-val,vmax=val)
        ax1.set_title('orig')
        ax1.set_xlabel('Y, m')
        ax1.set_ylabel('time,msec')
        fig.colorbar(im1, ax=ax1)  # Add colorbar for first subplot
        # Second subplot (like imagesc)
        im2 = ax2.imshow(data2.T, cmap='viridis', aspect='auto',vmin=-val,vmax=val)
        ax2.set_title('filtered')
        ax2.set_xlabel('Y, m')
        ax2.set_ylabel('time,msec')
        fig.colorbar(im2, ax=ax2)  # Add colorbar for second subplot
        # Adjust layout to prevent overlap
        plt.tight_layout()
        # Show the plot
        plt.show()

        aa=1

# print(os.getcwd())
# main();

if __name__ == "__main__":
    main()



#     # backup data and filter sismos
#     result = subprocess.run("pwd", shell=True, capture_output=True, text=True)
#     result = subprocess.run("mkdir ../../data/orig_data", shell=True, capture_output=True, text=True)
#     result = subprocess.run("scp ../../data/*sismos* ../../data/orig_data", shell=True, capture_output=True, text=True)





