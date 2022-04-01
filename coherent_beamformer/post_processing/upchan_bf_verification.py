# Verify whether the same results are seen from the cuFFT in the upchannelizer code

import matplotlib.pyplot as plt
from array import array
import numpy as np

from numpy.fft import fft, ifft

# Open binary file containing input data to cufft for comparison with nump fft
filename = "/datag/users/mruzinda/i/input_h_fft_bf.bin"

# Read file contents: np.fromfile(filename, dtype=float, count=- 1, sep='', offset=0)
contents_float = np.fromfile(filename, dtype=np.float32)

print(len(contents_float))
print(contents_float[0])

# Array dimensions
# 1k mode
#N_time = (4096*1024)/8 # 2^19
#N_coarse = 1 # 4
# 4k mode
N_time = (1024*1024)/8 # 2^17 
N_coarse = 4 # 4
# 32k mode
#N_time = (128*1024)/8 # 2^14
#N_coarse = 32 # 4

N_fine = N_time

N_pol = 2
N_win = 8
N_ant = 64 
N_iq = 2
N_beam = 61

# Reshape array to multidimensional one -> IQ X Polarization X Time samples X Time Windows X Coarse channels X Antenna
x = contents_float[0:(N_coarse*N_win*N_pol*N_ant*N_time*N_iq)].reshape(N_ant,N_coarse,N_win,N_time,N_pol,N_iq)

X_tmp = np.zeros(N_fine)
X = np.zeros(N_win*N_pol*N_ant*N_coarse*N_fine).reshape(N_win,N_pol,N_ant,N_coarse*N_fine)

# Combine coarse and fine channels
for c in range(0,N_coarse):
    for w in range(0,N_win):
        for p in range(0,N_pol):
            for a in range(0,N_ant):
                # FFT
                X_tmp = np.abs(fft(x[a,c,w,0:N_time,p,0] + 1j*x[a,c,w,0:N_time,p,1]))
                # FFT shift
                X[w,p,a,(c*N_fine):(N_fine+c*N_fine)] = np.concatenate((X_tmp[((N_fine/2)+1):(N_fine)],X_tmp[0:((N_fine/2)+1)]), axis=None)

# coeff_idx(a, p, b, f, Np, Nb)
# Generate coefficients - Antenna X Polarization X Beam X Coarse channel
coeff = np.ones(N_coarse*N_beam*N_pol*N_ant).reshape(N_coarse,N_beam,N_pol,N_ant)

#for b in range(0,N_beam):
#    coeff[2, 0, b, 2] = 1


# Beamform, convert to power and integrate over time samples
# coh_bf_idx(p, b, f, c, w, Np, Nb, Nc, Nf)
# pow_bf_idx(f, c, s, b, Nf, Nc, Ns)
x_pol = np.zeros(N_coarse*N_fine)
y_pol = np.zeros(N_coarse*N_fine)
bf_pow = np.zeros(N_win*N_coarse*N_fine).reshape(N_win, N_coarse*N_fine)
bf_sti = np.zeros(N_beam*N_coarse*N_fine).reshape(N_beam, N_coarse*N_fine)

for c in range(0,N_coarse):
    for b in range(0,N_beam):
        for w in range(0,N_win):
            x_pol[(c*N_fine):(N_fine+c*N_fine)] = np.dot(coeff[c,b,0,:],X[w,0,:,(c*N_fine):(N_fine+c*N_fine)])
            y_pol[(c*N_fine):(N_fine+c*N_fine)] = np.dot(coeff[c,b,1,:],X[w,1,:,(c*N_fine):(N_fine+c*N_fine)])
            bf_pow[w,(c*N_fine):(N_fine+c*N_fine)] = np.square(x_pol[(c*N_fine):(N_fine+c*N_fine)]) + np.square(y_pol[(c*N_fine):(N_fine+c*N_fine)])

        bf_sti[b,(c*N_fine):(N_fine+c*N_fine)] = np.sum(bf_pow[:,(c*N_fine):(N_fine+c*N_fine)], axis=0)

print("After for loops")

ant_idx = 0 # beam index to plot
pol_idx = 0 # polarization index to plot
time_idx = 0 # time sample index to plot
coarse_idx = 0 # coarse channel index to plot
# Plot intensity map of frequency vs. time
# "interpolation ='none'" removes interpolation which was there by default. 
# I'm only removing it for the sake of accurate analysis and diagnosis.
#plt.imshow(contents_array[0:N_win,0:N_coarse,beam_idx], extent=[1, N_coarse, 1, N_win], aspect='auto', interpolation='bicubic')
#plt.imshow(contents_array[coarse_idx,0:N_win,pol_idx,ant_idx,0:N_fine,iq_idx], extent=[1, N_fine, 1, N_win], aspect='auto', interpolation='none')
plt.imshow(X[0:N_win,pol_idx,ant_idx,0:N_coarse*N_fine], extent=[1, N_coarse*N_fine, 1, N_win], aspect='auto', interpolation='none')
plt.title('FFT Waterfall plot, antenna 0 (Frequency vs. time)')
plt.ylabel('Time samples')
plt.xlabel('Fine frequency channels')
plt.show()

# Plot of power spectrum
#plt.plot(contents_array[coarse_idx,time_idx,pol_idx,ant_idx,0:N_fine,iq_idx])
plt.plot(X[time_idx,pol_idx,ant_idx,0:N_coarse*N_fine])
plt.title('FFT at a time window')
plt.xlabel('Fine frequency channels')
plt.ylabel('Raw voltage (arb.)')
plt.show()

fig, axs = plt.subplots(2, 2)
fig.suptitle('FFT at different antennas')
axs[0, 0].plot(X[time_idx,pol_idx,0,0:N_coarse*N_fine])
axs[0, 0].set_title('Ant 1')
axs[0, 1].plot(X[time_idx,pol_idx,1,0:N_coarse*N_fine], 'tab:orange')
axs[0, 1].set_title('Ant 2')
axs[1, 0].plot(X[time_idx,pol_idx,2,0:N_coarse*N_fine], 'tab:green')
axs[1, 0].set_title('Ant 3')
axs[1, 1].plot(X[time_idx,pol_idx,57,0:N_coarse*N_fine], 'tab:red')
axs[1, 1].set_title('Ant 57')

# set the spacing between subplots
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.85,
                    wspace=0.4, 
                    hspace=0.6)

for ax in axs.flat:
    ax.set(xlabel='Frequency bins', ylabel='Raw voltage')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()
plt.show()

plt.imshow(bf_sti[:,:], extent=[1, N_coarse*N_fine, 1, N_beam], aspect='auto', interpolation='none')
plt.title('Power spectrum waterfall plot, antenna 0 (Frequency vs. time)')
plt.ylabel('Beam index')
plt.xlabel('Fine frequency channels')
plt.show()

# Plot of power spectrum
#plt.plot(contents_array[coarse_idx,time_idx,pol_idx,ant_idx,0:N_fine,iq_idx])
plt.plot(bf_sti[0,0:N_coarse*N_fine])
plt.title('Power spectum of a beam')
plt.xlabel('Fine frequency channels')
plt.ylabel('Power (arb.)')
plt.show()

fig, axs = plt.subplots(2, 2)
fig.suptitle('Power spectrum at different beams')
axs[0, 0].plot(bf_sti[0,0:N_coarse*N_fine])
axs[0, 0].set_title('Beam 1')
axs[0, 1].plot(bf_sti[1,0:N_coarse*N_fine], 'tab:orange')
axs[0, 1].set_title('Beam 2')
axs[1, 0].plot(bf_sti[2,0:N_coarse*N_fine], 'tab:green')
axs[1, 0].set_title('Beam 3')
axs[1, 1].plot(bf_sti[3,0:N_coarse*N_fine], 'tab:red')
axs[1, 1].set_title('Beam 4')

# set the spacing between subplots
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.85,
                    wspace=0.4, 
                    hspace=0.6)

for ax in axs.flat:
    ax.set(xlabel='Frequency bins', ylabel='Power (arb.)')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()
plt.show()
