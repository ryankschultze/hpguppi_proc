# This script reshapes the one dimensional array from the beamformer and plots the response for analysis
import matplotlib.pyplot as plt
from array import array
import numpy as np

# Open binary file containing beamformer output
filename = "/datag/users/mruzinda/o/output_d_fft_bf.bin"

# Read file contents: np.fromfile(filename, dtype=float, count=- 1, sep='', offset=0)
contents_float = np.fromfile(filename, dtype=np.float32)

print(len(contents_float))
print(contents_float[0])

# Array dimensions
N_beam = 61 # 64
N_time = 1 # STI windows
# 1k mode
#N_coarse = 1
#N_fine = (4096*1024)/8 # N_coarse*2^19
# 4k mode
N_coarse = 4 # 4
N_fine = (1024*1024)/8 # N_coarse*2^17 
# 32k mode
#N_coarse = 32
#N_fine = (128*1024)/8 # N_coarse*2^14

#define pow_bf_idx(f, c, s, b, Nf, Nc, Ns)
# Reshape array to 3D -> Fine channel X Coarse channel X Time samples X Beams
contents_array = contents_float[0:(N_time*N_coarse*N_fine*N_beam)].reshape(N_beam,N_time,N_coarse*N_fine)

beam_idx = 2 # beam index to plot
time_idx = 0 # time sample index to plot

if N_time > 1:
    # Plot intensity map of frequency vs. time
    # "interpolation ='none'" removes interpolation which was there by default. 
    # I'm only removing it for the sake of accurate analysis and diagnosis.
    #plt.imshow(contents_array[0:N_time,0:N_fine,beam_idx], extent=[1, N_fine, 1, N_time], aspect='auto', interpolation='bicubic')
    plt.imshow(contents_array[beam_idx,0:N_time,0:(N_coarse*N_fine)], extent=[1, (N_coarse*N_fine), 1, N_time], aspect='auto', interpolation='none')
    plt.title('Waterfall (Frequency vs. time)')
    plt.ylabel('Time samples')
    plt.xlabel('Frequency bins')
    plt.show()

#print(contents_array[0,0,(N_fine-10):(N_fine-1)])
#print(contents_array[0:N_beam,0,5])

# Plot of power spectrum
plt.plot(contents_array[beam_idx,time_idx,0:(N_coarse*N_fine)])
plt.title('Power spectrum')
plt.xlabel('Frequency bins')
plt.ylabel('Power (arb.)')
plt.show()

fig, axs = plt.subplots(2, 2)
fig.suptitle('Power spectra of individual beams')
axs[0, 0].plot(contents_array[0,time_idx,0:(N_coarse*N_fine)])
axs[0, 0].set_title('Beam 1')
axs[0, 1].plot(contents_array[1,time_idx,0:(N_coarse*N_fine)], 'tab:orange')
axs[0, 1].set_title('Beam 2')
axs[1, 0].plot(contents_array[2,time_idx,0:(N_coarse*N_fine)], 'tab:green')
axs[1, 0].set_title('Beam 3')
axs[1, 1].plot(contents_array[3,time_idx,0:(N_coarse*N_fine)], 'tab:red')
axs[1, 1].set_title('Beam 4')

# set the spacing between subplots
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.85,
                    wspace=0.4, 
                    hspace=0.6)

for ax in axs.flat:
    ax.set(xlabel='Frequency bins', ylabel='Power')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()
plt.show()

if N_time > 1:
    # Plot of power over time
    freq_idx = 5 # Frequency to plot
    plt.plot(contents_array[beam_idx,0:N_time,freq_idx])
    plt.title('Power over time at a particular frequency')
    plt.xlabel('Time samples')
    plt.ylabel('Power (arb.)')
    plt.show()

    fig, axs = plt.subplots(2, 2)
    fig.suptitle('Power over time of individual beams')
    axs[0, 0].plot(contents_array[0,0:N_time,freq_idx])
    axs[0, 0].set_title('Beam 1')
    axs[0, 1].plot(contents_array[1,0:N_time,freq_idx], 'tab:orange')
    axs[0, 1].set_title('Beam 2')
    axs[1, 0].plot(contents_array[2,0:N_time,freq_idx], 'tab:green')
    axs[1, 0].set_title('Beam 3')
    axs[1, 1].plot(contents_array[33,0:N_time,freq_idx], 'tab:red')
    axs[1, 1].set_title('Beam 33')

    # set the spacing between subplots
    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.85, 
                        wspace=0.4, 
                        hspace=0.6)

    for ax in axs.flat:
        ax.set(xlabel='Time samples', ylabel='Power')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    #for ax in axs.flat:
    #    ax.label_outer()
    plt.show()

# Check with incrementing set of simulated data and coefficients
chk_flag = 1
if(chk_flag==1):
    beam_idx = 64 # Change beam index to see what the output of the corresponding beam should be
    tmp_calc = 0
    for i in range(1,62):
        #tmp_calc = tmp_calc + beam_idx*i
        tmp_calc = tmp_calc + i
        
    print(tmp_calc*64*8) # Output of beamformer
    #print(tmp_calc*tmp_calc + tmp_calc*tmp_calc) # Output power with imaginary part set to zero
