from astropy import units as u
import setigen as stg

n_pols = 2
n_bits = 8
n_ant = 64
n_chan_perant = 128
n_time = 1024


file_bytesize = 1*(2**30)
block_bytesize = (n_time*n_ant*n_chan_perant*n_pols*n_bits)//8

n_obschan = n_ant*n_chan_perant
blocks_per_file = file_bytesize // block_bytesize

antenna = stg.voltage.Antenna(sample_rate=3e9*u.Hz,
                              fch1=6000e6*u.Hz,
                              ascending=True,
                              num_pols=n_pols)

antenna.x.add_noise(v_mean=0,
                    v_std=1)

antenna.x.add_constant_signal(f_start=6002.2e6*u.Hz,
                              drift_rate=-2*u.Hz/u.s,
                              level=0.002)

digitizer = stg.voltage.RealQuantizer(target_fwhm=32,
                                      num_bits=8)

filterbank = stg.voltage.PolyphaseFilterbank(num_taps=8,
                                             num_branches=n_obschan*2)

requantizer = stg.voltage.ComplexQuantizer(target_fwhm=32,
                                           num_bits=n_bits)

rvb = stg.voltage.RawVoltageBackend(antenna,
                                    digitizer=digitizer,
                                    filterbank=filterbank,
                                    requantizer=requantizer,
                                    start_chan=0,
                                    num_chans=n_obschan,
                                    block_size=block_bytesize,
                                    blocks_per_file=blocks_per_file,
                                    num_subblocks=32)

rvb.record(
	output_file_stem='../golden_synthesized_input',
	num_blocks=(2**30)//block_bytesize,
	length_mode='num_blocks',
	header_dict={
		'TELESCOP': 'SETIGEN',
    'NANTS'   : str(n_ant),
    'OBSID'   : 'SYNTHETIC',
    'piperblk': n_time
	},
	verbose=True
)