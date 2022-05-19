from astropy import units as u
import setigen as stg

n_ant = 8
n_chan_perant = 64
n_obschan = n_ant*n_chan_perant

antenna = stg.voltage.Antenna(sample_rate=3e9*u.Hz,
                              fch1=6000e6*u.Hz,
                              ascending=True,
                              num_pols=2)

antenna.x.add_noise(v_mean=0,
                    v_std=1)

antenna.x.add_constant_signal(f_start=6002.2e6*u.Hz,
                              drift_rate=-2*u.Hz/u.s,
                              level=0.002)

digitizer = stg.voltage.RealQuantizer(target_fwhm=32,
                                      num_bits=8)

filterbank = stg.voltage.PolyphaseFilterbank(num_taps=8,
                                             num_branches=1024)

requantizer = stg.voltage.ComplexQuantizer(target_fwhm=32,
                                           num_bits=8)

block_size = 128*(2**20)
file_size = 16*(2**30)
blocks_per_file = file_size // block_size
rvb = stg.voltage.RawVoltageBackend(antenna,
                                    digitizer=digitizer,
                                    filterbank=filterbank,
                                    requantizer=requantizer,
                                    start_chan=0,
                                    num_chans=n_obschan,
                                    block_size=block_size,
                                    blocks_per_file=blocks_per_file,
                                    num_subblocks=32)

rvb.record(
	output_file_stem='../golden_input',
	num_blocks=blocks_per_file//16,
	length_mode='num_blocks',
	header_dict={
		'TELESCOP': 'SETIGEN'
	},
	verbose=True
)