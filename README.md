# hpguppi_daq

To configure and compile the first time:
```
 $ cd src
 $ libtoolize
 $ autoconf
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/lib --with-libcoherent_beamformer=/home/mruzinda/hpguppi_proc/lib
 $ make
 ```

To configure and compile any other time after:
```
 $ cd src
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/lib --with-libcoherent_beamformer=/home/mruzinda/hpguppi_proc/lib
 $ make
 ```

If everything is configured as it should be, you can just:
```
 $ make
 ```

NOTE: libsla may be in a different directory depending on the machine.
NOTE: For the coherent beamformer library replace "/home/mruzinda" with the location of "hpguppi_proc/lib"

If python library is required:

```
 $ cd src
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/listen/lib --with-libcoherent_beamformer=/home/mruzinda/hpguppi_proc/lib --with-libpython3.7m=/opt/conda/lib
 $ make
 ```
