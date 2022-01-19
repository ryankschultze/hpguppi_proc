# hpguppi_daq

To compile:
```
 $ cd src
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/listen/lib --with-libcoherent_beamformer=/home/mruzinda/hpguppi_proc/lib
 $ make
 ```

NOTE: For the coherent beamformer library replace "/home/mruzinda" with the location of "hpguppi_proc/lib"

If python library is required:

```
 $ cd src
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/listen/lib --with-libcoherent_beamformer=/home/mruzinda/hpguppi_proc/lib --with-libpython3.7m=/opt/conda/lib
 $ make
 ```
