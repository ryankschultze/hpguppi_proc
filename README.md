# hpguppi_daq

To configure and compile the first time:
```
 $ cd src
 $ libtoolize
 $ autoconf
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/lib --with-libcoherent_beamformer=../lib
 $ make
 ```

To configure and compile any other time after:
```
 $ cd src
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/lib --with-libcoherent_beamformer=../lib
 $ make
 ```

If everything is configured as it should be, you can just:
```
 $ make
 ```

NOTE: libsla may be in a different directory depending on the machine.

If python library is required:

```
 $ cd src
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/listen/lib --with-libcoherent_beamformer=../lib --with-libpython3.7m=/opt/conda/lib
 $ make
 ```
