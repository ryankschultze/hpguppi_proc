# serial 1 upchannelizer_beamformer.m4
AC_DEFUN([AX_CHECK_UBF],
[AC_PREREQ([2.65])dnl
AC_ARG_WITH([upchannelizer_beamformer],
            AC_HELP_STRING([--with-upchannelizer_beamformer=DIR],
                           [Location of UBF files (..)]),
            [UBFDIR="$withval"],
            [UBFDIR=..])

orig_LDFLAGS="${LDFLAGS}"
LDFLAGS="${orig_LDFLAGS} -L${UBFDIR}/lib"
AC_CHECK_LIB([upchannelizer_beamformer], [init_upchan_beamformer],
             # Found
             AC_SUBST(UBF_LIBDIR,${UBFDIR}/lib),
             # Not found there, check UBFDIR
             AS_UNSET(ac_cv_lib_upchannelizer_beamformer_run_upchannelizer_beamformer)
             LDFLAGS="${orig_LDFLAGS} -L${UBFDIR}"
             AC_CHECK_LIB([upchannelizer_beamformer], [init_upchan_beamformer],
                          # Found
                          AC_SUBST(UBF_LIBDIR,${UBFDIR}),
                          # Not found there, error
                          AC_MSG_ERROR([UBF library not found])))
LDFLAGS="${orig_LDFLAGS}"

AC_CHECK_FILE([${UBFDIR}/include/upchannelizer_beamformer.h],
              # Found
              AC_SUBST(UBF_INCDIR,${UBFDIR}/include),
              # Not found there, check UBFDIR
              AC_CHECK_FILE([${UBFDIR}/upchannelizer_beamformer.h],
                            # Found
                            AC_SUBST(UBF_INCDIR,${UBFDIR}),
                            # Not found there, error
                            AC_MSG_ERROR([upchannelizer_beamformer.h header file not found])))

])
