## Patch version 0.1.1

This release fixes the installation ERROR on
r-devel-linux-x86_64-debian-clang, reported by Kurt Hornik on 2026-09-08. It
is the only change since 0.1.0.

RcppThread's thread pool uses an over-aligned `std::atomic`. clang emits a
call to `__atomic_compare_exchange` for it and gcc inlines it, so only that
one flavor failed.

The fix is a `configure` script that links against libatomic wherever it is
available. Where nothing calls into the library, `--as-needed` drops the flag
and the installed object is unchanged.

I reproduced the reported failure locally with clang 18.1.3 against 0.1.0, and
confirmed it is gone in 0.1.1. Under gcc 13.3.0 the linker drops the flag, and
neither version's `.so` depends on libatomic.

## Test environments

* local: Ubuntu 24.04, R 4.6.1, gcc 13.3.0
* local: Ubuntu 24.04, R 4.6.1, clang 18.1.3 (the configuration that fails
    for 0.1.0; verified against both 0.1.0 and 0.1.1)
* local: macOS Sequoia, R 4.6.0, Apple clang
* win-builder: R-devel, R-release, and R-oldrel

## R CMD check results

0 errors | 0 warnings | 1 note

The note appears only on the Ubuntu machines. macOS and all three win-builder
runs return OK.

* `checking compilation flags used ... NOTE: non-portable flag(s):
    '-mno-omit-leaf-frame-pointer'`. This comes from the Ubuntu R build's
    default `CXXFLAGS`, not from the package; `src/Makevars.in` sets only
    `$(SHLIB_OPENMP_CXXFLAGS)` and `-DARMA_64BIT_WORD=1`. It was present for
    0.1.0 as well and did not appear on any CRAN flavor.

`configure` and `cleanup` pass `checkbashisms`.
