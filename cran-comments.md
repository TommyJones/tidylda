## Patch version 0.1.1

This release fixes the installation ERROR on
r-devel-linux-x86_64-debian-clang, reported by Kurt Hornik on 2026-09-08 with
a correction deadline of 2026-09-29. It is the only change in this version.

The package linked but failed to load on that flavor after its update to the
LLVM 23 packages from Debian unstable:

```
unable to load shared object '.../tidylda.so':
  .../tidylda.so: undefined symbol: __atomic_compare_exchange
```

RcppThread's thread pool keeps each worker's loop range in an over-aligned
`std::atomic` (`alignas(64)` over an eight-byte payload). clang treats an
over-aligned atomic as not lock-free and emits a call to the size-generic
`__atomic_compare_exchange`, which lives in libatomic; gcc inlines the same
construct, which is why the other twelve flavors were unaffected. `-shared`
does not diagnose undefined symbols, so the failure surfaced at `dyn.load()`
rather than at link time.

The fix is a `configure` script that links a trivial program with `-latomic`
and adds the flag to `PKG_LIBS` wherever that succeeds. It is plain POSIX `sh`
rather than autoconf and it never fails the build. Where libatomic does not
exist the link fails, the flag is omitted, and the flags are identical to
0.1.0's; that covers macOS. The Windows build uses `src/Makevars.win` and is
unchanged.

The flag is not conditioned on the compiler. The Linux linkers on your check
flavors default to `--as-needed`, which keeps a library only if something
actually calls into it, so on the gcc flavors `-latomic` is dropped and the
resulting `.so` has no dependency on libatomic. Probing for the flag rather
than for the failing construct itself keeps the script from having to
reproduce, and track, a construct that lives in a dependency's header.

I reproduced the failure locally against clang 18.1.3, before and after the
fix. Building 0.1.0 with clang gives exactly the symbol reported on your
machine:

```
unable to load shared object '.../tidylda.so':
  .../tidylda.so: undefined symbol: __atomic_compare_exchange
```

Building 0.1.1 with the same compiler, `configure` reports "libatomic is
available; adding -latomic", `-latomic` appears in the link line, `readelf -d`
on the installed `tidylda.so` shows a dependency on `libatomic.so.1`, the
package loads, and the test suite passes (395 passed, 0 failed).

Building with gcc 13.3.0, `-latomic` is passed as well, but `--as-needed`
drops it: `readelf -d` shows no libatomic dependency and `nm -D` reports no
undefined `__atomic_*` symbol, so the installed object is as it was for 0.1.0.
The test suite passes there too (395 passed, 0 failed).

Nothing else has changed since 0.1.0. There are no changes to R code, to the
C++ sources, to documentation, or to results. The only other edit is seven
technical words added to `inst/WORDLIST` for the new NEWS entry.

## Test environments

* local: Ubuntu 24.04, R 4.6.1, gcc 13.3.0
* local: Ubuntu 24.04, R 4.6.1, clang 18.1.3 (the configuration that fails
    for 0.1.0; verified against both 0.1.0 and 0.1.1)

## R CMD check results

0 errors | 1 warning | 1 note

Both are properties of this machine rather than of the package:

* `checking top-level files ... WARNING: A complete check needs the
    'checkbashisms' script.` The script is not installed here. `configure` and
    `cleanup` were checked against the real `checkbashisms` from devscripts
    separately, and it reports no bashisms in either.
* `checking compilation flags used ... NOTE: non-portable flag(s):
    '-mno-omit-leaf-frame-pointer'`. This comes from the Ubuntu R build's
    default `CXXFLAGS`, not from the package; `src/Makevars.in` sets only
    `$(SHLIB_OPENMP_CXXFLAGS)` and `-DARMA_64BIT_WORD=1`. It was present for
    0.1.0 as well and did not appear on any CRAN flavor.
