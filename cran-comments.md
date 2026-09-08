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

The fix is a `configure` script that compiles and links that exact construct as
an executable, and adds `-latomic` only when the toolchain needs it. It is
plain POSIX `sh` rather than autoconf, it never fails the build, and it does
not add `-latomic` unconditionally, so platforms without libatomic (macOS) and
the Windows build, which uses `src/Makevars.win` and is unchanged, are
unaffected.

Nothing else has changed since 0.1.0. There are no changes to R code, to the
C++ sources, to documentation, or to results. The only other edit is seven
technical words added to `inst/WORDLIST` for the new NEWS entry.

## Test environments

* local: Ubuntu 24.04, R 4.6.1

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
