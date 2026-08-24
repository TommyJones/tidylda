// warp_rng.h --- random number generation for the warpLDA engine
//
// Implements roadmap D12 and D13. Two rules govern everything here:
//
//   1. R's RNG may only be touched from the main thread. Worker threads get
//      their own generators, seeded from a master value drawn from R's stream
//      at initialization so that set.seed() still governs the whole run.
//
//   2. Streams are derived from the WORK ITEM, not from the thread:
//
//          seed = f(master, iteration, pass, index)
//
//      where index is the document id in the doc pass and the word id in the
//      word pass. Document d therefore consumes the same random numbers no
//      matter which thread processes it, so results are reproducible
//      independent of thread count -- not merely at a fixed thread count.
//
// This is built in Phase 2, while the engine is still single-threaded, so that
// Phase 5 changes only how work is scheduled and never which stream a work item
// sees. "Phase 5 at threads = 1 reproduces Phase 2 bit for bit" is then a real
// regression check that separates a threading bug from an RNG-change bug.
//
// WHY THE SEEDS ARE HASHED (D13). xorshift-family generators produce CORRELATED
// streams from nearby seeds, so seeding thousands of per-document streams with
// consecutive integers would quietly bias the sampler in a way that looks like
// nothing at all until the benchmarks come back subtly wrong. splitmix64 is the
// standard companion for exactly this: it avalanches a counter into a
// well-separated 64-bit state.

#pragma once

#include <cstdint>

namespace warp {

// splitmix64. Public domain, Sebastiano Vigna. Used to expand a counter into
// generator state -- one call per output, no correlation between nearby inputs.
inline uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}


// xoshiro256++ 1.0. Public domain, Blackman and Vigna. Fast, small state,
// passes BigCrush. Seeded through splitmix64 as its authors require.
class Xoshiro256pp {
public:
  explicit Xoshiro256pp(uint64_t seed = 0) { reseed(seed); }

  void reseed(uint64_t seed) {
    // Each state word gets its own splitmix64 output so that a zero seed, or
    // two seeds differing in one bit, still give well-separated states.
    for (int i = 0; i < 4; i++) {
      seed = splitmix64(seed);
      s_[i] = seed;
    }
  }

  uint64_t next() {
    const uint64_t result = rotl(s_[0] + s_[3], 23) + s_[0];
    const uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

  // Uniform on [0, 1). 53 significant bits, the most a double can hold.
  double unif() {
    return static_cast<double>(next() >> 11) * 0x1.0p-53;
  }

  // Uniform on {0, ..., n-1}, without the modulo bias that `next() % n` has.
  //
  // Plain rejection sampling: draw until the value falls in the largest
  // multiple of n below 2^64, then reduce. Unbiased, and portable to anything
  // with 64-bit integers.
  //
  // Lemire's multiply-shift method is faster, but it needs a 128-bit integer
  // type, which is a compiler extension rather than standard C++. Guarding it
  // with #if would make the number of draws consumed depend on the compiler,
  // and therefore make a stream depend on the platform -- a bad trade for a
  // division. At N = 1.8e5 tokens and 200 iterations this costs well under a
  // second against fits that run for minutes.
  uint64_t below(uint64_t n) {
    if (n <= 1) return 0;
    const uint64_t threshold = (0 - n) % n;  // == 2^64 mod n
    uint64_t r;
    do {
      r = next();
    } while (r < threshold);
    return r % n;
  }

private:
  static uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  }
  uint64_t s_[4];
};


// Which pass a work item belongs to. Part of the seed so that document 7 in the
// doc pass and word 7 in the word pass never share a stream.
//
// `init` exists for the same reason and is not decorative: initialization runs
// once, conceptually at iteration 0, and the doc pass at iteration 0 already
// claims work_item_rng(master, 0, Pass::doc, d). Without a separate value the
// two would share a stream, so a token's starting topic and its very first
// proposal would be driven by the same uniform -- a correlation with no reason
// to exist and no obvious symptom.
enum class Pass : uint64_t { doc = 0, word = 1, init = 2 };


// The work-item seed of D12.
//
// The four inputs are mixed by hashing them in sequence rather than by packing
// them into disjoint bit ranges. Packing would make the seeds of adjacent
// documents differ in their low bits only -- exactly the nearby-seed case D13
// warns about -- and would also silently truncate on a corpus large enough to
// overflow whatever field width was chosen.
inline uint64_t work_item_seed(uint64_t master,
                               uint64_t iteration,
                               Pass pass,
                               uint64_t index) {
  uint64_t h = splitmix64(master);
  h = splitmix64(h ^ splitmix64(iteration));
  h = splitmix64(h ^ splitmix64(static_cast<uint64_t>(pass)));
  h = splitmix64(h ^ splitmix64(index));
  return h;
}


// Convenience: a generator already seeded for one work item.
inline Xoshiro256pp work_item_rng(uint64_t master,
                                  uint64_t iteration,
                                  Pass pass,
                                  uint64_t index) {
  return Xoshiro256pp(work_item_seed(master, iteration, pass, index));
}

}  // namespace warp
