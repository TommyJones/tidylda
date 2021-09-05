Here's what you need to do
1. Seed your RNG from R (still do this b/c C++ is faster)
2. Sample all your uniform or exponential RVs at once
  - consider doing this in parallel over how many threads you have
  - Only have to do the long jump thing once
3. Allocate those numbers in a structure that is cache access efficient
  - possibly best approach here is to have a nested vector of vectors of doubles of that conforms to [iterations][threads][doubles]
4. Rip through them as you sample


Another, possibly conceptually simpler thing to do
1. Seed the RNG once from R
2. For each iteration, sample `threads` number of integers
3. Pass those along to each thread to seed their RNG streams