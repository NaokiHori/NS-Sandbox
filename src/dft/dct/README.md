# `dct`

Discrete cosine transforms (`DCT`) of types II and III based on the algorithm presented by Lee in 1984.
The transforms are computed for inputs of identical size (`nitems`) across multiple runs (`repeat_for`), with support for thread parallelism.

## Time Complexity

The functions achieve a time complexity of `O(N log N)` for power-of-two sizes.
If the size of the sub-problem is not a power of two, the algorithm falls back to a naive `O(N^2)` approach.

## References

- Lee, "A New Algorithm to Compute the Discrete Cosine Transform", *IEEE T. Acoust. Speech*, 1984

