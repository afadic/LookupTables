Cubic multidimensional interpolation. This implementation is designed to work with any number of dimensions. It is implemented efficiently with minimal cache misses. The code can be evaluate splines as surrogate models for example in Fluent.

Anton Fadic

# Motivation
Motivation: Chemical kinetics are stiff systems of ODE, sometimes with additional
constraints such as algebraic conditions. Solving these systems numerically can be slow,
of the order of tens of ms per call. When coupling the chemical kinetics into CFD to study
the effect of transport phenomena, incorporating these kinetics directly **is slow**, especially
if many cells make these calls such as fine discretization of a surface of a catalyst. Another issue is that these kinetics **may** be unstable, leading to unreliable solutions if performed on the fly during a CFD iteration.

# Proposed solutions
Using precomputed data for a wide range of conditions and storing them in memory. One of the fastest
ways to do this is implementing them in C which sets a good baseline benchmark for these purposes. Cubic splines are a good choice due to their C1 continuity property, needed for CFD simulations. This is why linear interpolants are not used, as their derivatives are not continuous.

I compared two approaches:
- Catmull-Rom splines (CRInterpolation): Similar to Hermite splines, but with finite difference approximated derivatives at the ends. This approach is a tradeoff at build time compared to Hermite splines as it doesn't require the function derivatives at each node.
- Hermite splines: For a N-dimensional hermite spline, you need all the mixed partial derivatives. For N=3 you need 7, for N=4, 15 and for N=6 you need 63. For N dimensions you need 2^N-1 partial derivatives. This slows down significantly the building phase, but it is the most accurate approach.

# Findings
- The first and most interesting finding is that using Catmull-Rom interpolation for the ammonia oxidation case gave speedup factors of about 8x for a typical lookup table architecture. The important part here is that the comparison is done against the efficient implementation in Sundials IDA in C, compared against an efficient implementation of the lookup tables with minimal cache misses. As CFD applications are typically optimized themselves, optimized comparisons are needed to bring a realistic baseline and to set the right expectations.

However, this project led to many other interesting findings. Mainly I found that Lookup tables should be avoided for N>=4 because:
- Curse of dimensionality: Lookup tables are not practical for many dimensions (N>=4). This is because the memory requirements tend to explode, which may lead to inability to allocate RAM, and also to slowness.
- Larger memory requirements are detrimental for speeding up the simulation, even if search is O(1) which it is in this case (as it is implemented efficiently). The reason why this is O(1) is because I don't do a search for interpolating data.  I use a functional approach that points directly to the location in memory of what I need. 

- The reasons for their limitations are:
    - Memory access: typically the bottleneck is retrieving data from memory. This was confirmed with this project. This means that more memory calls reduce performance, and memory calls increase exponentially with the number of dimensions. **Modern CPUs can perform 100+ flops in the time it takes for one DRAM access**. 
    - Efficiency degradation at larger memory requirements: If multiple memory addresses are used to store data due to the OS being unable to allocate **contiguous memory**, then this fragmentation produces significant performance losses.
For N<4 lookup tables scale relatively fine.

# Recommendations
Due to these important caveats, is that I decided during my PhD that surrogate models should go in a different direction. The way that I explored is using Neural Networks, which you can find in my other repository. I showed that NN are superior to Lookup Tables. 

# Conclusions
I was able to point out early on the limitations of different lookup tables approaches, specifically pointing out the limitations at higher dimensions and the memory access slowness compared to CPU performance. I proposed a more efficient alternative to produce accurate results in an efficient manner.  
