This folder contains the files required to execute a multidimensional CR spline in C. 
The code can be reused to evaluate splines in Fluent.

Anton Fadic
14-07-2017

# Motivation
Motivation: Chemical kinetics are stiff systems of ODE, sometimes with additional
constraints such as algebraic conditions. Solving these systems numerically can be slow,
of the order of tens of ms per call. When coupling the chemical kinetics into CFD to study
the effect of transport phenomena, incorporating these kinetics directly **is slow**, especially
if many cells require to do these calls. Another issue is that these kinetics **may** be unstable,
leading to unreliable solutions if performed on the fly during a CFD iteration.

# Proposed solution
Using precomputed data for a wide range of conditions and storing them in memory. One of the fastest
ways to do this is implementing them in C which sets a good baseline benchmark for these purposes.

# Findings
This project led to many interesting observations for my own education. Mainly I found that Lookup tables should be avoided because:
- Curse of dimensionality: Lookup tables are not practical for many dimensions (N>=4). This is because the memory requirements tend to explode, which make lead to inability to allocate RAM, and also to slowness.
- Larger memory requirements are detrimental for speeding up the simulation, even if search is O(1) which it is in this case (as it is implemented efficiently). I am aware of other codes that are not O(1) in the search. The reason why this is O(1) is because I don't do a search for interpolating data. I use a functional approach that points directly to the location in memory of what I need.
- The reasons for their limitations are:
    - Memory access: typically the bottleneck is retrieving data from memory. This was confirmed with this project. This means that more memory calls reduce performance, and memory calls increase exponentially with the number of dimensions.
    - Efficiency degradation at larger memory requirements: If multiple memory addresses are used to store data due to the OS being unable to allocate **contiguous memory**, then this fragmentation produces significant performance losses.

# Recommendations
Due to these important caveats, is that I decided during my PhD that surrogate models should go in a different direction. The way that I explored is using Neural Networks, which you can find in my other repository. I showed that NN are superior to Lookup Tables. 