# Research

This code is all of the code used in my Undergraduate Honors Thesis research, ranging from prototypes and proofs of concept to the final k-d tree and octree.

The code is finished enough that I was able to get the results I needed, but still needs some work before it could be used as a full fluid simulation tool.
This involves finishing advection on the k-d tree, constraining facial velocity on T-junctions, and implementing level set advection.

It also needs way better programming practices. Currently I have no header files, no object oriented stuff for the things the k-d and quadtree have in common, and no separate test files or framework.
