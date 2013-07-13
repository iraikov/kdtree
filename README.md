
An implementation of the K-d tree spatial indexing data structure.
 
 http://en.wikipedia.org/wiki/K-d_tree

The k-d tree is a binary search tree in which every branching node
contains a k-dimensional point, and every leaf node contains a set of
points. Every branching node represents a splitting hyperplane that
divides the space into two parts, known as half-spaces.

Points to the left of the splitting hyperplane are contained in the
left subtree of the node and points right of the hyperplane are
contained in the right subtree. The splitting hyperplane is chosen so
as to be perpendicular to one of the axes in the k-dimensional
space. The axis at each branching level is chosen in a round-robin
fashion. For instance, in 3-D space, at level 0, the chosen axis is X,
so points are divided according to their X-coordinates; at level 1,
the chosen axis is Y, so the points are divided according to their
Y-coordinates; at the next branch level the chosen axis is Z, and so
on.

