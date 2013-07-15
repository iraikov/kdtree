kdtree
======

Standard ML implementation of the K-d tree spatial indexing data structure.
 
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


KDTREE         -Signature-


 type point
 type kdtree

 val empty: kdtree -> bool

Returns true if the given tree is empty, false otherwise.

 val app: (point -> unit) -> kdtree -> unit

Applies the given procedure to every point in the tree.

 val appi: (int * point -> unit)  -> kdtree -> unit

Applies the given procedure to every point in the tree and its index.

 val foldl: ((point * 'a) -> 'a) -> 'a -> kdtree -> 'a
 val foldli: ((int * point * 'a) -> 'a) -> 'a -> kdtree -> 'a
 val foldr: ((point * 'a) -> 'a) -> 'a -> kdtree -> 'a
 val foldri: ((int * point * 'a) -> 'a) -> 'a -> kdtree -> 'a

Folds over the tree.

 val ifoldr: ((int * 'a) -> 'a) -> 'a -> kdtree -> 'a

Folds over the tree using only point indices.

 val size: kdtree -> int

Returns the size of the tree.

 val toList: kdtree -> point list

Returns the list of points contained in the tree.

 val toIndexList: kdtree -> int list

Returns the list of points indices contained in the tree.

 val fromTensor: RTensor.tensor -> kdtree

Given a tensor of dimensionality [N,K], constructs a KD-tree.

 val nearestNeighbor: kdtree -> real list -> int option

 val nearNeighbors: kdtree -> real -> real list -> int list

 val remove: kdtree -> real -> real list -> kdtree

 val kNearestNeighbors: kdtree -> int -> real list -> int list


KDtreeFn         -Functor-


functor KDTreeFn (val K : int
                  val distanceSquared : (real list) * (real list) -> real): KDTREE 

