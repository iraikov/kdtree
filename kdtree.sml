(*
 An implementation of the K-d tree spatial indexing data structure.
 
 http://en.wikipedia.org/wiki/K-d_tree

 The k-d tree is a binary search tree in which every branching node
 contains a k-dimensional point, and every leaf node contains a set
 of points. Every branching node represents a splitting hyperplane
 that divides the space into two parts, known as half-spaces.

 Points to the left of the splitting hyperplane are contained in the
 left subtree of the node and points right of the hyperplane are
 contained in the right subtree. The splitting hyperplane is chosen
 so as to be perpendicular to one of the axes in the k-dimensional
 space. The axis at each branching level is chosen in a round-robin
 fashion. For instance, in 3-D space, at level 0, the chosen axis is
 X, so points are divided according to their X-coordinates; at level
 1, the chosen axis is Y, so the points are divided according to
 their Y-coordinates; at the next branch level the chosen axis is Z,
 and so on.


 This code is based on the Haskell kd-tree library implementation of
 K-D trees.

 Copyright 2012-2013 Ivan Raikov and the Okinawa Institute of
 Science and Technology.

 This program is free software: you can redistribute it and/or
 modify it under the terms of the GNU General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 A full copy of the GPL license can be found at
 <http://www.gnu.org/licenses/>.
*)


structure KDTree = 
struct

structure IntArraySort = ArrayMergeSortFn (IntArray)

datatype kdtree' =
         KdNode of { left: kdtree', i: int, right: kdtree', axis: int }
       | KdLeaf of { ii: IntVector.vector, axis: int }

type 'a kdtree = { N: int, P: RTensor.tensor, T: kdtree' }

exception Point

fun onepoint (N,P) i = RTensorSlice.fromto ([0,i],[N-1,i],P)

fun coord point = 
    let val base  = RTensorSlice.base point 
        val shape = RTensorSlice.shape point 
        val lo    = Range.first (RTensorSlice.range point)
    in
        case (shape,lo) of
            ([n,1],[0,p]) => (fn (i) => RTensor.sub(base,[i,p]))
          | _ => raise Point
    end


fun empty {N,P,T} =
    case T of
        KdNode _         => false
      | KdLeaf {ii,axis} => (IntVector.length ii)=0


fun app f {N,P,T} =
    let
        val p = onepoint (N,P)
        fun app' T =
            case T of 
                KdNode {left,i,right,axis} => (app' left; f (p i); app' right)
              | KdLeaf {ii,axis}           => IntVector.app (fn i => f (p i)) ii
    in
        app' T
    end


fun appi f {N,P,T} =
    let
        val p = onepoint (N,P)
        fun appi' T =
            case T of 
                KdNode {left,i,right,axis} => (appi' left; f (i, p i); appi' right)
              | KdLeaf {ii,axis}           => IntVector.app (fn i => f (i, p i)) ii
    in
        appi' T
    end


fun foldl f init {N,P,T} =
    let
        val p = onepoint (N,P)
        fun foldl' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init' = foldl' init left
                                                   val init''= f (p i, init')
                                               in foldl' init'' right end)
              | KdLeaf {ii,axis}           => IntVector.foldl (fn (i,ax) => f (p i, ax)) init ii
    in
        foldl' init T
    end


fun foldli f init {N,P,T} =
    let
        val p = onepoint (N,P)
        fun foldli' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init' = foldli' init left
                                                   val init''= f (i, p i, init')
                                               in foldli' init'' right end)
              | KdLeaf {ii,axis}           => IntVector.foldl (fn (i,ax) => f (i, p i, ax)) init ii
    in
        foldli' init T
    end
                                          

fun foldr f init {N,P,T} =
    let
        val p = onepoint (N,P)
        fun foldr' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init' = foldr' init right
                                                   val init''= f (p i, init')
                                               in foldr' init'' left end)
              | KdLeaf {ii,axis}           => IntVector.foldr (fn (i,ax) => f (p i, ax)) init ii
    in
        foldr' init T
    end


fun foldri f init {N,P,T} =
    let
        val p = onepoint (N,P)
        fun foldri' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init'  = foldri' init right
                                                   val init'' = f (i, p i, init')
                                               in foldri' init'' left end)
              | KdLeaf {ii,axis}           => IntVector.foldr (fn (i,ax) => f (i, p i, ax)) init ii
    in
        foldri' init T
    end
                                        

fun toList t = foldr (op ::) [] t


fun subtrees {N,P,T} =
    let
        fun subtrees' T = 
            case T of 
                KdNode {left,i,right,axis} => (subtrees' left) @ [T] @ (subtrees' right)
              | KdLeaf {ii,axis}           => [T]
    in
        subtrees' T
    end


fun findFrom cmp (a,from) = IntArraySlice.find cmp (IntArraySlice.slice (a,from,NONE))


(* Constructs a kd-tree from a tensor of points, starting with the given depth. *)

fun fromTensorWithDepth (N,P) depth =
    let 
        val onepoint'  = onepoint (N,P)
        val [sz,_]     = RTensor.shape P
        val bucketSize = 10 * (Int.max (Real.ceil (Math.log10 (Real.fromInt sz)),  1))

        val sub        = Unsafe.IntArray.sub

        fun findMedian (I, m, n, depth) =
            let 
                val axis   = Int.mod (depth, N)
                val _      = IntArraySort.sortRange
                                 (fn (x,y) => let val px = onepoint' x
                                                  val py = onepoint' y
                                              in
                                                  Real.compare (coord px axis, coord py axis)
                                              end)
                                 (I,(m,n))

                val median   = m+(Int.quot (n-m,2))
                val mediani  = sub (I,median)
                val medianp  = onepoint' mediani
                val medianc  = coord medianp axis

                val median' = findFrom 
                                  (fn (x) => 
                                      let val px = onepoint' (sub (I, x))
                                      in
                                          Real.< (medianc, coord px axis)
                                      end)
                                  (I,median)
                              
	    in 
                case median' of
                    NONE => (I,m)
                  | SOME i => (I,i)
            end


        fun fromTensorWithDepth' (I,m,n,depth) =
            (let
                 val k = n - m
             in
                 if (k <= bucketSize) orelse (k <= 1)
                 then 
                     KdLeaf {ii=IntVector.map
                                    (fn i => sub(I,i))
                                    (IntArraySlice.vector (IntArraySlice.slice (I, m, SOME n))), 
                             axis=Int.mod (depth, N)}
                 else 
                     (let 
                          val depth'      = depth+1
                          val (I',median) = findMedian (I,m,n,depth)
                      in
                          KdNode {left=fromTensorWithDepth' (I',m,median,depth'),
                                  i=sub(I',median), 
                                  right=fromTensorWithDepth' (I',median+1,n,depth'),
                                  axis=Int.mod (depth, N)}
                      end)
             end)

    in
        if sz=0 
        then KdLeaf {ii=IntVector.fromList [], axis=Int.mod (depth, N)}
        else fromTensorWithDepth' (IntArray.tabulate (sz, fn i => i), 0, sz-1, depth)
    end

end
