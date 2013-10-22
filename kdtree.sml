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


structure IntVectorUtils =
struct

fun minimumBy v cmpfn: int option =
    let
        fun recur 0 m = m
          | recur i m = (let val i' = IntVector.sub (v,i-1)
                             val r  = cmpfn (i',m) 
                         in 
                          (case r of LESS => recur (i-1) i'
                                   | _    => recur (i-1) m )
                         end)
        val n = IntVector.length v
    in
        if n = 0 then NONE else SOME (recur (n-1) (IntVector.sub (v,n-1)))
    end

fun filter v ffn: int list =
    IntVector.foldl (fn (i,ax) => if (ffn i) then i::ax else ax) [] v

fun toList v = 
    IntVector.foldl (fn (i,ax) => i::ax) [] v


end


structure ListUtils =
struct

fun minimumBy lst cmpfn: int option =
    let
        fun recur [] m = m
          | recur (i::rest) m = 
            (let 
                 val r  = cmpfn (i,m) 
             in 
                 (case r of LESS => recur (rest) i
                          | _    => recur (rest) m )
             end)
    in
        case lst of
            [] => NONE 
          | _  => SOME (recur (tl lst) (hd lst))
    end


end

signature KDTREE = 
sig

type point
type pointStorage
type kdtree

val empty: kdtree -> bool

val app: (point -> unit) -> kdtree -> unit
val appi: (int * point -> unit)  -> kdtree -> unit


val foldl: ((point * 'a) -> 'a) -> 'a -> kdtree -> 'a
val foldli: ((int * point * 'a) -> 'a) -> 'a -> kdtree -> 'a

val foldr: ((point * 'a) -> 'a) -> 'a -> kdtree -> 'a
val foldri: ((int * point * 'a) -> 'a) -> 'a -> kdtree -> 'a

val ifoldr: ((int * 'a) -> 'a) -> 'a -> kdtree -> 'a

val size: kdtree -> int

val toList: kdtree -> point list

val toIndexList: kdtree -> int list

val isValid: kdtree -> bool

val allSubtreesAreValid: kdtree -> bool

val make: pointStorage -> kdtree

val nearestNeighbor: kdtree -> real list -> int option

val nearNeighbors: kdtree -> real -> real list -> int list

val remove: kdtree -> real -> real list -> kdtree

val kNearestNeighbors: kdtree -> int -> real list -> int list

end

functor KDTreeFn (
                  type point
                  type pointStorage
                  val point: pointStorage -> int -> point
                  val pointList: point -> real list
                  val coord: point -> int -> real
                  val pointCoord: pointStorage -> (int * int) -> real
                  val pointStorageSize: pointStorage -> int
                  val K : int
                  val distanceSquared : (real list) * (real list) -> real): KDTREE = 
struct

structure IntArraySort = ArrayMergeSortFn (IntArray)

type point = point
type pointStorage = pointStorage

datatype kdtree' =
         KdNode of { left: kdtree', i: int, right: kdtree', axis: int }
       | KdLeaf of { ii: IntVector.vector, axis: int }

type kdtree = { P: pointStorage, T: kdtree' }

exception Point
exception IndexArray


fun compareDistance reltol probe (a,b) = 
    let
        val delta  = Real.- (distanceSquared (probe, pointList a), distanceSquared (probe, pointList b))
    in
        case reltol of
            NONE   => (if Real.< (delta, 0.0) then LESS else GREATER)
          | SOME r => (if Real.< (Real.* (delta,delta), r)
                       then EQUAL
                       else (if Real.< (delta, 0.0) 
                             then LESS else (if Real.> (delta, 0.0) then GREATER else EQUAL)))
    end



fun empty' T =
    case T of
        KdNode _         => false
      | KdLeaf {ii,axis} => (IntVector.length ii)=0
    

fun empty {P,T} = empty' T


fun app f {P,T} =
    let
        val p = point P
        fun app' T =
            case T of 
                KdNode {left,i,right,axis} => (app' left; f (p i); app' right)
              | KdLeaf {ii,axis}           => IntVector.app (fn i => f (p i)) ii
    in
        app' T
    end


fun appi f {P,T} =
    let
        val p = point P
        fun appi' T =
            case T of 
                KdNode {left,i,right,axis} => (appi' left; f (i, p i); appi' right)
              | KdLeaf {ii,axis}           => IntVector.app (fn i => f (i, p i)) ii
    in
        appi' T
    end


fun foldl f init {P,T} =
    let
        val p = point P
        fun foldl' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init' = foldl' init left
                                                   val init''= f (p i, init')
                                               in foldl' init'' right end)
              | KdLeaf {ii,axis}           => IntVector.foldl (fn (i,ax) => f (p i, ax)) init ii
    in
        foldl' init T
    end


fun foldli f init {P,T} =
    let
        val p = point P
        fun foldli' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init' = foldli' init left
                                                   val init''= f (i, p i, init')
                                               in foldli' init'' right end)
              | KdLeaf {ii,axis}           => IntVector.foldl (fn (i,ax) => f (i, p i, ax)) init ii
    in
        foldli' init T
    end
                                          

fun foldr f init {P,T} =
    let
        val p = point P
        fun foldr' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init' = foldr' init right
                                                   val init''= f (p i, init')
                                               in foldr' init'' left end)
              | KdLeaf {ii,axis}           => IntVector.foldr (fn (i,ax) => f (p i, ax)) init ii
    in
        foldr' init T
    end


fun foldri f init {P,T} =
    let
        val p = point P
        fun foldri' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init'  = foldri' init right
                                                   val init'' = f (i, p i, init')
                                               in foldri' init'' left end)
              | KdLeaf {ii,axis}           => IntVector.foldr (fn (i,ax) => f (i, p i, ax)) init ii
    in
        foldri' init T
    end


fun ifoldr f init {P,T} =
    let
        fun foldr' init T =
            case T of 
                KdNode {left,i,right,axis} => (let val init' = foldr' init right
                                                   val init''= f (i, init')
                                               in foldr' init'' left end)
              | KdLeaf {ii,axis}           => IntVector.foldr (fn (i,ax) => f (i, ax)) init ii
    in
        foldr' init T
    end


fun size {P,T} = pointStorageSize P

fun toList t = foldr (op ::) [] t

fun toIndexList t = ifoldr (op ::) [] t

fun subtrees {P,T} =
    let
        fun subtrees' T = 
            case T of 
                KdNode {left,i,right,axis} => (subtrees' left) @ [T] @ (subtrees' right)
              | KdLeaf {ii,axis}           => [T]
    in
        subtrees' T
    end


fun isValid {P,T} =
    case T of
        KdLeaf _ => true
      | KdNode {left,i,right,axis} =>
        let
            val x          = pointCoord P (i, axis)
            val leftValid  = (List.all (fn y => Real.< (coord y axis, x))
                                       (toList {P=P,T=left}))
            val rightValid = (List.all (fn y => Real.>= (coord y axis, x))
                                       (toList {P=P,T=right}))
        in
            leftValid andalso rightValid
        end

(* Checks whether the K-D tree property holds for the given tree and its subtreees *)

fun allSubtreesAreValid (t as {P,T}) = 
    List.all (fn (t') => isValid {P=P,T=t'}) (subtrees t)



fun findiFromTo cmp (a,from,to) = IntArraySlice.findi cmp (IntArraySlice.slice (a,from,SOME (to-from+1)))


(* Constructs a kd-tree from a storage of points, starting with the given depth. 
   If I is given, then only use the point indices contained in it, otherwise use all points. *)

fun makeWithDepth P I depth =
    let 
        val pointCoord'  = pointCoord P
        val sz           = pointStorageSize P
        val bucketSize   = 10 * (Int.max (Real.ceil (Math.log10 (Real.fromInt sz)),  1))

        val sub          = Unsafe.IntArray.sub

        fun findMedian (I, m, n, depth) =
            let 
                fun findGreaterCoord (ci,cc,to,axis) =
                    findiFromTo 
                        (fn (i,x) => 
                            let 
                                val cx = pointCoord' (sub (I, x), axis)
                            in
                                Real.< (cc, cx)
                            end)
                        (I,ci,to)

                val axis   = Int.mod (depth, K)

                val _      = IntArraySort.sortRange
                                 (fn (x,y) => 
                                     let 
                                         val cx = pointCoord' (x, axis)
                                         val cy = pointCoord' (y, axis)
                                     in
                                         Real.compare (cx, cy)
                                     end)
                                 (I,(m,n+1))

                val median   = m+(Int.quot (n-m,2))
                val medianc  = pointCoord' (sub (I,median), axis)
                val median'  = findGreaterCoord (median,medianc,n,axis)
                              
	    in 
                case median' of
                    SOME (i,_) => SOME (I,median+i)
                  | NONE => NONE
                    
            end


        fun makeWithDepth' (I,m,n,depth) =
            (let
                 val k = n - m
             in
                 if (k <= bucketSize) orelse (k <= 1)
                 then 
                     let
                         val ii  = IntArraySlice.vector (IntArraySlice.slice (I, m, SOME (k+1)))
                     in
                         KdLeaf {ii=ii, axis=Int.mod (depth, K)}
                     end
                 else 
                     (let 
                          val axis        = Int.mod (depth, K)
                          val depth'      = depth+1
                          val mmedian     = findMedian (I,m,n,depth)
                      in
                          case mmedian of
                              SOME (I',median) =>
                              (let
                                   val x = pointCoord' (sub (I',median), axis)
                                                                                                                                        
                                   val left  = makeWithDepth' (I',m,median-1,depth')
                                   val right = makeWithDepth' (I',median+1,n,depth')

                                           
                               in
                                   KdNode {left=left,i=sub(I',median),right=right,axis=axis}
                               end)
                            | NONE => (KdLeaf {ii=IntArraySlice.vector (IntArraySlice.slice (I, m, SOME (k+1))),  
                                               axis=Int.mod (depth, K)})
                      end)
             end)

    in
        if sz=0 
        then KdLeaf {ii=IntVector.fromList [], axis=Int.mod (depth, K)}
        else (case I of NONE => makeWithDepth' (IntArray.tabulate (sz, fn i => i), 0, sz-1, depth)
                      | SOME I' => if (IntArray.length I') <= sz 
                                   then makeWithDepth' (I', 0, sz-1, depth)
                                   else raise IndexArray)
    end
    
   fun make P = {P=P,T=(makeWithDepth P NONE 0)}


   (* Returns the index of the nearest neighbor of p in tree t. *)

   fun nearestNeighbor {P,T} probe =

       let
           val point'      = point P
           val pointCoord' = pointCoord P
           val compareDistance' = compareDistance (SOME 1E~16) probe

           fun compareDistance'' (i,j) =
               compareDistance' (point' i, point' j)

           fun findNearest (t1,t2,p,xp,xprobe) =
               let
                   val candidates' = 
                       case nearestNeighbor' t1 of
                           SOME best => [best,p]
                         | NONE      => [p]

                   val sphereIntersectsPlane = 
                       let
                           val delta = Real.- (xprobe, xp)
                           val candidate = point' (hd candidates')
                       in
                           Real.< (Real.* (delta,delta), 
                                   distanceSquared (probe, pointList candidate))
                       end

                   val candidates'' = if sphereIntersectsPlane
                                      then (case nearestNeighbor' t2 of
                                                SOME nn => candidates' @ [nn]
                                              | NONE => candidates')
                                      else candidates'

               in
                   IntVectorUtils.minimumBy (IntVector.fromList candidates'') compareDistance''
               end

           and nearestNeighbor' t =
               case t of

                   KdLeaf { ii, axis } => 
                   IntVectorUtils.minimumBy ii compareDistance''
                   
                 | KdNode { left, i, right, axis } => 
                   let 
                       val xprobe = List.nth (probe,axis)
                       val xp     = pointCoord' (i,axis)
                   in
                       if Real.< (xprobe, xp)
                       then findNearest (left, right, i, xp, xprobe)
                       else findNearest (right, left, i, xp, xprobe)


                   end
       in
           nearestNeighbor' T
       end


   (* Returns all neighbors within distance r from p in tree t. *)
   fun nearNeighbors {P,T} radius probe =
       let
           val point'      = point P
           val pointCoord' = pointCoord P
           val r2          = Real.* (radius, radius)
           fun filterIndices ii = IntVectorUtils.filter 
                                      ii (fn (i) => (Real.<= (distanceSquared (probe, pointList (point' i)), r2)))

           fun nearNeighbors' t =
               case t of
                   
                   KdLeaf { ii, axis } => filterIndices ii
                                      
	         | KdNode { left, i, right, axis } =>
                   
	           (let 
                        val maybePivot = filterIndices (IntVector.fromList [i])
	            in
		        if (empty' left) andalso (empty' right)
                        then maybePivot
                        else 
		            (let 
                                 val xprobe = List.nth (probe,axis)
                                 val xp     = pointCoord' (i, axis)
                             in
			         if (Real.<= (xprobe, xp))
                                 then
			             (let 
                                          val nearest = maybePivot @ (nearNeighbors' left)
                                      in
				          (if Real.> (Real.+ (xprobe, (Real.abs radius)), xp)
                                           then (nearNeighbors' right) @ nearest
				           else nearest)
                                      end)
                                 else
			             (let 
                                          val nearest = maybePivot @ (nearNeighbors' right)
                                      in
				          if Real.< (Real.- (xprobe, (Real.abs radius)), xp)
                                          then (nearNeighbors' left) @ nearest
                                          else nearest
                                      end)
                             end)
                    end)
       in
           nearNeighbors' T
       end
  
  (* Removes the point p from t. *)
   fun remove  {P,T} tol pkill =

       let
           val point'      = point P
           val pointCoord' = pointCoord P
           val tol2        = Real.* (tol, tol)
           fun filterIndices ii = IntVectorUtils.filter 
                                      ii (fn (i) => (Real.> (distanceSquared (pkill, pointList (point' i)), tol2)))

	   fun remove' t =
               case t of
                   
                   KdLeaf { ii, axis } => 
                   KdLeaf { ii = IntVector.fromList (filterIndices ii), axis = axis }
                                      
	         | KdNode { left, i, right, axis } =>
                   if (Real.> (distanceSquared (pkill, pointList (point' i)), tol2))
                   then 
                       (let
                            val I = IntArray.fromList ((toIndexList {P=P,T=left}) @ (toIndexList {P=P,T=right}))
                        in
                            makeWithDepth P (SOME I) axis
                        end)
                   else
		       (if (Real.< (List.nth (pkill,axis), pointCoord' (i, axis)))
                        then 
                            KdNode { left = remove' left, i=i, right=right, axis=axis }
                        else 
                            KdNode { left = left, i=i, right=remove' right, axis=axis })
       in
           { P=P, T=remove' T }
       end
  

  (* Returns the k nearest points to p within tree. *)
  fun kNearestNeighbors {P,T} k probe =
      let
          val point'      = point P
                            
          val compareDistance' = 
              compareDistance (SOME 1E~16) probe

          fun compareDistance'' (i,j) =
              compareDistance' (point' i, point' j)
              
          fun kNearestNeighbors' t =
              (case t of
                  KdLeaf { ii, axis } => 
                  (let
                       fun recur (res, ii, k) =
                           if ((k <= 0) orelse (List.null ii))
                           then res
                           else (let 
                                     val nearest = valOf (ListUtils.minimumBy ii compareDistance'')
                                 in
                                     recur (nearest::res, List.filter (fn (i) => not (i = nearest)) ii, k-1)
                                 end)
                   in
                       recur ([], IntVectorUtils.toList ii, k)
                   end)
                   
	        | KdNode { left, i, right, axis } =>
                  if (k <= 0) 
                  then []
                  else (let
                            val nearest = nearestNeighbor {P=P,T=t} probe
                        in
                            case nearest of
                                NONE => []
                              | SOME n => 
                                let
                                    val {P=_, T=t'} = remove {P=P,T=t} 1E~16 (pointList (point' n))
                                in
                                    n :: (kNearestNeighbors {P=P,T=t'} (k - 1) probe)
                                end
                        end))
      in
          kNearestNeighbors' T
      end
                       

      fun inBounds probe bMin bMax =
          ((ListPair.all (fn(x,y) => Real.<= (x,y)) (probe,bMax)) andalso
           (ListPair.all (fn(x,y) => Real.>= (x,y))) (bMin,probe))
          
(*
      fun rangeSearch' {P,T} bMin bMax d
kdtRangeSearchRec Empty _ _ _ = []
kdtRangeSearchRec (Leaf lpos ldat) bMin bMax d =
  if (vecDimSelect lpos d) > (vecDimSelect bMin d) &&
     (vecDimSelect lpos d) < (vecDimSelect bMax d) &&
     (kdtInBounds lpos bMin bMax) then [(lpos,ldat)]
                                  else []
kdtRangeSearchRec (Node left npos ndata right) bMin bMax d =
  if (vecDimSelect npos d) < (vecDimSelect bMin d) then
    kdtRangeSearchRec right bMin bMax (d+1)
  else
    if (vecDimSelect npos d) > (vecDimSelect bMax d) then
      kdtRangeSearchRec left bMin bMax (d+1)
    else
      if (kdtInBounds npos bMin bMax) then
        (npos,ndata) : ((kdtRangeSearchRec right bMin bMax (d+1))++
                        (kdtRangeSearchRec left bMin bMax (d+1)))
      else
        (kdtRangeSearchRec right bMin bMax (d+1))++
        (kdtRangeSearchRec left bMin bMax (d+1))

kdtRangeSearch :: (KDTreeNode a) -> Vec3 -> Vec3 -> [(Vec3,a)]
kdtRangeSearch t bMin bMax =
  kdtRangeSearchRec t bMin bMax 0
*)

end
