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


functor KDTreeFn (val N : int
                  val distanceSquared : (real list) * (real list) -> real) = 
struct

structure IntArraySort = ArrayMergeSortFn (IntArray)

datatype kdtree' =
         KdNode of { left: kdtree', i: int, right: kdtree', axis: int }
       | KdLeaf of { ii: IntVector.vector, axis: int }

type 'a kdtree = { P: RTensor.tensor, T: kdtree' }

exception Point

fun point P i = RTensorSlice.fromto ([i,0],[i,N-1],P)

fun pointList p = List.rev (RTensorSlice.foldl (op ::) [] p)

fun coord point = 
    let val base  = RTensorSlice.base point 
        val shape = RTensorSlice.shape point 
        val lo    = Range.first (RTensorSlice.range point)
        val hi    = Range.last (RTensorSlice.range point)
    in
        case (shape,lo,hi) of
            ([1,n],[p,0],[p',n']) => 
            (if ((p=p') andalso (n=n'+1))
             then (fn (i) => RTensor.sub(base,[p,i]))
             else raise Point)
          | _ => raise Point
    end


fun pointCoord P (i,c) = RTensor.sub (P,[i,c])


fun compareDistance reltol probe (a,b) = 
    let
(*
        val _ = (print "compareDistance: probe = " ; TensorFile.realListLineWrite TextIO.stdOut probe)
        val _ = (print "compareDistance: a = " ; TensorFile.realListLineWrite TextIO.stdOut (pointList a))
        val _ = (print "compareDistance: b = " ; TensorFile.realListLineWrite TextIO.stdOut (pointList b))

        val _ = (print ("compareDistance: dist^2(probe,a) = " ^ (Real.toString (distanceSquared (probe, pointList a))) ^ "\n"))
        val _ = (print ("compareDistance: dist^2(probe,b) = " ^ (Real.toString (distanceSquared (probe, pointList b))) ^ "\n"))
*)
        val delta  = Real.- (distanceSquared (probe, pointList a), distanceSquared (probe, pointList b))
(*
        val _ = print ("compareDistance: delta = " ^ (Real.toString delta) ^ "\n")
*)
                   
    in
        case reltol of
            NONE   => (if Real.< (delta, 0.0) then LESS else GREATER)
          | SOME r => (if Real.< (Real.* (delta,delta), r)
                       then EQUAL
                       else (if Real.< (delta, 0.0) 
                             then LESS else (if Real.> (delta, 0.0) then GREATER else EQUAL)))
    end



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



fun empty {P,T} =
    case T of
        KdNode _         => false
      | KdLeaf {ii,axis} => (IntVector.length ii)=0


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

fun size {P,T} = hd (RTensor.shape P)

fun toList t = foldr (op ::) [] t


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
            val x = pointCoord P (i, axis)
            val leftValid = (List.all (fn y => Real.< (coord y axis, x))
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


(* Constructs a kd-tree from a tensor of points, starting with the given depth. *)

fun fromTensorWithDepth P depth =
    let 
        val pointCoord'  = pointCoord P
        val [sz,_]     = RTensor.shape P
        val bucketSize = 10 * (Int.max (Real.ceil (Math.log10 (Real.fromInt sz)),  1))

        val sub        = Unsafe.IntArray.sub

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

                val axis   = Int.mod (depth, N)

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


        fun fromTensorWithDepth' (I,m,n,depth) =
            (let
                 val k = n - m
             in
                 if (k <= bucketSize) orelse (k <= 1)
                 then 
                     let
                         val ii  = IntArraySlice.vector (IntArraySlice.slice (I, m, SOME (k+1)))
                     in
                         KdLeaf {ii=ii, axis=Int.mod (depth, N)}
                     end
                 else 
                     (let 
                          val axis        = Int.mod (depth, N)
                          val depth'      = depth+1
                          val mmedian     = findMedian (I,m,n,depth)
                      in
                          case mmedian of
                              SOME (I',median) =>
                              (let
                                   val x = pointCoord' (sub (I',median), axis)
                                                                                                                                        
                                   val left  = fromTensorWithDepth' (I',m,median-1,depth')
                                   val right = fromTensorWithDepth' (I',median+1,n,depth')

                                           
                               in
                                   KdNode {left=left,i=sub(I',median),right=right,axis=axis}
                               end)
                            | NONE => (KdLeaf {ii=IntArraySlice.vector (IntArraySlice.slice (I, m, SOME (k+1))),  
                                               axis=Int.mod (depth, N)})
                      end)
             end)

    in
        if sz=0 
        then KdLeaf {ii=IntVector.fromList [], axis=Int.mod (depth, N)}
        else fromTensorWithDepth' (IntArray.tabulate (sz, fn i => i), 0, sz-1, depth)
    end
    
   fun fromTensor P = {P=P,T=(fromTensorWithDepth P 0)}

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
                   minimumBy (IntVector.fromList candidates'') compareDistance''
               end

           and nearestNeighbor' t =
               case t of

                   KdLeaf { ii, axis } => 
                   minimumBy ii compareDistance''
                   
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

end
