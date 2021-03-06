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


 This code is based on the Haskell KdTree library by Isaac Trotts and
 on Haskell kd-tree code by Matthew Sottile.

 Copyright 2013-2016 Ivan Raikov.

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


fun putStrLn (str) = 
    (TextIO.output (TextIO.stdOut, str);
     TextIO.output (TextIO.stdOut, "\n"))

structure IntArraySort = ArrayMergeSortFn (IntArray)

structure ListUtils =
struct

fun filterMap f lst = 
    let
        fun recur ([],ax) = List.rev ax
          | recur (x::lst,ax) = 
            case f x of
                SOME item => recur (lst, item::ax)
              | NONE  => recur (lst, ax)
    in
        recur (lst, [])
    end

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

(* A data structure for storing and looking up K-dimensional points *)

signature KPOINT_SPACE =
sig
    type point
    type pointSpace

    val K : int
    val empty: pointSpace
    val point: pointSpace -> int -> point
    val pointList: point -> real list
    val coord: point -> int -> real
    val pointCoord: pointSpace -> (int * int) -> real
    val size: pointSpace -> int
    val insert: real list * pointSpace -> int * pointSpace
end


(* Point space based on tensors *)
functor TensorPointSpaceFn (val K : int): KPOINT_SPACE =
struct

    exception Point
    exception NotImplemented

    type point = RTensorSlice.slice

    type pointSpace = RTensor.tensor
                        
    val K = K

    val empty = RTensor.new ([1,K], 0.0)
            
    fun point P i = RTensorSlice.fromto ([i,0],[i,K-1],P)
                    
    fun pointList p = List.rev (RTensorSlice.foldl (op ::) [] p)
                      
    fun coord point = 
        let val base  = RTensorSlice.base point 
            val shape = hd (RTensorSlice.shapes point)
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
    fun size P = hd (RTensor.shape P)
                 
    fun insert (point, P) = raise NotImplemented
end


(* Point space based on maps *)
functor MapPointSpaceFn (val K : int): KPOINT_SPACE =
struct

    exception Point
    exception NotImplemented

    type point = real list

    structure H = RedBlackMapFn (struct
			          type ord_key = int
			          val compare = Int.compare
			          end)
					    
    type pointSpace = point H.map

    fun insert (m, k, a) = H.insert(m, k, a)
    fun find (m, k) = H.lookup(m, k)
                        
    val K = K
            
    fun point P i = find (P,i)
                    
    fun pointList p = p
                      
    fun coord point i = List.nth (point, i)
        
    fun pointCoord P (i,c) = coord (find (P,i)) c
    fun size P = H.numItems(P)
                 
    val empty = H.empty
    fun insert (point, P) = let val i = size P in (i, H.insert(P, i, point)) end
end


(* Extensible point space based on a splay tree of tensors *)
(*
functor SplayTensorPointSpaceFn (val K : int): KPOINT_SPACE =
struct

    open SplayTree

    exception Point

    type point = RTensorSlice.slice

    type pointSpaceNode = int ref * int * RTensor.tensor

    val nodeSize = Word.toInt (Word.<<(0wx1, 0wx10))
                        
    val K = K
            
    fun point (SplayObj {(sz,lb,te), right, left}) i = 
        if (i >= lb) 
        then (if i < (lb + sz) then RTensorSlice.fromto ([i-lb,0],[i-lb,K-1],te) else point right i) 
        else point left i
      | point (SplayNil) i = 
        raise Point
                                                  
    fun pointCoord (SplayObj {(sz,lb,te), right, left}) (i,c) = 
        if (i >= lb) 
        then (if i < (lb + sz) then RTensor.sub (te,[i-lb,c]) else pointCoord right i) 
        else pointCoord left i
      | pointCoord (SplayNil) (i,c) =
        raise Point

    fun pointList p = List.rev (RTensorSlice.foldl (op ::) [] p)
                      
    fun coord point = 
        let val base  = RTensorSlice.base point 
            val shape = hd (RTensorSlice.shapes point)
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
                             
    fun size (SplayObj {(sz,lb,te), right, left}) = sz + (size left) + (size right)
      | size (SplayNil) = 0

    fun capacity (SplayObj {(sz,lb,te), right, left}) = 
        (hd (RTensor.shape te)) + (capacity left) + (capacity right)
      | capacity (SplayNil) = 0

    fun height (SplayObj {value=(sz,lb,te), right, left}) = 1 + Int.max(height left,height right)
      | height (SplayNil) = 0

    fun hasCapacityLeft P = (size P) < (capacity P)

    fun new (lb,cap) = SplayObj {value=(ref 0,lb,RTensor.new([cap,K])),right=SplayNil,left=SplayNil}

    fun insert (point, P) =
        if hasCapacityLeft P 
        then insert (point, P)
        else (let
                 val sz   = size P
                 val node = new (sz+1,nodeSize)
                 val P'   = splay(,join(P,node))
             in
                 
             end)

                 
end
*)

signature KDTREE = 
sig

type point
type pointSpace
type kdtree

structure S: KPOINT_SPACE

val pointSpace: kdtree -> S.pointSpace

val isEmpty: kdtree -> bool

val empty: kdtree

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


val nearestNeighbor: kdtree -> real list -> int option

val nearNeighbors: kdtree -> real -> real list -> int list

val remove: kdtree -> real -> real list -> kdtree

val kNearestNeighbors: kdtree -> int -> real list -> int list

val rangeSearch: kdtree -> (int * int) -> int list

val fromPoints: pointSpace -> kdtree

val addPoint: (int * int) -> real list * kdtree -> kdtree

end

functor KDTreeFn (
                  structure S : KPOINT_SPACE
                  val distance : (real list) * (real list) -> real
                 ): KDTREE = 
struct

structure S = S

type point = S.point
type pointSpace = S.pointSpace

datatype kdtree' =
         KdNode of { left: kdtree', i: int, right: kdtree', axis: int }
       | KdLeaf of { ii: int list, axis: int }

type kdtree = { P: pointSpace, T: kdtree' list }

exception Point
exception IndexArray

val K = S.K

fun pointSpace { P: pointSpace, T: kdtree' list } = P

fun compareDistance reltol probe (a,b) = 
    let
        val delta  = Real.- (distance (probe, S.pointList a), 
                             distance (probe, S.pointList b))
    in
        case reltol of
            NONE   => (if Real.< (delta, 0.0) then LESS else GREATER)
          | SOME r => (if Real.< (Real.* (delta,delta), r)
                       then EQUAL
                       else (if Real.< (delta, 0.0) 
                             then LESS else (if Real.> (delta, 0.0) then GREATER else EQUAL)))
    end

val empty = { P = S.empty, T = [KdLeaf { ii=[], axis=0 }] }

fun isEmpty' t =
    case t of
        KdNode _         => false
      | KdLeaf {ii,axis} => (List.length ii)=0

fun isEmpty {P,T} = List.all isEmpty' T


fun app f {P,T} =
    let
        val p = S.point P
        fun app' t =
            case t of 
                KdNode {left,i,right,axis} => (app' left; f (p i); app' right)
              | KdLeaf {ii,axis}           => List.app (fn i => f (p i)) ii
    in
        List.app app' T
    end


fun appi f {P,T} =
    let
        val p = S.point P
        fun appi' t =
            case t of 
                KdNode {left,i,right,axis} => (appi' left; f (i, p i); appi' right)
              | KdLeaf {ii,axis}           => List.app (fn i => f (i, p i)) ii
    in
        List.app appi' T
    end


fun foldl f init {P,T} =
    let
        val p = S.point P
        fun foldl' init t =
            case t of 
                KdNode {left,i,right,axis} => 
                (let val init' = foldl' init left
                     val init''= f (p i, init')
                 in foldl' init'' right end)
              | KdLeaf {ii,axis} => 
                List.foldl (fn (i,ax) => f (p i, ax)) init ii
    in
        List.foldl (fn (t, ax) => foldl' ax t) init T
    end


fun foldli f init {P,T} =
    let
        val p = S.point P
        fun foldli' init T =
            case T of 
                KdNode {left,i,right,axis} => 
                (let val init' = foldli' init left
                     val init''= f (i, p i, init')
                 in foldli' init'' right end)
              | KdLeaf {ii,axis}           => 
                List.foldl (fn (i,ax) => f (i, p i, ax)) init ii
    in
        List.foldl (fn (t, ax) => foldli' ax t) init T
    end
                                          

fun foldr f init {P,T} =
    let
        val p = S.point P
        fun foldr' init t =
            case t of 
                KdNode {left,i,right,axis} => 
                (let val init' = foldr' init right
                     val init''= f (p i, init')
                 in foldr' init'' left end)
              | KdLeaf {ii,axis} => 
                List.foldr (fn (i,ax) => f (p i, ax)) init ii
    in
        List.foldr (fn (t, ax) => foldr' ax t) init T
    end


fun foldri f init {P,T} =
    let
        val p = S.point P
        fun foldri' init t =
            case t of 
                KdNode {left,i,right,axis} => 
                (let val init'  = foldri' init right
                     val init'' = f (i, p i, init')
                 in foldri' init'' left end)
              | KdLeaf {ii,axis}  => 
                List.foldr (fn (i,ax) => f (i, p i, ax)) init ii
    in
        List.foldr (fn(t, ax) => foldri' ax t) init T
    end


fun ifoldr f init {P,T} =
    let
        fun foldr' init t =
            case t of 
                KdNode {left,i,right,axis} => 
                (let val init' = foldr' init right
                     val init''= f (i, init')
                 in foldr' init'' left end)
              | KdLeaf {ii,axis} => 
                List.foldr (fn (i,ax) => f (i, ax)) init ii
    in
        List.foldr (fn(t, ax) => foldr' ax t) init T
    end


fun size {P,T} = S.size P

fun toList t = foldr (op ::) [] t

fun toIndexList t = ifoldr (op ::) [] t

fun subtrees {P,T} =
    let
        fun subtrees' t ax = 
            case t of 
                KdNode {left,i,right,axis} => 
                (subtrees' left (t :: (subtrees' right ax)))
              | KdLeaf {ii,axis} => 
                t :: ax
    in
        List.foldl (fn(t, ax) => subtrees' t ax) [] T
    end

fun isValid' P t =
    case t of
        KdLeaf _ => true
      | KdNode {left,i,right,axis} =>
        let
            val x          = S.pointCoord P (i, axis)
            val leftValid  = (List.all (fn y => Real.< (S.coord y axis, x))
                                       (toList {P=P,T=[left]}))
            val rightValid = (List.all (fn y => Real.>= (S.coord y axis, x))
                                       (toList {P=P,T=[right]}))
        in
            leftValid andalso rightValid
        end

fun isValid {P,T} = List.all (isValid' P) T


(* Checks whether the K-D tree property holds for the given tree and its subtreees *)

fun allSubtreesAreValid (t as {P,T}) = 
    let
        val ts: kdtree' list = subtrees t
    in
        List.all (isValid' P) ts
    end


fun findiFromTo cmp (a,from,to) = IntArraySlice.findi cmp (IntArraySlice.slice (a,from,SOME (to-from+1)))

fun islice a (from,to) = IntArraySlice.slice (a,from,SOME (to-from+1))

(* Constructs a kd-tree from a point space, starting with the given depth. 
   If I is given, then only use the point indices contained in it, otherwise use all points. *)

fun fromPointsWithDepth P I depth =
    let 
        val sub          = IntArray.sub
        val pointCoord'  = S.pointCoord P
        val sz           = case I of NONE => S.size P
                                           | SOME I' => IntArray.length I'
        val bucketSize   = 10 * (Int.max (Real.ceil (Math.log10 (Real.fromInt sz)),  1))

        fun findMedian (I, m, n, depth) =
            let 
                fun findGreaterCoord (ci,cc,to,axis) =
                    findiFromTo 
                        (fn (i,x) => 
                            let 
                                val cx = pointCoord' (x, axis)
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


        fun fromPointsWithDepth' (I,m,n,depth) =
            (let
                 val k = n - m
             in
                 if (k <= bucketSize) orelse (k <= 1)
                 then 
                     let
                         val ii  = IntArraySlice.foldr (op ::) [] 
                                                       (IntArraySlice.slice (I, m, SOME (k+1)))
                     in
                         KdLeaf {ii=ii, axis=Int.mod(depth,K)}
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
                                                                                                                                        
                                   val left  = fromPointsWithDepth' (I',m,median-1,depth')
                                   val right = fromPointsWithDepth' (I',median+1,n,depth')
                               in
                                   KdNode {left=left,i=sub(I',median),
                                           right=right,axis=axis}
                               end)
                            | NONE => (KdLeaf {ii=IntArraySlice.foldr (op ::) []
                                                  (IntArraySlice.slice (I, m, SOME (k+1))),  
                                               axis=Int.mod (depth, K)})
                      end)
             end)
    in
        if sz=0 
        then KdLeaf {ii=[], axis=Int.mod (depth, K)}
        else (case I of NONE => fromPointsWithDepth' (IntArray.tabulate (sz, fn i => i), 0, sz-1, depth)
                      | SOME I' => if (IntArray.length I') <= sz 
                                   then fromPointsWithDepth' (I', 0, sz-1, depth)
                                   else raise IndexArray)
    end
    
   fun fromPoints P = {P=P,T=[(fromPointsWithDepth P NONE 0)]}

   (* Returns the index of the nearest neighbor of p in tree t. *)

   fun nearestNeighbor {P,T} probe =

       let
           val point'      = S.point P
           val pointCoord' = S.pointCoord P
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
                                   distance (probe, S.pointList candidate))
                       end

                   val candidates'' = if sphereIntersectsPlane
                                      then (case nearestNeighbor' t2 of
                                                SOME nn => candidates' @ [nn]
                                              | NONE => candidates')
                                      else candidates'

               in
                   ListUtils.minimumBy candidates'' compareDistance''
               end

           and nearestNeighbor' t =
               case t of

                   KdLeaf { ii, axis } => 
                   ListUtils.minimumBy ii compareDistance''
                   
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
           ListUtils.minimumBy (ListUtils.filterMap nearestNeighbor' T) compareDistance''
       end


   (* Returns all neighbors within distance r from p in tree t. *)
   fun nearNeighbors {P,T} radius probe =
       let
           val point'      = S.point P
           val pointCoord' = S.pointCoord P
           val r2          = Real.* (radius, radius)
           fun filterIndices ii = List.filter (fn (i) => (Real.<= (distance (probe, S.pointList (point' i)), r2))) ii

           fun nearNeighbors' t =
               case t of
                   
                   KdLeaf { ii, axis } => filterIndices ii
                                      
	         | KdNode { left, i, right, axis } =>
                   
	           (let 
                        val maybePivot = filterIndices [i]
	            in
		        if (isEmpty' left) andalso (isEmpty' right)
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
           List.foldl (fn (t,ax) => (nearNeighbors' t) @ ax) [] T
       end
  
  (* Removes the point p from t. *)
   fun remove  {P,T} tol pkill =

       let
           val point'      = S.point P
           val pointCoord' = S.pointCoord P
           val tol2        = Real.* (tol, tol)
           fun filterIndices ii = List.filter (fn (i) => (Real.> (distance (pkill, S.pointList (point' i)), tol2))) ii

	   fun remove' t =
               case t of
                   
                   KdLeaf { ii, axis } => 
                   KdLeaf { ii = filterIndices ii, axis = axis }
                                      
	         | KdNode { left, i, right, axis } =>
                   if (Real.> (distance (pkill, S.pointList (point' i)), tol2))
                   then 
                       (let
                            val I = IntArray.fromList ((toIndexList {P=P,T=[left]}) @ (toIndexList {P=P,T=[right]}))
                        in
                            fromPointsWithDepth P (SOME I) axis
                        end)
                   else
		       (if (Real.< (List.nth (pkill,axis), pointCoord' (i, axis)))
                        then 
                            KdNode { left = remove' left, i=i, right=right, axis=axis }
                        else 
                            KdNode { left = left, i=i, right=remove' right, axis=axis })
       in
           { P=P, T=List.map remove' T }
       end
  

  (* Returns the k nearest points to p within tree. *)
  fun kNearestNeighbors {P,T} k probe =
      let
          val point'      = S.point P
                            
          val compareDistance' = 
              compareDistance (SOME 1E~16) probe

          fun compareDistance'' (i,j) =
              compareDistance' (point' i, point' j)

          fun merge ([], ys) = ys
            | merge (xs, []) = xs
            | merge (l as x::xs, r as y:: ys) =
              case compareDistance'' (x,y) of
                  LESS => x::merge(xs, r)
                | _ => y::merge(l, ys)
              
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
                       recur ([], ii, k)
                   end)
                   
	        | KdNode { left, i, right, axis } =>
                  if (k <= 0) 
                  then []
                  else (let
                            val nearest = nearestNeighbor {P=P,T=[t]} probe
                        in
                            case nearest of
                                NONE => []
                              | SOME n => 
                                let
                                    val T' = remove {P=P,T=[t]} 1E~16 (S.pointList (point' n))
                                in
                                    n :: (kNearestNeighbors T' (k - 1) probe)
                                end
                        end))
      in
          List.take (List.foldl (fn (t,ax) => merge (kNearestNeighbors' t, ax)) [] T, k)
      end
                       

      fun inBounds (xMin, xMax) probe  =
          (Real.<= (probe,xMax) andalso Real.>= (xMin,probe))


      fun rangeSearch' (P,t) (bMin, bMax) ax =
          let
              val pointCoord' = S.pointCoord P
          in
              case t of
                  KdLeaf { ii, axis } => 
                  (let
                      val xMin = pointCoord' (bMin, axis)
                      val xMax = pointCoord' (bMax, axis)
                      val inBounds' = inBounds (xMin, xMax)
                  in
                      if Real.< (xMin, xMax) andalso 
                         Real.>= (xMax, pointCoord' (List.hd ii, axis))  andalso
                         Real.<= (xMin, pointCoord' (List.last ii, axis))
                      then (List.filter (fn (x) =>  inBounds' (pointCoord' (x, axis))) ii) @ ax
                      else ax
                  end)
                | KdNode { left, i, right, axis } =>
                  (let
                      val xMin = pointCoord' (bMin, axis)
                      val xMax = pointCoord' (bMax, axis)
                      val xVal = pointCoord' (i, axis)
                      val inBounds' = inBounds (xMin, xMax)
                  in
                      if Real.< (xMax, xVal)
                      then rangeSearch' (P,left) (bMin,bMax) ax
                      else (if Real.> (xMin, xVal)
                            then rangeSearch' (P,right) (bMin,bMax) ax
                            else (let val ax' = if inBounds' xVal then i::ax else ax
                                  in rangeSearch' (P,right) (bMin,bMax) 
                                                  (rangeSearch' (P,t) (bMin,bMax) ax')
                                  end))
                                                                
                  end)
          end

      fun rangeSearch {P,T} (bMin, bMax) = 
          List.foldl (fn(t,ax) => rangeSearch' (P,t) (bMin, bMax) ax) [] T


      fun joinTrees P (t,u) = 
          let
              val I = IntArray.fromList ((toIndexList {P=P,T=[t]}) @ (toIndexList {P=P,T=[u]}))
          in
              fromPointsWithDepth P (SOME I) 0
          end


      fun addPoint' (nodesize,leafsize) (P,T,j,n,lst) = 
          let
              val pointCoord   = S.pointCoord P
          in
              case T of
                  (KdLeaf {ii,axis})::ts =>
                  let
                      val t' = if List.null ii 
                               then KdLeaf {ii=[j],axis=axis}
                               else (let 
                                        val ii' = ListMergeSort.sort (fn(i,j) => pointCoord (i,axis) < pointCoord (j,axis)) (j :: ii)
                                        val len = List.length ii'
                                        fun remove (i, x::lst, ax) =
                                            if x=i then List.revAppend (ax, lst) else remove (i, lst, x::ax)
                                          | remove (i, [], ax) =
                                            List.rev ax 
                                    in
                                        if len > leafsize 
                                        then (let
                                                 val axis'   = Int.mod (axis+1, K)
                                                 val median  = List.nth (ii', Int.quot(len, 2))
                                                 val medianc = pointCoord (median, axis)
                                                 val (ii_left, ii_right) = List.partition (fn(i) => pointCoord (i,axis) < medianc) ii'
                                                 val left  = KdLeaf {ii=ii_left, axis=axis'}
                                                 val right = KdLeaf {ii=remove (median, ii_right, []), axis=axis'}
                                             in
                                                 KdNode {left=left,i=median,right=right,axis=axis}
                                             end)
                                        else (KdLeaf {ii=ii',axis=axis})
                                    end)
                  in
                      (t'::lst) @ ts
                  end
                | (t as KdNode {left,i,right,axis})::ts => 
                  if n > nodesize
                  then (case lst of 
                            u::rest => addPoint' (nodesize,leafsize) (P,[KdLeaf{ii=[],axis=0}],j,n-1,joinTrees P (t,u)::rest)
                          | [] => addPoint' (nodesize,leafsize) (P,[KdLeaf{ii=[],axis=0}],j,n,[t]))
                  else addPoint' (nodesize,leafsize) (P,ts,j,n+1,t::lst)
                | [] => (KdLeaf {ii=[j],axis=0}) :: lst
          end
                           
                           
      fun addPoint (nodesize,leafsize) (point, {P,T}) =
          let
              val (pidx,P') = S.insert (point, P)
          in
              {P=P',T=addPoint' (nodesize,leafsize) (P',T,pidx,1,[])}
          end

end
