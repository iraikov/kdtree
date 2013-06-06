
exception Point

fun distanceSquared3D ([x1,x2,x3], [y1,y2,y3]) = 
    let 
        val d1 =  Real.- (x1, y1)
        val d2 =  Real.- (x2, y2)
        val d3 =  Real.- (x3, y3)
    in
        List.foldl Real.+ 0.0 [Real.* (d1,d1), Real.* (d2,d2), Real.* (d3,d3)]
    end
  | distanceSquared3D (_, _) = raise Point

fun distance3D (x,y) = Math.sqrt(distanceSquared3D (x,y))

structure KDTree = KDTreeFn (val N = 3
                              val distanceSquared = distanceSquared3D)

functor KdTreeTestFn (val distance : (real list) * (real list) -> real) =
struct

   fun sortPoints (origin,pts) =
       ListMergeSort.sort
           (fn (x,y) => Real.> (distance (origin,x), distance (origin,y)))
           pts

           
   fun check (t)  = 
       (if not (KDTree.isValid t) then raise Fail "invalid KDTree" else ();
        if not (KDTree.allSubtreesAreValid t) then raise Fail "invalid subtree of a KDTree" else ())

end

fun realRandomTensor (xseed,yseed) shape =
    let 
        val length = Index.length shape
        val seed   = Random.rand (xseed,yseed)
        val a      = RTensor.Array.array(length, Random.randReal seed)
        fun loop 0 = RTensor.fromArray(shape, a)
          | loop j = (RTensor.Array.update(a, length-j, Random.randReal seed);
                      loop (j-1))
    in 
        loop (length - 1)
    end


structure KdTreeTest = KdTreeTestFn (val distance = distance3D)

val N  = 3

val P  = realRandomTensor (13,17) [1000000,N]

(*val _ = (print "P = "; TensorFile.realTensorWrite (TextIO.stdOut) P)*)

val t  = KDTree.fromTensor P

val _ = print "tree constructed\n"

val _  = KdTreeTest.check t
