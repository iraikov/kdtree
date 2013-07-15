
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

structure KDTree = KDTreeFn (val K = 3
                             val distanceSquared = distanceSquared3D)

functor KdTreeTestFn (val distance : (real list) * (real list) -> real) =
struct


   fun sortPoints (origin,pts) =
       ListMergeSort.sort
           (fn (x,y) => Real.> (distance (origin,x), distance (origin,y)))
           pts

           
   fun check (t)  = 
       (if not ((List.length (KDTree.toList t)) = (KDTree.size t)) then raise Fail "invalid KDTree size" else ();
        if not (KDTree.isValid t) then raise Fail "invalid KDTree" else ();
        if not (KDTree.allSubtreesAreValid t) then raise Fail "invalid subtree of a KDTree" else ())


   fun testNearestNeighbor (sorted,t as {P,T},x) =
       let
           val nn = valOf (KDTree.nearestNeighbor t x)
           val nn' = [RTensor.sub (P, [nn,0]),
                      RTensor.sub (P, [nn,1]),
                      RTensor.sub (P, [nn,2])]

           val _ = (print "nn' = "; TensorFile.realListLineWrite TextIO.stdOut nn')
           val _ = (print "distance(x,nn') = "; TensorFile.realWrite TextIO.stdOut (distance3D (x,nn')))
           val _ = (print "hd sorted = "; TensorFile.realListLineWrite TextIO.stdOut (hd sorted))
           val _ = (print "distance(x,hd sorted) = "; TensorFile.realWrite TextIO.stdOut (distance3D (x,(hd sorted))))

       in
           if not (ListPair.all (fn (x,y) => Real.>= (1E~16, Real.- (x,y) )) (nn', (hd sorted)))
           then raise Fail "nearestNeighbor" else ()
       end
       

   fun testNearNeighbors (sorted,t as {P,T},x,r) =
       let
           val nns  = KDTree.nearNeighbors t r x
           val nns' = sortPoints (x, List.map (fn (nn) => [RTensor.sub (P, [nn,0]),
                                                           RTensor.sub (P, [nn,1]),
                                                           RTensor.sub (P, [nn,2])])
                                              nns)
           val (sss,_) = List.partition (fn (p) => Real.<= (distance3D (p,x), r))
                                        sorted
       in
           
           if not (ListPair.all
                       (fn (nn,ss) => (ListPair.all (fn (x,y) => Real.>= (1E~16, Real.- (x,y) )) (nn,ss)))
                       (nns', sss))
           then raise Fail "nearNeighbors" else ()
       end
       

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

val M  = 1000000
val P  = realRandomTensor (13,17) [M,N]

fun sortPoints distance (origin,pts) =
    ListMergeSort.sort
        (fn (x,y) => Real.> (distance (origin,x), distance (origin,y)))
        pts

val pts = let 
             fun recur (i, lst) =
                 if (i < M) 
                 then (let val p = [RTensor.sub (P, [i,0]),
                                    RTensor.sub (P, [i,1]),
                                    RTensor.sub (P, [i,2])] in
                           recur (i+1,p :: lst) end)
                 else lst
          in
              recur (0, [])
          end



val t  = KDTree.fromTensor P

val _ = print "tree constructed\n"

val _ = print ("tree size = " ^ (Int.toString (KDTree.size t)) ^ "\n")
val _ = print ("length of tree list = " ^ (Int.toString (List.length (KDTree.toList t))) ^ "\n")

val _  = KdTreeTest.check t

val _ = print "consistency check passed\n"

val seed   = Random.rand (19,21)

val Ntrials = 100

val _  = let  fun recur (i) =
               if i > 0
               then
                   (let
                        val xi = Int.mod (Random.randNat seed, N)
                        val x  = [Real.+ (0.1, RTensor.sub (P, [xi,0])),
                                  Real.- (RTensor.sub (P, [xi,1]), 0.1),
                                  Real.+ (0.1, RTensor.sub (P, [xi,2]))]
                                 
                        val sorted = sortPoints distance3D (x,pts)
                    in 
                        print ("test trial " ^ (Int.toString i) ^ "\n");
                        KdTreeTest.testNearestNeighbor (sorted,t,x);
                        print "nearest neighbor check passed\n";
                        KdTreeTest.testNearNeighbors (sorted,t,x,0.3);
                        print "near neighbors check passed\n";
                        recur (i-1)
                    end)
               else ()
         in
             recur (Ntrials)
         end

