
fun timing (action) = 
    let
        val timer = Timer.startCPUTimer ()
        val result = action ()
        val times = Timer.checkCPUTimer timer
    in
        (result, Time.+ (#usr times, #sys times))
    end

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

val _ = print "starting KDTree tests\n"

structure TensorKDTree = KDTreeFn (structure S = TensorPointSpaceFn (val K = 3)
                                   val distance = distanceSquared3D)
structure MapKDTree = KDTreeFn (structure S = MapPointSpaceFn (val K = 3)
                                val distance = distanceSquared3D)

functor KdTreeTestFn (structure KDTree: KDTREE
                      val distance : (real list) * (real list) -> real) =
struct


   fun sortPoints (origin,pts) =
       ListMergeSort.sort
           (fn (x,y) => Real.> (distance (origin,x), distance (origin,y)))
           pts

           
   fun check (t)  = 
       (if not ((List.length (KDTree.toList t)) = (KDTree.size t)) then raise Fail "invalid KDTree size" else ();
        if not (KDTree.isValid t) then raise Fail "invalid KDTree" else ();
        if not (KDTree.allSubtreesAreValid t) then raise Fail "invalid subtree of a KDTree" else ())


   fun testNearestNeighbor (sorted,t,x) =
       let
           val nn = valOf (KDTree.nearestNeighbor t x)
           val nn' = KDTree.S.pointList (KDTree.S.point (KDTree.pointSpace t) nn)

           val _ = (print "nn' = "; TensorFile.realListLineWrite TextIO.stdOut nn')
           val _ = (print "distance(x,nn') = "; TensorFile.realWrite TextIO.stdOut (distance3D (x,nn')))
           val _ = (print "hd sorted = "; TensorFile.realListLineWrite TextIO.stdOut (hd sorted))
           val _ = (print "distance(x,hd sorted) = "; TensorFile.realWrite TextIO.stdOut (distance3D (x,(hd sorted))))

       in
           if not (ListPair.all (fn (x,y) => Real.>= (1E~16, Real.- (x,y) )) (nn', (hd sorted)))
           then raise Fail "nearestNeighbor" else ()
       end
       

   fun testNearNeighbors (sorted,t,x,r) =
       let
           val P = KDTree.pointSpace t
           val nns  = KDTree.nearNeighbors t r x
           val nns' = sortPoints (x, List.map (fn (nn) => KDTree.S.pointList (KDTree.S.point P nn))
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


structure TensorKdTreeTest = KdTreeTestFn (structure KDTree = TensorKDTree
                                           val distance = distance3D)

structure MapKdTreeTest = KdTreeTestFn (structure KDTree = MapKDTree
                                        val distance = distance3D)

val M  = 1000000
val N  = 3


val _ = print ("constructing tensor point space...\n")
val P  = realRandomTensor (13,17) [M,N]

val _ = print ("constructing map point space...\n")
val (PM,ti) = timing (fn () =>
                         Loop.foldi (0, M, fn (i,pm) => 
                                              let 
                                                  val p = [RTensor.sub (P, [i,0]),
                                                           RTensor.sub (P, [i,1]),
                                                           RTensor.sub (P, [i,2])]
                                              in 
                                                  #2(MapKDTree.S.insert(p,pm) )
                                              end,
                                     MapKDTree.S.empty))
val _ = print ("map point space constructed (" ^ (Time.toString ti) ^ " s)\n")



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



val _ = print ("constructing tensor tree...\n")
val (tt,ti)  = timing (fn () => TensorKDTree.fromPoints P)
val _ = print ("tensor tree constructed (" ^ (Time.toString ti) ^ " s)\n")

val _ = print ("constructing map tree...\n")
val (mt,ti)  = timing (fn () => MapKDTree.fromPoints PM)
val _ = print ("map tree constructed (" ^ (Time.toString ti) ^ " s)\n")


val _ = print ("tensor tree size = " ^ (Int.toString (TensorKDTree.size tt)) ^ "\n")
val _ = print ("map tree size = " ^ (Int.toString (TensorKDTree.size tt)) ^ "\n")

val _ = print ("length of tree list = " ^ (Int.toString (List.length (TensorKDTree.toList tt))) ^ "\n")

val _  = TensorKdTreeTest.check tt
val _ = print "tensor tree consistency check passed\n"

val _  = MapKdTreeTest.check mt
val _ = print "map tree consistency check passed\n"

val seed   = Random.rand (19,21)

val Ntrials = 2

val _  = let  fun recur (i) =
               if i > 0
               then
                   (let
                        val xi = Int.mod (Random.randNat seed, N)
                        val x  = [Real.+ (0.1, RTensor.sub (P, [xi,0])),
                                  Real.- (RTensor.sub (P, [xi,1]), 0.1),
                                  Real.+ (0.1, RTensor.sub (P, [xi,2]))]
                                 
                        val sorted  = sortPoints distance3D (x,pts)
                        val _       = print ("test trial " ^ (Int.toString i) ^ "\n")
                        val (_,ti)  = timing (fn () => TensorKdTreeTest.testNearestNeighbor (sorted,tt,x))
                        val _       = print ("(tensor tree) nearest neighbor check passed (" ^ (Time.toString ti) ^ " s)\n")
                        val (_,ti)  = timing (fn () => MapKdTreeTest.testNearestNeighbor (sorted,mt,x))
                        val _       = print ("(map tree) nearest neighbor check passed (" ^ (Time.toString ti) ^ " s)\n")
                        val (_,ti)  = timing (fn () => TensorKdTreeTest.testNearNeighbors (sorted,tt,x,0.3))
                        val _       = print ("(tensor tree) near neighbors check passed (" ^ (Time.toString ti) ^ "s)\n")
                        val (_,ti)  = timing (fn () => MapKdTreeTest.testNearNeighbors (sorted,mt,x,0.3))
                        val _       = print ("(map tree) near neighbors check passed (" ^ (Time.toString ti) ^ "s)\n")
                    in 
                        recur (i-1)
                    end)
               else ()
         in
             recur (Ntrials)
         end

