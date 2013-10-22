
fun timing (action) = 
    let
        val timer = Timer.startCPUTimer ()
        val result = action ()
        val times = Timer.checkCPUTimer timer
    in
        (result, Time.+ (#usr times, #sys times))
    end

exception Point

val K = 3

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

fun pointStorageSize P = hd (RTensor.shape P)

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


structure KDTree = KDTreeFn (type point = RTensorSlice.slice
                             type pointStorage = RTensor.tensor
                             val point = point
                             val pointList = pointList
                             val coord = coord
                             val pointCoord = pointCoord
                             val pointStorageSize = pointStorageSize
                             val K = K
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



val (t,ti)  = timing (fn () => KDTree.make P)

val _ = print ("tree constructed (" ^ (Time.toString ti) ^ " s)\n")

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
                                 
                        val sorted  = sortPoints distance3D (x,pts)
                        val _       = print ("test trial " ^ (Int.toString i) ^ "\n")
                        val (_,ti)  = timing (fn () => KdTreeTest.testNearestNeighbor (sorted,t,x))
                        val _       = print ("nearest neighbor check passed (" ^ (Time.toString ti) ^ " s)\n")
                        val (_,ti)  = timing (fn () => KdTreeTest.testNearNeighbors (sorted,t,x,0.3))
                        val _       = print ("near neighbors check passed (" ^ (Time.toString ti) ^ "s)\n")
                    in 
                        recur (i-1)
                    end)
               else ()
         in
             recur (Ntrials)
         end

