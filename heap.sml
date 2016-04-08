(* Taken from Purely Functional Data Structures, by Chris Okasaki.           *)

signature HEAP =
sig

type 'a heap

val empty    : ('a * 'a -> order) -> 'a heap
val add      : 'a -> 'a heap -> 'a heap
val isEmpty  : 'a heap -> bool
val top      : 'a heap -> 'a            (* raises Empty *)
val remove   : 'a heap -> 'a * 'a heap  (* raises Empty *)
val numItems : 'a heap -> int
val app      : ('a -> unit) -> 'a heap -> unit
val toList   : 'a heap -> 'a list

end

structure Heap :> HEAP =
struct

datatype 'a node = E | T of int * 'a * 'a node * 'a node;

datatype 'a heap = Heap of ('a * 'a -> order) * int * 'a node;

fun rank E = 0
  | rank (T (r, _, _, _)) = r;

fun makeT (x, a, b) =
  if rank a >= rank b then T (rank b + 1, x, a, b) else T (rank a + 1, x, b, a);

fun merge f =
  let
    fun mrg (h, E) = h
      | mrg (E, h) = h
      | mrg (h1 as T (_, x, a1, b1), h2 as T (_, y, a2, b2)) =
      (case f (x, y) of GREATER => makeT (y, a2, mrg (h1, b2))
       | _ => makeT (x, a1, mrg (b1, h2)))
  in
    mrg
  end;

fun empty f = Heap (f, 0, E);

fun add x (Heap (f, n, a)) = Heap (f, n + 1, merge f (T (1, x, E, E), a));

fun isEmpty (Heap (_, _, E)) = true
  | isEmpty (Heap (_, _, T _)) = false;

fun top (Heap (_, _, E)) = raise Empty
  | top (Heap (_, _, T (_, x, _, _))) = x;

fun remove (Heap (_, _, E)) = raise Empty
  | remove (Heap (f, n, T (_, x, a, b))) = (x, Heap (f, n - 1, merge f (a, b)));

fun numItems (Heap (_, n, _)) = n;

fun app f =
  let
    fun ap [] = ()
      | ap (E :: rest) = ap rest
      | ap (T (_, d, a, b) :: rest) = (f d; ap (a :: b :: rest))
  in
    fn Heap (_, _, a) => ap [a]
  end;

fun toList' res h =
    if isEmpty h then rev res
    else let val (x, h) = remove h in toList' (x :: res) h end;

fun toList h = toList' [] h

end


structure MedianHeap =
struct

structure H = Heap

datatype 'a heap = { maxh: 'a H.heap, minh: 'a H.heap }

fun empty f = {maxh=Heap (maxf f, 0, E), minh=Heap (f, 0, E)}

fun add x (h as {maxh, minh}) = 
    let
        val h' = case H.compare (x,maxh) of
                     GREATER => {maxh=maxh, minh=H.add (x,minh)}
                   | _ => {maxh=H.add (x,maxh), minh=minh}
    in
        balance h'
    end
        
        


end
