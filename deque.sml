

signature QUEUE =
sig
  type 'a Queue

  val empty   : 'a Queue
  val isEmpty : 'a Queue -> bool

  val enqueue  : 'a * 'a Queue  -> 'a Queue
  val hd    : 'a Queue -> 'a         (* raises Empty if queue is empty *)
  val tl    : 'a Queue -> 'a Queue   (* raises Empty if queue is empty *)

  val map    : ('a -> 'b) -> 'a Queue -> 'b Queue   (* raises Empty if queue is empty *)
end


structure HoodMelvilleQueue : QUEUE =
struct
  datatype 'a RotationState =
      Idle
    | Reversing of int * 'a list * 'a list * 'a list * 'a list
    | Appending of int * 'a list * 'a list
    | Done of 'a list

  type 'a Queue = int * 'a list * 'a RotationState * int * 'a list

  fun exec (Reversing (ok, x :: f, f', y :: r, r')) =
      Reversing (ok+1, f, x :: f', r, y :: r')
    | exec (Reversing (ok, [], f', [y], r')) = Appending (ok, f', y :: r')
    | exec (Appending (0, f', r')) = Done r'
    | exec (Appending (ok, x :: f', r')) = Appending (ok-1, f', x :: r')
    | exec state = state

  fun invalidate (Reversing (ok, f, f', r, r')) =
      Reversing (ok-1, f, f', r, r')
    | invalidate (Appending (0, f', x :: r')) = Done r'
    | invalidate (Appending (ok, f', r')) = Appending (ok-1, f', r')
    | invalidate state = state

  fun exec2 (lenf, f, state, lenr, r) =
        case exec (exec state) of
            Done newf => (lenf, newf, Idle, lenr, r)
	  | newstate => (lenf, f, newstate, lenr, r)

  fun check (q as (lenf, f, state, lenr, r)) =
        if lenr <= lenf then exec2 q
        else let val newstate = Reversing (0, f, [], r, [])
             in exec2 (lenf+lenr, f, newstate, 0, []) end

  val empty = (0, [], Idle, 0, [])
  fun isEmpty (lenf, f, state, lenr, r) = (lenf = 0)

  fun enqueue (x, (lenf, f, state, lenr, r)) = 
      check (lenf,f,state,lenr+1,x::r)

  fun hd (lenf, [], state, lenr, r) = raise Empty
    | hd (lenf, x :: f, state, lenr, r) = x

  fun tl (lenf, [], state, lenr, r) = raise Empty
    | tl (lenf, x :: f, state, lenr, r) =
      check (lenf-1, f, invalidate state, lenr, r)

  fun map1 f que ax =
      if isEmpty que then ax else map1 f (tl que) (enqueue (f (hd que), ax))

  fun map f que = map1 f que empty

                       
end