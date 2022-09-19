
:- module(findnsols, [findnsols/4]).

:- meta_predicate
        findnsols(+, ?, :, -),
        findnsols(+, ?, :, -, ?).

findnsols(N, Template, Generator, List) :-
        findnsols(N, Template, Generator, List, []).

findnsols(N, Template, Generator, List, Tail) :-
        N > 0, !,
        findall(Template, maxsols(N, Generator), List, Tail).
findnsols(_, _, _, Tail, Tail).

maxsols(N, Generator) :-
        State = count(0),
        Generator,
        arg(1, State, C0),
        C1 is C0+1,
        (   C1 == N
        ->  !
        ;   nb_setarg(1, State, C1)
        ).
