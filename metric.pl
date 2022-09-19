% -*- mode: prolog; -*-

:- module(metric, [d/3]).

:- use_module(genetic_code).

% ------------------------------------------------------------------------------

d_rec(                  [],                   [], N1, N1) :- !.
d_rec([(AC, CODONS1) | T1], [(AC, CODONS2) | T2], N1, N3) :- !,
	subtract(CODONS1, CODONS2, List1),
	subtract(CODONS2, CODONS1, List2),
	length(List1, L1),
	length(List2, L2),
	N2 is (L1 + L2)/2.0 + N1,
	d_rec(T1, T2, N2, N3).
% otherwise fail
d_rec(_, _, _, error("Invalid genetic code!")).

d(V1, V2, N) :-
	(V1 = book(B1) -> implements_code(B1, CODE1) ; CODE1 = V1),
	(V2 = book(B2) -> implements_code(B2, CODE2) ; CODE2 = V2),
	d_rec(CODE1, CODE2, 0, N).

/** <examples>

?- create_codebook(ecoli, B1), implements_code(B1, C1), generate_code(1, B1, C2), d(C1, C2, N).

*/
