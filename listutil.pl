% -*- mode: prolog; -*-

:- module(listutil,
	  [ replace/4
	  , delete/3
	  , insert/3 ]).

% ------------------------------------------------------------------------------

replace(_, _, [], []) :- !.
replace(O1, R, [O2|T1], [R|T2]) :-
	copy_term(O1, O1C),
	copy_term(O2, O2C),
	O1C = O2C,
	replace(O1, R, T1, T2), !.
replace(O, R, [H|T1], [H|T2]) :-
	replace(O, R, T1, T2).

delete(_, [], []) :- !.
delete(O1, [O2|T],    T2 ) :-
	copy_term(O1, O1C),
	copy_term(O2, O2C),
	O1C = O2C,
	delete(O1, T, T2), !.
delete(O, [H|T1], [H|T2]) :-
	delete(O, T1, T2).

insert(O, T, [O|T]).
