% -*- mode: prolog; -*-

:- module(pairings_io,
	  [ xenoprint/2
	  , xenoprint/3
          ]).

:- use_module(xenoprint).
:- use_module(pairings).

% interface
% ------------------------------------------------------------------------------

xenoprint:xenoprint(pairings, PAIRS) :- !,
	print_pairs(current_output, PAIRS).

xenoprint:xenoprint(STREAM, pairings, PAIRS) :- !,
	print_pairs(STREAM, PAIRS).

% ------------------------------------------------------------------------------

print_pairs(STREAM, PAIRS) :-
	format(STREAM, "\nAC  | Codons\n----|---------------", []),
	forall(member((AA, C_LIST), PAIRS),
	       (format(STREAM, "\n~w | ", [AA]),
		member((AA, C_LIST), PAIRS),
		forall(member(C, C_LIST),
		       format(STREAM, "~w ", [C])))),
	format(STREAM, '\n', []).
