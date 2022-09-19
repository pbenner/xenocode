% -*- mode: prolog; -*-

:- module(xenocode_io,
	  [ xenoprint/2
	  , xenoprint/3
          ]).

:- use_module(xenoprint).
:- use_module(genetic_code).

% interface
% ------------------------------------------------------------------------------

xenoprint:xenoprint(genetic_code, book(BOOK)) :- !,
	implements_code(BOOK, CODE),
	print_code(current_output, CODE).

xenoprint:xenoprint(genetic_code, code(CODE)) :- !,
	print_code(current_output, CODE).

xenoprint:xenoprint(STREAM, genetic_code, book(BOOK)) :- !,
	implements_code(BOOK, CODE),
	print_code(STREAM, CODE).

xenoprint:xenoprint(STREAM, genetic_code, code(CODE)) :- !,
	print_code(STREAM, CODE).

% ------------------------------------------------------------------------------

% ?- create_codebook(ecoli, B), generate_code(2, B, CODE).
print_title(STREAM, N,SIZE) :-
	format(STREAM, "\n# "),
	format(STREAM, "~w reachable codes in ~w", [SIZE, N]),
	(N>1 -> format(STREAM, " steps", []); format(STREAM, " step", [])),
	format(STREAM, "\n").

print_code(STREAM, CODE) :-
	format(STREAM, "\n## Code", []),
	format(STREAM, "\nAA  | Codons\n----|---------------", []),
	forall(member((AA, C_LIST), CODE),
	       (format(STREAM, "\n~w | ", [AA]),
		member((AA,C_LIST),CODE),
		forall(member(C,C_LIST),
		       (format(STREAM, "~w ", [C]))))),
	format(STREAM, '\n', []).

% TODO: implement stream or rename predicate to avoid
% name clashes with print_code/2
print_code(_STREAM, CODE,CODE_NEW) :-
	format("\n## Code"),
	format("\nAA  | Codons\n----|---------------"),
	forall((member((AA,CODONS_OLD),CODE),
		member((AA,CODONS_NEW),CODE_NEW),
		intersection(CODONS_OLD,CODONS_NEW,COMMON),
		subtract(CODONS_NEW,CODONS_OLD,GAINED),
		subtract(CODONS_OLD,CODONS_NEW,LOST),
		format('\n'),format(AA),format(' | '),
		forall(member(C,LOST),
		       (format("<s>"), format(C), format('</s> '))),
		forall(member(C,GAINED),
		       (format("__"), format(C), format('__ '))),
		forall(member(C,COMMON),
		       (format(C), format(' '))))),
	format('\n').

print_code_changes(CODE,CODE_NEW) :-
	format("\n## Code"),
	format("\nAA  | Codons\n----|---------------"),
	forall((member((AA,CODONS_OLD),CODE),
		member((AA,CODONS_NEW),CODE_NEW),
		intersection(CODONS_OLD,CODONS_NEW,COMMON),
		subtract(CODONS_NEW,CODONS_OLD,GAINED),
		subtract(CODONS_OLD,CODONS_NEW,LOST),
		length(GAINED,NGAIN),
		length(LOST,NLOST),
		((NGAIN=0,NLOST=0) -> true;
		(format('\n'),format(AA),format(' | '),
		forall(member(C,LOST),
		       (format("<s>"), format(C), format('</s> '))),
		forall(member(C,GAINED),
		       (format("__"), format(C), format('__ '))),
		forall(member(C,COMMON),
		       (format(C), format(' '))))))),
	format('\n').

diff_code(CODE1,CODE2,DIFF) :-
	setof((C,AA_OLD,AA_NEW),
	      (member((AA_OLD,CODONS_OLD),CODE1),
	       member(C,CODONS_OLD),
	       member((AA_NEW,CODONS_NEW),CODE2),
	       member(C,CODONS_NEW),
	       AA_OLD \= AA_NEW),
	      DIFF).
