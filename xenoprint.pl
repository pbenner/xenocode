% -*- mode: prolog; -*-

:- module(xenoprint,
	  [ xenoprint/2
	  , xenoprint/3
	  ]).

:- multifile xenoprint/2.
:- multifile xenoprint/3.

xenoprint(string(Str), A, B) :- !,
	new_memory_file(Handle),
	open_memory_file(Handle, write, S1),
	xenoprint(S1, A, B),
	close(S1),
	open_memory_file(Handle, read, S2),
	read_string(S2, _, Str),
	close(S2).
