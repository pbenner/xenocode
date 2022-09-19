% -*- mode: prolog; -*-

:- module(alphabet,
	  [ nucleobase/2
	  , nucleobase_triplet/2 ]).

% nucleotide definitions on DNA and RNA
% ------------------------------------------------------------------------------

nucleobase(dna, a).
nucleobase(dna, c).
nucleobase(dna, g).
nucleobase(dna, t).

nucleobase(rna, a).
nucleobase(rna, c).
nucleobase(rna, g).
nucleobase(rna, u).

nucleobase_triplet(Type, T) :-
	(nonvar(T) -> atom_chars(T, [T1, T2, T3]) ; true),
	nucleobase(Type, T1),
	nucleobase(Type, T2),
	nucleobase(Type, T3),
	atom_chars(T, [T1, T2, T3]).
