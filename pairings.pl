% -*- mode: prolog; -*-

:- module(pairings, [pairings/2]).

:- use_module(alphabet).
:- use_module(listutil).
:- use_module(query).

% Generate unambiguous (C, AC) pairs, considering possible AC modifications
% but ignore whether or not a suitable tRNA exists
% ------------------------------------------------------------------------------

pairs_with_modification(C, AC, BOOK) :-
	% bypass the query interface in order to ignore the
	% tRNA database
	query:pairs_with(C, ACM, BOOK),
	\+ query:stop_codon(C, BOOK),
	query:matures(foobar, AC, ACM, BOOK),
	nucleobase_triplet(rna, AC).

codons_pairing_with_ac(AC, CODONS, BOOK) :-
	setof(C, pairs_with_modification(C, AC, BOOK), CODONS).

generate_ac_list(AC_LIST, BOOK) :-
	setof(AC, C^pairs_with_modification(C, AC, BOOK), AC_LIST).

remove_ambiguous(    [],     []) :- !.
remove_ambiguous(PAIRS1, PAIRS3) :-
	PAIRS1 = [(AA1, CODON_LIST1) | TAIL],
	member(CODON, CODON_LIST1),
	member((AA2, CODON_LIST2), TAIL),
	member(CODON, CODON_LIST2), !,
	% remove one of the anticodons
	( delete((AA1, CODON_LIST1), PAIRS1, PAIRS2)
	; delete((AA2, CODON_LIST2), PAIRS1, PAIRS2)),
	remove_ambiguous(PAIRS2, PAIRS3).
remove_ambiguous(PAIRS1, PAIRS2) :-
	PAIRS1 = [H | T1],
	PAIRS2 = [H | T2],
	remove_ambiguous(T1, T2).

% generates a list with elements (AC, CODON_LIST), where each
% C in CODON_LIST pairs with AC
pairings(PAIRS2, BOOK) :-
	% get list of anticodons
	generate_ac_list(AC_LIST, BOOK),
	% for each anticodon AC get the list of valid codons
	setof((AC, CODONS), (member(AC, AC_LIST),
			     codons_pairing_with_ac(AC, CODONS, BOOK)),
	      PAIRS1),
	remove_ambiguous(PAIRS1, PAIRS2).

/** <examples>

% remove all codons containing a c and compute all unambiguous valid pairs,
% print all lengths of PAIRS (i.e. the number of amino acids that could in
% principle be coded)
?- create_codebook(ecoli, B1),
	findall(L, (xfreecode:filter_codon_list(c, B1, B2),
		    pairings(PAIRS, B2),
		    length(PAIRS, L)),
		List),
	list_to_set(List, S),
	writeln(S).

*/

% generate a set of valid anticodons to complement the x-free code
% ------------------------------------------------------------------------------

complement_ac_list(USED_AC_LIST, FREE_AC_LIST, BOOK) :-
	% get list of anticodons
	generate_ac_list(AC_LIST, BOOK),
	% for each anticodon AC get the list of valid codons
	setof((AC, CODONS), (member(AC, AC_LIST),
			     member(AC, USED_AC_LIST),
			     codons_pairing_with_ac(AC, CODONS, BOOK)),
	      USED_PAIRS),
	setof((AC, CODONS), (member(AC, AC_LIST),
			     \+ member(AC, USED_AC_LIST),
			     codons_pairing_with_ac(AC, CODONS, BOOK)),
	      FREE_PAIRS1),
	setof((AC, CODONS1), (member((AC, CODONS1), FREE_PAIRS1),
			      % the codes should not overlap
			      \+ (member((_, CODONS2), USED_PAIRS),
				  member(C, CODONS1),
				  member(C, CODONS2))),
	      FREE_PAIRS2),
	remove_ambiguous(FREE_PAIRS2, FREE_PAIRS3),
	setof(AC, CODONS^member((AC, CODONS), FREE_PAIRS3), FREE_AC_LIST).
