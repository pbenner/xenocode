% -*- mode: prolog; -*-

:- module(xfreecode,
	  [ filter_codon_list/3
	  , xfree_code/4
	  ]).

:- use_module(codebook).
:- use_module(query).
:- use_module(xenocode).
:- use_module(pairings).

% generate x-free codes
% ------------------------------------------------------------------------------

xfree_codon(X, C, BOOK) :-
	(var(C) -> query([c], [C], BOOK) ; true),
	atom_chars(C, [C1, C2, C3]),
	C1 \= X, C2 \= X, C3 \= X.

xfree_anticodon(X, AC, TRNA, BOOK) :-
	(var(AC)
	 -> L = [a, c, g, u],
	    member(AC1, L),
	    member(AC2, L),
	    member(AC3, L),
     	    atom_chars(AC, [AC1, AC2, AC3])
	  ; true),
	once(((modifies_ac(_, TRNA, AC, ACM, BOOK) -> true; ACM = AC),
	      pairs_with(C, ACM, BOOK),
	      xfree_codon(X, C, BOOK))).

xfree_code(X, BOOK) :-
        forall(amino_acid(AA, BOOK),
	       (is_valid_code(C, _, _, AA, BOOK),
	        xfree_codon(X, C, BOOK))).

% TRNA has an anticodon that pairs with at least one x-free anticodon
xfree_trna(X, TRNA, BOOK) :-
	ac_on_trna(AC, TRNA, BOOK),
	once((is_valid_code(C, AC, TRNA, _, BOOK),
	      xfree_codon(X, C, BOOK))).

% remove all codons from the codebook that contain X
% ------------------------------------------------------------------------------

filter_codon_list_rec(_, [], []) :- !.
filter_codon_list_rec(X, [C | TAIL1], [C | TAIL2]) :-
	atom_chars(C, [C1, C2, C3]),
	C1 \= X, C2 \= X, C3 \= X, !,
	filter_codon_list_rec(X, TAIL1, TAIL2).
filter_codon_list_rec(X, [_ | TAIL1], TAIL2) :-
	filter_codon_list_rec(X, TAIL1, TAIL2).

filter_codon_list(X, BOOK1, BOOK2) :-
	get_chapter(codons, BOOK1, CHAPTER1),
	filter_codon_list_rec(X, CHAPTER1, CHAPTER2),
	replace_chapter(codons, BOOK1, CHAPTER2, BOOK2).

% ------------------------------------------------------------------------------

% get a list of TRNAs where the x-free ones are at the end
sorted_trna_list(X, BOOK, TRNA_LIST_NONXFREE, TRNA_LIST_XFREE, TRNA_LIST) :-
	setof(TRNA, AC^ac_on_trna(AC, TRNA, BOOK), TRNA_LIST_FULL),
	setof(TRNA,    xfree_trna( X, TRNA, BOOK), TRNA_LIST_XFREE),
	subtract(TRNA_LIST_FULL, TRNA_LIST_XFREE, TRNA_LIST_NONXFREE),
	append(TRNA_LIST_NONXFREE, TRNA_LIST_XFREE, TRNA_LIST).
sorted_trna_list(X, BOOK, TRNA_LIST) :-
	sorted_trna_list(X, BOOK, _, _, TRNA_LIST).

% ------------------------------------------------------------------------------

% keep only one tRNA for each amino acid
xfree_code_delete_trnas_rec(_,        [],           [],     [], BOOK1, BOOK1) :-
	!.
xfree_code_delete_trnas_rec(X, TRNA_LIST, FREE_AC_LIST, RECIPE, BOOK1, BOOK3) :-
	% N: number of deleted tRNAs
	% M: number of released anti-codons
	TRNA_LIST = [TRNA | TAIL1],
	% get the AA associated with TRNA
	aa_is_loaded_on(AA, TRNA, _, BOOK1),
	% count the number of TRNAs that carry AA
	setof((AA,T), AC_^T^aa_is_loaded_on(AA, T, AC_, BOOK1), LIST),
	length(LIST, L),
	% delete TRNA if there are others that carry AA
	(L > 1 -> delete_trna(TRNA, BOOK1, BOOK2),
		  RECIPE = [delete(TRNA) | TAIL3],
	          format("deleting tRNA ~w\n", TRNA),
		  % if the AC associated with TRNA pairs with and x-free codon,
	          % append AC to the FREE_AC_LIST
		  ((is_valid_code(C, AC, TRNA, _, BOOK1), xfree_codon(X, C, BOOK1))
		   -> FREE_AC_LIST = [ AC | TAIL2], format("tRNA has valid anticodon: ~w\n", AC)
 		    ; FREE_AC_LIST = TAIL2)
	        ; BOOK1 = BOOK2,
		  RECIPE = TAIL3,
		  FREE_AC_LIST = TAIL2),
	xfree_code_delete_trnas_rec(X, TAIL1, TAIL2, TAIL3, BOOK2, BOOK3).

xfree_code_delete_trnas(X, FREE_AC_LIST, RECIPE, BOOK1, BOOK2) :-
	sorted_trna_list(X, BOOK1, TRNA_LIST),
	xfree_code_delete_trnas_rec(X, TRNA_LIST, FREE_AC_LIST, RECIPE, BOOK1, BOOK2).

% ------------------------------------------------------------------------------

% deactivate synthasis if it requires anticodons that pair
% ONLY with codons containing X
xfree_code_deactivate_synth_rec(_,         [],     [], BOOK1, BOOK1) :- !.
xfree_code_deactivate_synth_rec(X, SYNTH_LIST, RECIPE, BOOK1, BOOK2) :-
	SYNTH_LIST = [SYNTH | TAIL1],
	recognizes(SYNTH, _, AC, BOOK1),
	% SYNTH is already deactive
	var(AC), !,
	xfree_code_deactivate_synth_rec(X, TAIL1, RECIPE, BOOK1, BOOK2).
xfree_code_deactivate_synth_rec(X, SYNTH_LIST, RECIPE, BOOK1, BOOK3) :-
	SYNTH_LIST = [SYNTH | TAIL1],
	(forall((recognizes(SYNTH, _, AC, BOOK1),
		 pairs_with(C, AC, BOOK1)),
		\+ xfree_codon(X, C, BOOK1))
	 -> deactivate_synth(SYNTH, BOOK1, BOOK2), RECIPE = [deactivate(SYNTH) | TAIL2], format("deactivating ~w\n", SYNTH)
	  ; BOOK1 = BOOK2, RECIPE = TAIL2),
	xfree_code_deactivate_synth_rec(X, TAIL1, TAIL2, BOOK2, BOOK3).

xfree_code_deactivate_synth(X, RECIPE, BOOK1, BOOK2) :-
	setof(SYNTH, AA^aa_is_loaded_by(AA, SYNTH, BOOK1), SYNTH_LIST),
	xfree_code_deactivate_synth_rec(X, SYNTH_LIST, RECIPE, BOOK1, BOOK2).

% ------------------------------------------------------------------------------

% mutate the given TRNA
xfree_code_mutate_ac_(X, AC, TRNA, BOOK1, BOOK2) :-
	xfree_anticodon(X, AC, TRNA, BOOK1),
	% AC should not be part of the codebook
	\+ query([ac], [AC], BOOK1),
	% mutate the TRNA and make sure the code is not
	% ambiguous
	mutate_trna_anticodon(_, AC, TRNA, BOOK1, BOOK2),
	% consistency checks
	once(is_valid_code(_, AC, TRNA, _, BOOK2)),
	implements_code(BOOK2, CODE),
	\+ code_is_ambiguous(CODE).

xfree_code_mutate_ac(X, AC, TRNA, RECIPE, BOOK1, BOOK2) :-
	xfree_code_mutate_ac_(X, AC, TRNA, BOOK1, BOOK2),
	RECIPE = [mutate(TRNA, AC)].
	
% ------------------------------------------------------------------------------

% mutate anticodons to pair with x-free codons
xfree_code_rec(X, _, [], BOOK1, BOOK1) :-
	% terminate if the code is x-free
	print_code(BOOK1),
	xfree_code(X, BOOK1),
	!.
xfree_code_rec(X, AC_LIST, TRNA_LIST, RECIPE, BOOK1, BOOK3) :-
	format("TRNA LIST: ~w\n", [TRNA_LIST]),
	% make sure we have enough free anticodons
	length(TRNA_LIST, L1),
	length(  AC_LIST, L2),
	L2 >= L1, !,
	% select an anticodon
	select(AC, AC_LIST, TAIL1),
	% get the next TRNA
	TRNA_LIST  = [TRNA | TAIL2],
	% mutate the anticodon on TRNA to match C
	xfree_code_mutate_ac(X, AC, TRNA, RECIPE_TMP, BOOK1, BOOK2),
	% add mutation to the recipe
	append(RECIPE_TMP, TAIL3, RECIPE),
	% recursive call
	xfree_code_rec(X, TAIL1, TAIL2, TAIL3, BOOK2, BOOK3).

xfree_code(X, RECIPE, BOOK1, BOOK5) :-
	var(X) -> fail;
	xfree_code_deactivate_synth(X, RECIPE1, BOOK1, BOOK2),
	% keep one tRNA for each AA
	xfree_code_delete_trnas(X, RELEASED_AC_LIST, RECIPE2, BOOK2, BOOK3),
	% get a list of the remaining tRNAs
	sorted_trna_list(X, BOOK3, TRNA_LIST_NONXFREE, _, _TRNA_LIST),
	format("released x-free anticodons  : ~w\n", [RELEASED_AC_LIST]),
	format("tRNAs that require mutations: ~w\n", [TRNA_LIST_NONXFREE]), !,
	filter_codon_list(X, BOOK3, BOOK4),
	setof(AC, TRNA^ac_on_trna(AC, TRNA, BOOK4), USED_AC_LIST),
	complement_ac_list(USED_AC_LIST, FREE_AC_LIST, BOOK4),
	format("free ac list: ~w\n", [FREE_AC_LIST]),
	% mutate tRNA anticodons (use only those tRNAs with non-xfree anticodons)
	xfree_code_rec(X, FREE_AC_LIST, TRNA_LIST_NONXFREE, RECIPE3, BOOK4, BOOK5),
	% construct the full recipe
	flatten([RECIPE1, RECIPE2, RECIPE3], RECIPE).

% ------------------------------------------------------------------------------

/** <examples>

?- create_codebook(ecoli, B1), xfree_code(c, RECIPE, B1, B2), print_code(B2), writeln(RECIPE).

*/
