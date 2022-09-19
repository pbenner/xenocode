% -*- mode: prolog; -*-

:- module(xenocode,
	  [ delete_trna/3
          , delete_trna/4
          , insert_trna/4
	  , mutate_trna_anticodon/5
	  , deactivate_synth/3
	  , recode/5
	  , recode/6
	  , generate_code/4
	  , generate_code/5
	  , generate_codebook/4
	  , generate_codebook/5
	  , generate_codes/4
	  , generate_codes/5
	  , generate_codes_mt/4
          ]).

:- use_module(query).
:- use_module(codebook).
:- use_module(genetic_code).
:- use_module(listutil).
:- use_module(library(thread)).

% tools for altering the genetic code
% ------------------------------------------------------------------------------

% true if the two codons have different nucleotides at
% all three positions
pairwise_distinct(C1, C2) :-
	atom_chars(C1, [N1, N2, N3]),
	atom_chars(C2, [M1, M2, M3]),
	N1 \= M1, N2 \= M2, N3 \= M3.

pairwise_distinct_2(C1, C2) :-
	atom_chars(C1, [N1, N2, _]),
	atom_chars(C2, [M1, M2, _]),
	N1 \= M1, N2 \= M2.

delete_trna(AC, TRNA, BOOK1, BOOK3) :-
	delete_entry(ac_on_trna, (AC, TRNA), BOOK1, BOOK2),
	(is_clone(AC, TRNA, BOOK2)
	 -> delete_entry(trna_clones, (AC, TRNA), BOOK2, BOOK3)
	  ; BOOK2 = BOOK3).

delete_trna(TRNA, BOOK1, BOOK2) :-
	delete_trna(_, TRNA, BOOK1, BOOK2) .

insert_trna(AC, TRNA, BOOK1, BOOK2) :-
	\+ query([ac, trna], [AC, TRNA], BOOK1),
	insert_entry(ac_on_trna, (AC, TRNA), BOOK1, BOOK2).

% take a tRNA and add a clone to the codebook with a new
% anticodon AC
insert_cloned_trna(AC, TRNA, BOOK1, BOOK4) :-
        forall(query([ac, trna], [AC_OLD, TRNA], BOOK1),
	       AC \= AC_OLD),
	add_chapter(trna_clones, BOOK1, BOOK2),
	insert_entry(ac_on_trna,  (AC, TRNA), BOOK2, BOOK3),
	insert_entry(trna_clones, (AC, TRNA), BOOK3, BOOK4).

is_clone(AC, TRNA, BOOK) :-
	get_chapter(trna_clones, BOOK, CHAPTER),
	member((AC, TRNA), CHAPTER).

mutate_trna_anticodon(AC1, AC2, TRNA, BOOK1, BOOK2) :-
       replace_entry(ac_on_trna, (AC1, TRNA), (AC2, TRNA), BOOK1, BOOK2).

deactivate_synth(SYNTH, BOOK1, BOOK2) :-
	get_chapter(recognizes, BOOK1, CHAPTER1),
	once(member((SYNTH, TRNA_FAMILY, _), CHAPTER1)),
	replace((SYNTH, TRNA_FAMILY, _), (SYNTH, TRNA_FAMILY, _), CHAPTER1, CHAPTER2),
	list_to_set(CHAPTER2, CHAPTER3),
	replace_chapter(recognizes, BOOK1, CHAPTER3, BOOK2), !.
deactivate_synth(_, BOOK1, BOOK1).

% modify the genetic code
% ------------------------------------------------------------------------------

% code an amino acid (AA) with a new codon (C), i.e.
% - release the codon (C) by removing the corresponding tRNA
% - create a clone (CLONE) of the tRNA that is loaded with the target amino
%   acid (AA)
% - alter the anticodon of the cloned tRNA (CLONE) to match the target codon (C)
% - check if the code is still valid
% requirement:
% - the Manhatten distance between the old and new anticodon on the clone
%   should be three
% - the Manhattan distance between the new codon and each cognate codon for
%   AA should be three

recode_once(BOOK_INIT, BOOK1, BOOK3, RECIPE) :-
	RECIPE = [AC, AC_OLD, TRNA, TRNA_OLD, AA, AA_OLD],
	query([aa, ac, trna], [AA_OLD, AC, TRNA_OLD], BOOK1),
	% delete the old TRNA
	delete_trna(AC, TRNA_OLD, BOOK1, BOOK2),
	% get a new TRNA associated with AA and its current anticodon
	query([aa, ac, trna], [AA, AC_OLD, TRNA], BOOK_INIT),
	% assure that we selected another TRNA
	TRNA \= TRNA_OLD,
	% the new anticodon should have a different nucleotide
	% at every position
	pairwise_distinct(AC_OLD, AC),
	% there should be one new codon for aa at distance >2 from all cognate codons
	% clone the TRNA with new anticodon AC
	insert_cloned_trna(AC, TRNA, BOOK2, BOOK3).

recode_rec(0, BOOK_INIT, BOOK, BOOK, _, []) :- !,
	% for all clones (AC, TRNA) in the new codebook BOOK
	forall(is_clone(AC, TRNA, BOOK),
	       (% get the amino acid
		once(query([aa, ac, c, trna], [AA, AC, _, TRNA], BOOK)),
		% and the cognate codons from BOOK_INIT
		setof(C, query([aa, c], [AA, C], BOOK_INIT), COG_LIST),
		once((% find a codon C for AA in the new codebook BOOK
		      query([aa, ac, c, trna], [AA, AC, C, TRNA], BOOK),
		      % which is different to all cognate codons in COG_LIST
		      forall(member(C_COG, COG_LIST),
			     pairwise_distinct_2(C, C_COG)))))),
	codebook_is_complete(BOOK),
	codebook_is_unambiguous(BOOK).

recode_rec(N, BOOK_INIT, BOOK1, BOOK3, TRNA_LIST, RECIPE) :-
	M is N - 1,
	% select a TRNA including its AC for deletion
	append(_, [(AC, TRNA_OLD) | TRNA_TAIL], TRNA_LIST),
	%select((AC, TRNA_OLD), TRNA_LIST, TRNA_TAIL),
	% and fix it in the recipe
	RECIPE_HEAD = [AC,_,_,TRNA_OLD,_,_],
	% generate a new code
	recode_once(BOOK_INIT, BOOK1, BOOK2, RECIPE_HEAD),
	% update recipe list
	RECIPE = [RECIPE_HEAD | RECIPE_TAIL],
	% tail-recursive call	
	recode_rec(M, BOOK_INIT, BOOK2, BOOK3, TRNA_TAIL, RECIPE_TAIL).

recode(N, BOOK_INIT, BOOK_START, BOOK_NEW, RECIPE) :-
	setof((AC, TRNA), query([ac, trna], [AC, TRNA], BOOK_START), TRNA_LIST),
	recode_rec(N, BOOK_INIT, BOOK_START, BOOK_NEW, TRNA_LIST, RECIPE).

% accept options
% - deactivate(SYNTH) : deactivate the synthetase SYNTH before recoding

recode(N, BOOK_INIT, BOOK_START, BOOK_NEW, RECIPE, []) :- !,
       recode(N, BOOK_INIT, BOOK_START, BOOK_NEW, RECIPE).

recode(N, BOOK_INIT, BOOK_START, BOOK_NEW, RECIPE, [deactivate(SYNTH) | TAIL]) :-
	deactivate_synth(SYNTH, BOOK_START, BOOK_TMP),
	recode(N, BOOK_INIT, BOOK_TMP, BOOK_NEW, RECIPE, TAIL).

% ?- create_codebook(ecoli, B), recode(2, B, RECIPE).

% debugging:
% ?- create_codebook(ecoli, B), gtrace, recode(2, B, RECIPE).

% (CODEBOOK, RECIPE) -> NEW_CODEBOOK
% ------------------------------------------------------------------------------

implement_recipe_rec(        _, BOOK1, BOOK1,     []) :- !.

implement_recipe_rec(BOOK_ORIG, BOOK1, BOOK3, RECIPE) :-
	RECIPE = [RECIPE_HEAD | RECIPE_TAIL],
	recode_once(BOOK_ORIG, BOOK1, BOOK2, RECIPE_HEAD),
	implement_recipe_rec(BOOK_ORIG, BOOK2, BOOK3, RECIPE_TAIL).

implement_recipe(BOOK1, BOOK3, RECIPE) :-
	implement_recipe_rec(BOOK1, BOOK1, BOOK3, RECIPE).

% generate all possible codebooks/codes
% ------------------------------------------------------------------------------

generate_codebook(N, BOOK_INIT, BOOK_START, BOOK_NEW, OPTIONS) :-
	recode(N, BOOK_INIT, BOOK_START, BOOK_NEW, _, OPTIONS).
generate_codebook(N, BOOK_INIT, BOOK_START, BOOK_NEW) :-
	generate_codebook(N, BOOK_INIT, BOOK_START, BOOK_NEW, []).

generate_codebooks(N, BOOK_INIT, BOOK_START, BOOK_LIST, OPTIONS) :-
	setof(B, generate_codebook(N, BOOK_INIT, BOOK_START, B, OPTIONS), BOOK_LIST).
generate_codebooks(N, BOOK_INIT, BOOK_START, BOOK_LIST) :-
	generate_codebooks(N, BOOK_INIT, BOOK_START, BOOK_LIST, []).

generate_code(N, BOOK_INIT, BOOK_START, CODE, OPTIONS) :-
	recode(N, BOOK_INIT, BOOK_START, BOOK_NEW, _, OPTIONS),
	implements_code(BOOK_NEW, CODE).
generate_code(N, BOOK_INIT, BOOK_START, CODE) :-
	generate_code(N, BOOK_INIT, BOOK_START, CODE, []).

generate_codes(N, BOOK_INIT, BOOK_START, CODE_LIST, OPTIONS) :-
	setof(C, generate_code(N, BOOK_INIT, BOOK_START, C, OPTIONS), CODE_LIST).
generate_codes(N, BOOK_INIT, BOOK_START, CODE_LIST) :-
	generate_codes(N, BOOK_INIT, BOOK_START, CODE_LIST, []).

% multithreaded recoding
% ------------------------------------------------------------------------------

generate_codes_mt(N, BOOK_INIT, BOOK_START, CODE_LIST) :-
	once(generate_codebooks(1, BOOK_INIT, BOOK_START, BOOKS1)),
	concurrent_maplist(generate_codes(N-1,BOOK_INIT),BOOKS1,CODE_LLIST),
	foldl(union,CODE_LLIST,[],CODE_LIST).
