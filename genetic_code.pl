% -*- mode: prolog; -*-

:- module(genetic_code,
	  [ implements_code/2
	  , code_is_ambiguous/1
	  , code_is_unambiguous/1
	  , codebook_is_complete/1
	  , codebook_is_ambiguous/1
	  , codebook_is_unambiguous/1
	  ]).

:- use_module(query).

% convert a codebook to its genetic code
% ------------------------------------------------------------------------------

aa_codons(AA, C_LIST2, BOOK) :-
        setof(C, query([aa, c], [AA, C], BOOK), C_LIST1),
	sort(C_LIST1, C_LIST2).

implements_code(BOOK, CODE2) :-
	setof((AA, C_LIST),
	      aa_codons(AA, C_LIST, BOOK),
	      CODE1),
	sort(CODE1, CODE2).

% code properties
% ------------------------------------------------------------------------------

code_is_ambiguous(CODE) :-
	member((AA1, CODONS1), CODE),
	member(C, CODONS1),
	member((AA2, CODONS2), CODE),
	member(C, CODONS2),
	AA1 \= AA2.

code_is_unambiguous(CODE) :-
	\+ code_is_ambiguous(CODE).

% codebook properties
% ------------------------------------------------------------------------------

% true if there exists at least one valid codon for every
% amino acid
codebook_is_complete_aa(BOOK) :-
        forall(query([aa], [AA], BOOK),
	       query([aa, ac, c, trna], [AA, _, _, _], BOOK)).
% true if every codon codes for some amino-acid
codebook_is_complete_c(BOOK) :-
	forall(query([c], [C], BOOK),
	       query([aa, ac, c, trna], [_, _, C, _], BOOK)).

codebook_is_complete(BOOK) :-
        codebook_is_complete_aa(BOOK),
	codebook_is_complete_c(BOOK).

codebook_is_ambiguous(BOOK) :-
        implements_code(BOOK, CODE),
	code_is_ambiguous(CODE).

codebook_is_unambiguous(BOOK) :-
        implements_code(BOOK, CODE),
	code_is_unambiguous(CODE).
