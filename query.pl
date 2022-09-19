% -*- mode: prolog; -*-

:- module(query,
	  [ create_codebook/2
	  , query/3 ]).

:- use_module(codebook).

% Glossary
% ------------------------------------------------------------------------------
% - C     codon
% - AC    anticodon (before modifications)
% - ACM   anticodon (after  modifications)
% - AA    amino acid

% database utilities
% ------------------------------------------------------------------------------

trna_name(TRNA, AA, ID) :-
	atomic_list_concat([AA, tRNA, ID], '_', TRNA).

trna_family_name(TRNA_FAMILY, AA) :-
	atomic_list_concat([AA, tRNA], '_', TRNA_FAMILY).

% database access
% ------------------------------------------------------------------------------

codon(C, BOOK) :-
	get_chapter(codons, BOOK, CHAPTER),
	member(C, CHAPTER).

stop_codon(C, BOOK) :-
	get_chapter(stop_codons, BOOK, CHAPTER),
	member(C, CHAPTER).

amino_acid(AA, BOOK) :-
	get_chapter(amino_acids, BOOK, CHAPTER),
	member(AA, CHAPTER).

watson_crick_pair(C1, C2, BOOK) :-
	get_chapter(watson_crick_pairs, BOOK, CHAPTER),
	member((C1, C2), CHAPTER).

wobble_pair(C1, C2, BOOK) :-
	get_chapter(wobble_pairs, BOOK, CHAPTER),
	member((C1, C2), CHAPTER).

ac_on_trna(AC, TRNA, BOOK) :-
	get_chapter(ac_on_trna, BOOK, CHAPTER),
	member((AC, TRNA), CHAPTER).

aa_is_loaded_by(AA, SYNTH, BOOK) :-
	get_chapter(aa_is_loaded_by, BOOK, CHAPTER),
	member((AA, SYNTH), CHAPTER).

modifies_ac(ENZYME, TRNA, AC, ACM, BOOK) :-
	(nonvar( AC) -> atom_chars( AC,  ACL) ; true),
	(nonvar(ACM) -> atom_chars(ACM, ACML) ; true),
	get_chapter(modifies_ac, BOOK, CHAPTER),
	% variables in the chapter might get instantiated
	% by searching for a match; hence, the chapter
	% needs to be copied before
	copy_term(CHAPTER, CHAPTER_COPY),
	member((ENZYME, TRNA, ACL, ACML), CHAPTER_COPY),
	atom_chars( AC,  ACL),
	atom_chars(ACM, ACML).

matures(TRNA, AC, ACM, BOOK) :-
	nonvar(AC), nonvar(TRNA), !,
	(modifies_ac(_, TRNA, AC, ACM_C, BOOK)
	 -> ACM_C = ACM
          ; AC = ACM).
matures(TRNA, AC, ACM, BOOK) :-
	var(AC), nonvar(ACM), !,
	(modifies_ac(_, TRNA_C, AC_C, ACM, BOOK)
	 -> TRNA_C = TRNA, AC_C = AC
          ; AC = ACM).
matures(_, _, _, _) :-
	throw(error(instantiation_error, _)).

/** <examples>

 ?- create_codebook(ecoli, B), query:matures(ile_tRNA_X, cau, ACM, B).
 ACM = lau.
 ?- create_codebook(ecoli, B), query:matures(ile_tRNA_X, cau, cau, B).
 false.
 ?- create_codebook(ecoli, B), query:matures(ile_tRNA_Y, cau, cau, B).
 true.
 ?- create_codebook(ecoli, B), query:matures(foobar, acg, icg, B).
 true.
 ?- create_codebook(ecoli, B), query:matures(foobar, acg, igg, B).
 false.
 ?- create_codebook(ecoli, B), query:matures(foobar, acg, ACM, B).
 ACM = icg.
 ?- create_codebook(ecoli, B), query:matures(TRNA, AC, lau, B).
 TRNA = ile_tRNA_X,
 AC = cau.
 ?- create_codebook(ecoli, B), query:matures(TRNA, AC, cgg, B).
 AC = cgg.
 ?- create_codebook(ecoli, B), query:matures(TRNA, AC, ACM, B).
 ERROR: Arguments are not sufficiently instantiated
 ?- create_codebook(ecoli, B), query:matures(TRNA, acg, ACM, B).
 ERROR: Arguments are not sufficiently instantiated

*/

% Example:
% initial_codebook(B), modifies_ac(_, auu, X, B), modifies_ac(_, agc, Y, B).

recognizes(SYNTH, TRNA_FAMILY, AC, BOOK) :-
	get_chapter(recognizes, BOOK, CHAPTER),
	copy_term(CHAPTER, CHAPTER_COPY),
	member((SYNTH, TRNA_FAMILY, AC), CHAPTER_COPY).

% database queries
% ------------------------------------------------------------------------------

pairs_with(C, ACM, BOOK) :-
	(nonvar(  C) -> atom_chars(  C, [  C1,  C2,  C3]) ; true),
	(nonvar(ACM) -> atom_chars(ACM, [ACM1,ACM2,ACM3]) ; true),
	% check pairing
	watson_crick_pair(C1, ACM3, BOOK),
	watson_crick_pair(C2, ACM2, BOOK),
	wobble_pair(C3, ACM1, BOOK),
	% C should be a valid codon
	codon(C, BOOK),
	% mRNA 5' -- C1 C2 C3 -- 3'
	atom_chars(  C, [  C1,  C2,  C3]),
	% tRNA 5' -- A1 A2 A3 -- 3'
	atom_chars(ACM, [ACM1,ACM2,ACM3]).

trna_family(TRNA_FAMILY, TRNA, BOOK) :-
	(nonvar(TRNA)
	 -> trna_name(TRNA, AA, _)
	  ; true),
	get_chapter(ac_on_trna, BOOK, CHAPTER),
	member((_, TRNA), CHAPTER),
	trna_name(TRNA, AA, _),
	trna_family_name(TRNA_FAMILY, AA).

aa_is_loaded_on(AA, TRNA, AC, BOOK) :-
	aa_is_loaded_by(AA, SYNTH, BOOK),
	trna_family(TRNA_FAMILY, TRNA, BOOK),
	recognizes(SYNTH, TRNA_FAMILY, AC, BOOK).

is_valid_code(C, AC, TRNA, AA, BOOK) :-
	% AC: anticodon before modifications
	% ACM: modified anticodon
	ac_on_trna(AC, TRNA, BOOK),
	(modifies_ac(_, TRNA, AC, ACM, BOOK) -> true; AC = ACM),
	% this check needs only be executed once
	pairs_with(C, ACM, BOOK),
	not(stop_codon(C, BOOK)),
	aa_is_loaded_on(AA, TRNA, ACM, BOOK).

% Examples:
% ?- is_valid_code(gca, ugc, ala_tRNA_T, ala).
% ?- is_valid_code(C, AC, TRNA, AA).

% general query interface
% ------------------------------------------------------------------------------

query([aa], [AA], BOOK) :- !,
	amino_acid(AA, BOOK).

query([c], [C], BOOK) :- !,
	codon(C, BOOK).

query([stop_c], [C], BOOK) :- !,
	stop_codon(C, BOOK).

query([trna], [TRNA], BOOK) :- !,
	ac_on_trna(_, TRNA, BOOK).

query([ac], [AC], BOOK) :- !,
	ac_on_trna(AC, _, BOOK).

query([ac, trna], [AC, TRNA], BOOK) :- !,
	ac_on_trna(AC, TRNA, BOOK).

query([aa, ac, trna], [AA, AC, TRNA], BOOK) :- !,
	ac_on_trna(AC, TRNA, BOOK),
	matures(TRNA, AC, ACM, BOOK),
	aa_is_loaded_on(AA, TRNA, ACM, BOOK).

query([trna, trna_family], [TRNA, TRNA_FAMILY], BOOK) :- !,
	trna_family(TRNA_FAMILY, TRNA, BOOK).

query([aa, ac, c, trna], [AA, AC, C, TRNA], BOOK) :- !,
	is_valid_code(C, AC, TRNA, AA, BOOK).

query([aa, c], [AA, C], BOOK) :- !,
	setof((AA, C), AC^TRNA^(query([aa, ac, c, trna], [AA, AC, C, TRNA], BOOK)),
	      LIST),
	member((AA, C), LIST).

query([acm, c], [ACM, C], BOOK) :- !,
	pairs_with(C, ACM, BOOK).

query([ac, c], [AC, C], BOOK) :- !,
	pairs_with(C, ACM, BOOK),
	\+ stop_codon(C, BOOK),
	ac_on_trna(AC, TRNA, BOOK),
	matures(TRNA, AC, ACM, BOOK).

query([ac, acm, trna], [AC, ACM, TRNA], BOOK) :- !,
	ac_on_trna(AC, TRNA, BOOK),
	matures(TRNA, AC, ACM, BOOK).

query([ac, acm], [AC, ACM], BOOK) :- !,
	ac_on_trna(AC, TRNA, BOOK),
	matures(TRNA, AC, ACM, BOOK).

query([aa, synth], [AA, SYNTH], BOOK) :- !,
	aa_is_loaded_by(AA, SYNTH, BOOK).

query([ac, synth, trna_family], [AC, SYNTH, TRNA_FAMILY], BOOK) :- !,
	recognizes(SYNTH, TRNA_FAMILY, AC, BOOK).

query([synth, trna], [SYNTH, TRNA], BOOK) :- !,
	recognizes(SYNTH, TRNA_FAMILY, _, BOOK),
	trna_family(TRNA_FAMILY, TRNA, BOOK).
