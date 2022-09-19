% -*- mode: prolog; -*-

:- module(codebook,
	  [ create_codebook/2
	  , add_chapter/3
	  , get_chapter/3
	  , replace_chapter/4
	  , replace_entry/5
	  , delete_entry/4
	  , insert_entry/4 ]).

:- use_module(listutil).

% construct a codebook
% ------------------------------------------------------------------------------

chapter_list(NAMESPACE, CHAPTER_LIST) :-
	findall(CHAPTER, NAMESPACE:codebook(CHAPTER, _), CHAPTER_LIST).

create_codebook(NAMESPACE, BOOK) :-
	chapter_list(NAMESPACE, CHAPTER_LIST),
	bagof((CHAPTER_NAME, CHAPTER),
	      (member(CHAPTER_NAME, CHAPTER_LIST),
	       NAMESPACE:codebook(CHAPTER_NAME, CHAPTER)), BOOK).

% database utilities
% ------------------------------------------------------------------------------

add_chapter(WHICH_CHAPTER, BOOK1, BOOK2) :-
	\+ member((WHICH_CHAPTER, _), BOOK1)
	-> BOOK2 = [(WHICH_CHAPTER, []) | BOOK1]
	 ; BOOK2 = BOOK1.

get_chapter(WHICH_CHAPTER, BOOK, CHAPTER) :-
	member((WHICH_CHAPTER, CHAPTER), BOOK), !.

replace_chapter(_, [], _, []) :- !.
replace_chapter(WHICH_CHAPTER, BOOK, NEW_CHAPTER, NEW_BOOK) :-
	    BOOK = [(WHICH_CHAPTER,           _)|T],
	NEW_BOOK = [(WHICH_CHAPTER, NEW_CHAPTER)|T],
	!.
replace_chapter(WHICH_CHAPTER, BOOK, NEW_CHAPTER, NEW_BOOK) :-
	    BOOK = [H|T1],
	NEW_BOOK = [H|T2],
	replace_chapter(WHICH_CHAPTER, T1, NEW_CHAPTER, T2).

replace_entry(WHICH_CHAPTER, FROM, TO, BOOK, NEW_BOOK) :-
	get_chapter(WHICH_CHAPTER, BOOK, CHAPTER),
	replace(FROM, TO, CHAPTER, NEW_CHAPTER),
	replace_chapter(WHICH_CHAPTER, BOOK, NEW_CHAPTER, NEW_BOOK), !.

delete_entry(WHICH_CHAPTER, ENTRY, BOOK, NEW_BOOK) :-
	get_chapter(WHICH_CHAPTER, BOOK, CHAPTER),
	delete(ENTRY, CHAPTER, NEW_CHAPTER),
	replace_chapter(WHICH_CHAPTER, BOOK, NEW_CHAPTER, NEW_BOOK), !.

insert_entry(WHICH_CHAPTER, ENTRY, BOOK, NEW_BOOK) :-
	get_chapter(WHICH_CHAPTER, BOOK, CHAPTER),
	insert(ENTRY, CHAPTER, NEW_CHAPTER),
	replace_chapter(WHICH_CHAPTER, BOOK, NEW_CHAPTER, NEW_BOOK), !.

% codebook properties
% ------------------------------------------------------------------------------

