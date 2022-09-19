
% change working directory to examples so that files can be easily included
% ------------------------------------------------------------------------------

:- working_directory(CWD, CWD),
    concat(CWD, "examples/", NCWD),
    working_directory(CWD, NCWD).

% allow loading databases
% ------------------------------------------------------------------------------

:- multifile sandbox:safe_primitive/1.

:- use_module(ecoli).

sandbox:safe_primitive(ecoli:codebook(_, _)).
sandbox:safe_primitive(query:matures(_,_,_,_)).
sandbox:safe_primitive(query:pairs_with(_,_,_)).
sandbox:safe_primitive(query:stop_codon(_,_)).
sandbox:safe_primitive(xenoprint:xenoprint(_,_)).
sandbox:safe_primitive(xenoprint:xenoprint(_,_,_)).
