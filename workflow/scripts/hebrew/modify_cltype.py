# modify default clause type values for BHSA

import re

cl_type_res = [
    ('Way.*', 'W'),
    ('^[xX].*', 'x'),
    ('^W[xX].*', 'Wx'),
    ('^W.*', 'W'),
    ('^_W_.+', 'Wx'),
    ('^_W_$', 'W'),
]

cl_type_res = [(re.compile(search), replace) for search, replace in cl_type_res]

def search_replace_re(string, patterns, default=None):
    """Match a regex string"""
    for search, replace in patterns:
        if search.match(string):
            return replace
    return default

def simplify_cl_type(clause_atom, prec_lexs, api):
    """Simplify a clausetype string into (x|X|Ø)Verb"""
    
    typ = api.F.typ.v(clause_atom)
    
    # apply to verbs missing X|x data
    if typ in {'MSyn', 'CPen', 'Voct', 'InfC', 'Ellp', 'Ptcp'}:
        return search_replace_re(prec_lexs, cl_type_res, 'Ø')
        
    # apply to other types
    return search_replace_re(typ, cl_type_res, 'Ø')
