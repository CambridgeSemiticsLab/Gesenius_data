"""
Tag clause arguments
"""

def clause_objects(verb, clause_atom, clause, api):
    """Search a given clause for any marked objects."""
    
    F, E, L = api.F, api.E, api.L

    clause_phrases = L.d(clause, 'phrase')
    daughters = E.mother.t(clause_atom)
    daught_relas = set(F.rela.v(d) for d in daughters)
    daught_codes = set(F.code.v(d) for d in daughters)
    phrase_functs = set(F.function.v(p) for p in clause_phrases)
    obj_phrases = list(p for p in clause_phrases
                            if F.function.v(p) == 'Objc')
    
    # we can count direct speech clauses following amar as direct objects
    #amar_obj = F.lex.v(verb) == '>MR[' and 999 in daught_codes 
    
    data = {}

    # evaluate conditions for object presence
    data['has_objc'] = any([
        bool(obj_phrases),
        'PreO' in phrase_functs,
        'PtcO' in phrase_functs,
        'Objc' in daught_relas,
        #amar_obj,
    ])

    # retrieve object positionality for phrasal objects
    if data['has_objc']:
        if obj_phrases:
            pred_phrase = L.u(verb, 'phrase')[0]
            phrases = sorted([pred_phrase] + obj_phrases)

            # build order label with allowance for double objects
            order_label = ''
            for ph in phrases:
                if ph == pred_phrase:
                    order_label += 'V'
                else:
                    order_label += 'O'
            data['objc_pos'] = order_label
        else: 
            data['objc_pos'] = 'VO?'
    else:
        data['objc_pos'] = ''

    return data


def get_loca_assocs(api):
    """Identifies lexemes statistically associated with location."""
    
    F, E, L = api.F, api.E, api.L
    lexs = set()
    for cmpl in F.function.s('Loca'):
        heads = E.nhead.t(cmpl)
        for head in heads:
            loc_assoc = F.funct_assoc.v(head) or 0
            if loc_assoc >= 1.5:
                lexs.add(F.lex.v(head))
    return lexs

def clause_locas(verb, loca_lexs, api):
    """Tag location data within a clause"""

    F, E, L = api.F, api.E, api.L

    data = {
        'has_loca': False,
        'loca_type': [],
        'loca_heads': [],
    }

    clause = L.u(verb, 'clause')[0]
    clause_phrases = L.d(clause, 'phrase')
    
    # check for location
    loca_ph = [p for p in clause_phrases if F.function.v(p) == 'Loca']
    
    # attempt to find a locative complement
    if not loca_ph:
        cmpl_phrs = [p for p in clause_phrases if F.function.v(p) == 'Cmpl']
        for cp in cmpl_phrs:
            for head in E.nhead.t(cp):
                if F.lex.v(head) in loca_lexs:
                    loca_ph.append(cp)

    # detect presence of location
    if loca_ph:
        data['has_loca'] = bool(loca_ph)

    # tag various features on the location phrases
    for lp in loca_ph:
        data['loca_type'].append(F.typ.v(lp))
        for head in E.nhead.t(lp):
            data['loca_heads'].append(F.lex.v(head))
        
    
    data['loca_type'] = '|'.join(data['loca_type'])
    data['loca_heads'] = '|'.join(data['loca_heads'])

    return data

