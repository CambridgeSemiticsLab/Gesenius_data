def is_imperative(token, debug=False):
    """Test whether a verb is an imperative."""

    ancestors = list(token.ancestors)
    subtree = list(token.subtree)
    ancest_tags = set(t.tag_ for t in ancestors if t.tag_.startswith('VB'))
    ancest_lemmas = set(t.lemma_ for t in ancestors)
    subtree_deps = set(t.dep_ for t in subtree)
    head_deps = set(t.dep_ for t in token.head.children)
    ancest_lemmas = set(t.lemma_ for t in ancestors)
    clause_pos = subtree.index(token)
    modal_set = {
        'let', 'may', 
        'shall', 'must', 'can'
    }
    auxs = {'aux', 'auxpass'}
    
    ancestor_rules = [
        ancest_tags.issubset({'VB'}),
        token.dep_ != 'conj',
    ]
    if ancestors:
        ancestor_rules.append(is_imperative(token.head))
    position_rules = [
        clause_pos == 0,
        not ancestors,
    ]
    
    # check for boundary
    rules = [
        token.tag_ == 'VB',
        token.dep_ not in auxs,
        any(ancestor_rules),
        not ancest_lemmas & modal_set,
        not head_deps & auxs,
        any(position_rules),
        token.lemma_ not in modal_set,
    ]
    # specify that word must occur at 
    # a major boundary
    if token.i != 0:
        pre_t = token.doc[token.i-1]
        punct_rules = any([
            pre_t.is_punct,
            pre_t.lemma_ == 'and',
            (pre_t.is_title and pre_t.tag_ == 'RB'),
            pre_t.lemma == 'please',
        ])
        rules.append(punct_rules)
    else:
        rules.append(token.is_title)
    
    
    if debug:
        print(rules)
        print('anc rules:', ancestor_rules)
        print('pos rules:', position_rules)
    
    if all(rules):
        return True
    else:
        return False

def test_span_impv(span):
    """Test whether span truly contains an imperative"""
    for token in span:
        if is_imperative(token):
            return True
    return False
