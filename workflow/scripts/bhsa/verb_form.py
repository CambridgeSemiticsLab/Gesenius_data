from positions import PositionsTF

def get_verbform(node, api):
    """Remap BHSA verb tense values to custom values
    
    Args:
        node: int, representing BHSA node id
        api: an instance of Text-Fabric with BHSA loaded
    """
    
    tense_map = {
            'impf': 'yqtl',
            'perf': 'qtl',
            'ptca': 'ptcp',
    }
    
    bhsa_tense = api.F.vt.v(node)
    verb_form = tense_map.get(bhsa_tense, bhsa_tense)
    P = PositionsTF(node, 'clause', api)
    
    # adjust weqatal
    if verb_form == 'qtl' and P.get(-1, 'lex') == 'W':
        return 'wqtl'
        
    return verb_form
