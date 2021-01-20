import re

tam_re = re.compile(r'(.*)\((.*)\.(.*)\.(.*)\)')

def split_TAM(TAM_tag):
    """Split TAM tag and return parts as dict
    
    ! DEPRECATED !
    Starting with yiqtol, we have deprecated
    this function since tense tags have now 
    changed. The function process_TAM now
    fulfills the old role.
    """
    tam_match = tam_re.match(TAM_tag)
    if tam_match:
        name, tense, aspect, modality = tam_match.groups()
        return {
            'tense': tense or '',
            'aspect': aspect or '',
            'modality': modality or '',
            'TAM': f'{tense}.{aspect}.{modality}',
            'TAMtag': name.strip(),
        }
    else:
        return {
            'tense': '',
            'aspect': '',
            'modality': '',
            'TAM': '',
            'TAMtag': '',
        }

def process_TAM(full_tam):
    """Take in a TAM tag and build TAM features."""

    tam_simp_map = {
        'PRES do not': 'PRES',
        'PRES question': 'PRES',
        'PRES wh-question': 'PRES',
        'PRES PERF question': 'PRES PERF',
        'PRES PERF wh-question': 'PRES PERF',
        'PAST did not': 'PAST',
        'PAST question': 'PAST',
        'PAST wh-question': 'PAST',
        'PAST PROG keep': 'PAST PROG',
        'FUT question': 'FUT',
        'FUT wh-question': 'FUT',
        'PRES do-support': 'PRES',
        'PAST do-support': 'PAST',
        'IMPV do not': 'IMPV',
        'MOD quest shall': 'MOD shall',
        'MOD quest should': 'MOD should',
        'MOD quest would': 'MOD would',
        'MOD quest can': 'MOD can',
        'MOD quest could': 'MOD could',
    }

    features = {
        'TAM': full_tam,
        'TAMsimp': tam_simp_map.get(full_tam, full_tam)
    }

    return features
