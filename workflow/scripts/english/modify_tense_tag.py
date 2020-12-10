import re

tam_re = re.compile(r'(.*)\((.*)\.(.*)\.(.*)\)')
    
def split_TAM(TAM_tag):
    """Split TAM tag and return parts as dict"""
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
