"""
Matcher object rules for Spacy (version 2.x)
for matching tenses of interest in French
"""

rules = [ 

    # -- passé simple --
    (   
        'passé_simp',
        [               
            {'TAG': {'REGEX': 'VERB__Mood=Ind.*Tense=Past\|VerbForm=Fin'}},
        ]   
    ),
    
    # -- passé composé --
    (   
        'passé_comp',        
        [   
            {'TAG': {'REGEX': 'AUX__Mood=Ind.*Tense=Pres\|VerbForm=Fin'}, 'LEMMA': {'IN': ['avoir', 'être']}},
            {'TAG': {'REGEX': 'VERB.*Tense=Past\|VerbForm=Part'}},
        ]   
    ),
    
    # -- l'imparfait --
    (   
        'imparfait',
        [               
            {'TAG': {'REGEX': 'VERB__Mood=Ind.*Tense=Imp\|VerbForm=Fin'}},
        ],   
    ),
#    (
#        'imparfait',
#        [               
#            {'TAG': {'REGEX': 'AUX__Mood=Ind.*Tense=Imp\|VerbForm=Fin'}},
#        ]   
#    ), 
]
