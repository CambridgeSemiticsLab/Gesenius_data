# use Spacy Matcher rules identify English tense constructions
# overlapping results will be filtered out and the longest matching span
# will be kept in its place

# these patterns can be inserted between verb auxiliaries
# and their heads to represent any number of interrupting
# adverbial modifiers:
advb_pronouns = {'TAG': {'IN':['RB', 'PRP']}, 'OP': '*'}
advbs = {'TAG': {'IN':['RB']}, 'OP': '*'}
non_verbs = {'TAG': {'NOT_IN':['VB', 'VBD', 'VBG', 'VBN', 'VBP', 'VBZ']}, 'OP': '*'}
pres_modals = ['let', 'may', 'shall', 'must', 'can', 'should']
past_modals = ['would', 'could']

rules = [

# -- present tense --
    (
        'PRES', 
        [
            {'TAG':{'IN':['VBZ', 'VBP']}, 'DEP': {'NOT_IN': ['aux']}},
        ]
    ),
    (
        'PRES does not',
        [
            {'LOWER': 'does'},
            {'LOWER': 'not'},
            non_verbs,
            {'TAG': 'VB'},
        ]
    ),
    (
        'PRES PROG', 
        [
            {'TAG': {'IN':['VBZ', 'VBP']}, 'LEMMA':'be'},
            non_verbs,
            {'TAG':'VBG', 'LEMMA': {'NOT_IN':['go']}},
        ]
    ),
    (
        'PRES PERF',
        [
            {'TAG': {'IN': ['VBZ', 'VBP']}, 'LEMMA': {'REGEX': 'have'}},
            non_verbs,
            {'TAG': 'VBN', 'DEP': {'NOT_IN': ['aux']}},
        ]
    ),
    (
        'PRES PERF PROG',
        [
            {'TAG': {'IN': ['VBZ', 'VBP']}, 'LEMMA': 'have'},
            non_verbs,
            {'TAG': 'VBN', 'LEMMA': 'be'},
            {'TAG': 'VBG'},
        ]
    ),

    # -- past --
    (
        'PAST',
        [
            {'TAG': 'VBD', 'DEP': {'NOT_IN':['aux']}},
        ]
    ),
    (
        'PAST did not', 
        [
            {'LOWER': 'did'},
            {'LOWER': 'not'},
            non_verbs,
            {'TAG': {'IN': ['VB', 'VBP']}},
        ]
    ),
    (
        'PAST PERF',
        [
            {'TAG': {'IN': ['VBD']}, 'LEMMA': 'have'},
            non_verbs,
            {'TAG': 'VBN', 'DEP': {'NOT_IN': ['aux']}},
        ]
    ),
    (
        'PAST PERF PROG',
        [
            {'TAG': {'IN': ['VBD']}, 'LEMMA': 'have'},
            non_verbs,
            {'TAG': 'VBN', 'LEMMA': 'be'},
            {'TAG': 'VBG'},
        ]
    ),
    (
        'PAST PROG',
        [
            {'TAG':'VBD', 'LEMMA': {'IN': ['be', 'keep']}},
            non_verbs,
            {'TAG': 'VBG'},
        ]
    ),

    # -- future --
    (
        'FUT',
        [
            {'TAG': 'MD', 'LEMMA': {'REGEX':'[wW]ill'}, 'DEP': 'aux'},
            non_verbs,
            {'TAG': 'VB', 'DEP': {'NOT_IN': ['aux']}},
        ]
    ),
    (
        'FUT going-to',
        [
            {'TAG': {'IN':['VBZ', 'VBP']}, 'LEMMA':'be'}, 
            non_verbs,
            {'LOWER': 'going'},
            {'LOWER': 'to'},
            {'TAG': 'VB'},
        ]
    ),
    (
        'FUT PERF',
        [
            {'TAG': 'MD', 'LEMMA': 'will'},
            non_verbs,
            {'TAG': {'IN': ['VB']}, 'LEMMA': 'have'},
            non_verbs,
            {'TAG': 'VBN', 'DEP': {'NOT_IN': ['aux']}},
        ]
    ),
    (
        'FUT PERF PROG',
        [
            {'TAG': 'MD', 'LEMMA': 'will'},
            non_verbs,
            {'TAG': {'IN': ['VB']}, 'LEMMA': 'have'},
            non_verbs,
            {'TAG': 'VBN', 'LEMMA': 'be'},
            {'TAG': 'VBG'},
        ]
    ),

    # -- modals --
    (
        'MOD',
        [
            {'TAG': {'IN':['VB', 'MD']}, 'LEMMA': {'IN': pres_modals}},
            non_verbs,
            {'TAG': 'VB'},
        ]
    ),
    (
        'MOD past',
        [
            {'TAG': {'IN':['VB', 'MD']}, 'LEMMA': {'IN': past_modals}},
            non_verbs,
            {'TAG': 'VB'},
        ]
    ),
     (
        'MOD lest',
        [
            {'lower': 'lest'},
            non_verbs,
            {'TAG': {'IN': ['VBP', 'VB']}},
        ],
    ),

    # -- habituals --
    (
        'FUT-IN-PAST',
        [
            {'LOWER': 'would', 'DEP': {'IN': ['aux']}},
            non_verbs,
            {'TAG':'VB'}
        ]
    ),
    (
        'HAB used to',
        [
            {'TAG': 'VBD', 'lemma': 'use'},
            {'LOWER': 'to'},
            {'TAG': 'VB'},
        ],
    ),

    # -- emphatics --
    (
        'PRES do-support',
        [
            {'LOWER': 'does'},
            {'TAG': 'VB'},
        ],
    ),
    (
        'PAST do-support',
        [
            {'LOWER': 'did'},
            {'TAG': 'VB'},
        ],
    ),

    # -- participles and infinitives --
    (
        'PRES PART',
        [
            {'TAG': 'VBG', 'DEP': {'NOT_IN':['aux']}},
        ]
    ),
    (
        'PRES PART',
        [
            {'TAG': {'IN': ['JJ']}, 'LOWER':{'REGEX':'^.+ing$'}},
        ]
    ),
    (
        'TO INF',
        [
            {'LOWER': 'to'},
            {'TAG':'VB'}
        ]
    ),
   
    # -- imperatives --
    # the imperative in English consists of the base
    # form of a verb and is distinguished syntactically.
    # This makes it a complex tense to automatically recognize.
    # Generally an imperative verb is a base form that occurs
    # at the front of a sentence or clause; adverbial elements
    # may instead precede the verb
    (
        'VB',
        [ 
            {'TAG': 'VB'}
        ]
    ),
    (
        'IMPV',
        [
            {'TAG': 'VB', 'DEP':{'NOT_IN':['aux']}, 'IS_SENT_START': True},
        ]
    ),
    (
        'IMPV',
        [
            {'TAG': 'VB', 'DEP':{'NOT_IN':['aux']}, 'IS_TITLE': True},
        ]
    ),
    (
        'IMPV',
        [
            {'TAG': 'RB', 'IS_SENT_START': True},
            {'TAG': 'VB', 'DEP':{'NOT_IN':['aux']}} 
        ]
    ),
    (
        'IMPV',
        [
            {'TAG': 'RB', 'IS_TITLE': True},
            {'TAG': 'VB', 'DEP':{'NOT_IN':['aux']}} 
        ]
    ),
    (
        'IMPV',
        [
            {'IS_PUNCT': True},
            {'TAG': 'VB', 'DEP':{'NOT_IN':['aux']}, 'IS_SENT_START': False, 'IS_TITLE': False} 
        ]
    ),
    (
        '?IMPV',
        [
            {'IS_PUNCT': True},
            {'TAG': 'VBP', 'DEP':{'NOT_IN':['aux']}, 'IS_SENT_START': False, 'IS_TITLE': False} 
        ]
    ),

    (
        'IMPV',
        [
            {'TAG': 'VBP', 'IS_SENT_START': True, 'LEMMA': {'NOT_IN': ['have', 'be', 'do']}},
        ],
    ),
    (
        'IMPV',
        [
            {'TAG': 'VBP', 'IS_TITLE': True, 'LEMMA': {'NOT_IN': ['have', 'be', 'do']}},
        ],
    ),
    (
        '?IMPV',
        [
            {'TAG': 'VBP', 'IS_TITLE': True, 'LOWER': {'IN': ['have', 'do', 'be']}},
            {'LOWER': {'NOT_IN': ['you', 'i', 'not', 'any', 'men', ]}},
        ],
    ),

    (
        'IMPV not',
        [
            {'LOWER': 'do'},
            {'LEMMA': 'not'},
            non_verbs,
            {'TAG': 'VB'}
        ]
    ),
]
