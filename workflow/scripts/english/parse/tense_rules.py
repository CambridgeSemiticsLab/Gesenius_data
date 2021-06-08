# use Spacy Matcher rules identify English tense constructions
# overlapping results will be filtered out and the longest matching span
# will be kept in its place

# these patterns can be inserted between verb auxiliaries
# and their heads to represent any number of interrupting
# adverbial modifiers:
advb_pronouns = {'TAG': {'IN':['RB', 'PRP']}, 'OP': '*'}
advbs = {'TAG': {'IN':['RB']}, 'OP': '*'}
non_verbs = {'TAG': {'NOT_IN':['MD', 'VB', 'VBD', 'VBG', 'VBN', 'VBP', 'VBZ']}, 'OP': '*'}
not_aux = {'NOT_IN': ['aux']}

ambiguous_pasts = [
    'put', 'set', 'cut', 'lay', 'cast', 
    'spread', 'spit', 'read', 'rid', 'darken',
]

rules = [
# -- present tense --
    (
        '?PRES',
        [
            {'TAG': 'VBP', 'DEP': not_aux},
        ]
    ),
    (
        'PRES',
        [
            {'LEMMA': '-PRON-'},
            {'TAG': 'VBP', 'DEP': not_aux, 'LEMMA': {'NOT_IN': ['be']}},
        ]
    ),
    (
        'PRES',
        [
            {'DEP': 'nsubj', 'LEMMA': {'NOT_IN': ['-PRON-']}},
            {'TAG': 'VBP', 'DEP': not_aux, 'LEMMA': {'NOT_IN': ['be']}, 'DEP': {'NOT_IN': ['aux', 'ccomp']}},
        ]
    ),
    (
        '?PRES',
        [
            {'DEP': 'nsubj', 'LEMMA': {'NOT_IN': ['-PRON-']}},
            {'TAG': 'VBP', 'DEP': not_aux, 'LEMMA': {'NOT_IN': ['be']}, 'DEP': {'IN': ['ccomp']}},
        ]
    ),

    (
        'PRES', 
        [
            {'TAG': {'IN': ['VBZ', 'VBP']}, 'DEP': not_aux, 'LEMMA': {'IN': ['be']}},
        ]
    ),
    (
        'PRES',
        [
            {'TAG': {'IN': ['VBZ']}, 'DEP': not_aux, 'LEMMA': {'NOT_IN': ['be']}},
        ]
    ),
    (
        'PRES do not',
        [
            {'DEP': 'nsubj'},
            advbs,
            {'TEXT': 'do', 'IS_SENT_START': False, 'DEP': 'aux'},
            {'LEMMA': 'not'},
            non_verbs,
            {'TAG': 'VB'},
        ]
    ),
    (
        'PRES do not',
        [
            {'TEXT': 'does', 'DEP': 'aux'},
            {'LEMMA': 'not'}, 
            non_verbs,
            {'TAG': 'VB'},
        ],
    ),
    (
        'PRES question',
        [
            {'TEXT': {'IN': ['Do', 'Does']}},
            {'LEMMA': '-PRON-'},
            non_verbs,
            {'TAG': {'IN': ['VBP', 'VB']}, 'DEP': not_aux},
        ],
    ),
    (
        'PRES wh-question',
        [
            {'TAG': 'WRB'},
            {'TEXT': {'IN': ['does', 'do']}},
            non_verbs,
            {'TAG': {'IN': ['VBP', 'VB']}, 'DEP': not_aux},
        ],
    ),
    (
        'PRES PROG', 
        [
            {'TEXT': {'IN': ['are', 'am', 'is']}},
            {'TAG':'VBG', 'LEMMA': {'NOT_IN':['go']}},
        ]
    ),
    (
        'PRES PERF',
        [
            {'TAG': {'IN': ['VBZ', 'VBP']}, 'LEMMA': {'REGEX': 'have'}},
            non_verbs,
            {'TAG': 'VBN', 'DEP': not_aux},
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
    (
        'PRES PERF question',
        [
            {'TEXT': {'IN': ['Have', 'Has']}},
            non_verbs,
            {'TAG': 'VBN', 'DEP': not_aux},
        ],
    ),
    (
        'PRES PERF wh-question',
        [
            {'TAG': 'WRB'},
            {'TEXT': {'IN': ['have', 'has']}},
            non_verbs,
            {'TAG': 'VBN', 'DEP': not_aux},
        ],
    ),

    # -- past --
    (
        'PAST',
        [
            {'TAG': 'VBD', 'DEP': {'NOT_IN':['aux']}, 'LEMMA': {'NOT_IN': ambiguous_pasts}},
        ]
    ),
    (
        '?PAST',
        [
            {'TAG': 'VBD', 'DEP': {'NOT_IN':['aux']}, 'LEMMA': 'put'},
        ]

    ),
    (
        'PAST did not', 
        [
            {'LOWER': 'did'},
            {'LEMMA': 'not'},
            non_verbs,
            {'TAG': {'IN': ['VB', 'VBP']}},
        ]
    ),
    (
        'PAST question',
        [
            {'TEXT': 'Did'},
            non_verbs,
            {'TAG': {'IN': ['VB']}},
        ]
    ),
    (
        'PAST wh-question',
        [
            {'TAG': 'WRB'},
            {'TEXT': 'did'},
            non_verbs,
            {'TAG': {'IN': ['VB']}},
        ]
    ),

    (
        'PAST PERF',
        [
            {'TAG': {'IN': ['VBD']}, 'LEMMA': 'have', 'IS_TITLE': False},
            non_verbs,
            {'TAG': 'VBN', 'DEP': not_aux},
        ]
    ),
    (
        'PAST PERF PROG',
        [
            {'TAG': {'IN': ['VBD']}, 'LEMMA': 'have'},
            {'TEXT': 'been'},
            {'TAG': 'VBG'},
        ]
    ),
    (
        'PAST PROG',
        [
            {'TEXT': {'IN': ['was', 'were']}},
            {'TAG': 'VBG'},
        ]
    ),
    (
        'PAST PROG keep',
        [
            {'TAG':'VBD', 'LEMMA': 'keep'},
            {'TAG': 'VBG'},
        ]
    ),

    # -- future --
    (
        'FUT',
        [
            {'LEMMA': 'will', 'DEP': 'aux', 'IS_TITLE': False, 'IS_SENT_START': False},
            non_verbs,
            {'TAG': 'VB', 'DEP': not_aux},
        ]
    ),
    (
        'FUT PROG',
        [
            {'LEMMA': 'will', 'DEP': 'aux'},
            {'LEMMA': 'be'},
            {'TAG':'VBG'},
        ],
    ),
    (
        'FUT question',
        [
            {'TEXT': 'Will'},
            {'LEMMA': '-PRON-'},
            non_verbs,
            {'TAG': 'VB'}
        ],
    ),
    (
        '?FUT question',
        [
            {'TEXT': 'will', 'IS_SENT_START': True},
            {'LEMMA': '-PRON-'},
            non_verbs,
            {'TAG': 'VB'}
        ],
    ),
    (
        'FUT wh-question',
        [
            {'TAG': 'WRB'},
            {'TEXT': 'will'},
            {'LEMMA': '-PRON-'},
            non_verbs,
            {'TAG': 'VB'}
        ],
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
            {'TAG': 'MD', 'LOWER': 'will'},
            non_verbs,
            {'TAG': {'IN': ['VB']}, 'LEMMA': 'have'},
            non_verbs,
            {'TAG': 'VBN', 'DEP': not_aux},
        ]
    ),
    (
        'FUT PERF PROG',
        [
            {'TAG': 'MD', 'LOWER': 'will'},
            non_verbs,
            {'TAG': {'IN': ['VB']}, 'LEMMA': 'have'},
            non_verbs,
            {'TAG': 'VBN', 'LEMMA': 'be'},
            {'TAG': 'VBG'},
        ]
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
            {'TAG': 'VBD', 'LEMMA': 'use'},
            {'LOWER': 'to'},
            {'TAG': 'VB'},
        ],
    ),

    # -- emphatics --
    (
        'PRES do-support',
        [
            {'TEXT': {'IN': ['does', 'do']}},
            {'TAG': 'VB', 'LEMMA': {'NOT_IN': ['evil']}},
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
        '?VB',
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
        'IMPV do not',
        [
            {'TEXT': 'Do'},
            {'LEMMA': 'not'},
            non_verbs,
            {'TAG': 'VB'}
        ],
    ),
    (
        'IMPV do not',
        [
            {'TEXT': 'do', 'IS_SENT_START': True},
            {'LEMMA': 'not'},
            non_verbs,
            {'TAG': 'VB'}
        ],
    ),
    (
        'IMPV do not',
        [
            {'IS_PUNCT': True},
            {'TEXT': 'do'},
            {'LEMMA': 'not'},
            non_verbs,
            {'TAG': 'VB'},
        ],        
    ),
    (
        'IMPV do not',
        [
            {'TAG': {'IN': ['CC']}, 'IS_TITLE': True},
            advbs,
            {'TEXT': 'do'},
            {'LEMMA': 'not'},
            non_verbs,
            {'TAG': 'VB'}
        ],
    ),
    (
        'IMPV do not',
        [
            {'IS_PUNCT': True},
            {'TAG': {'IN': ['CC']}},
            advbs,
            {'TEXT': 'do'},
            {'LEMMA': 'not'},
            non_verbs,
            {'TAG': 'VB'}
        ],
    ),

    # -- modals --
     (
        'MOD lest',
        [
            {'lower': 'lest'},
            non_verbs,
            {'TAG': {'IN': ['VBP', 'VB']}},
        ],
    ),
    (   
        'MOD need not',
        [   
            {'TEXT': 'need'},
            {'TEXT': 'not'},
            {'TAG': 'VB'},
        ],  
    ),
    (
        'MOD that be',
        [
            {'TEXT': 'that'},
            non_verbs,
            {'TEXT': 'be'},
        ]
    ), 
    (
        'MOD is to',
        [
            {'TEXT': {'IN': ['is', 'are']}},
            {'LEMMA': 'not', 'OP': '?'},
            {'TEXT': 'to'},
            {'TAG': 'VB'},
        ],
    ),
    (
        'MOD was to',
        [
            {'TEXT': {'IN': ['was', 'were']}},
            {'TEXT': 'to'},
            {'TAG': 'VB'},
        ],
    ),
]

# add a series of modal forms
modal_words = ['let', 'may', 'must', 'might']
quest_modals = ['shall', 'should', 'would', 'can', 'could']
modal_words = modal_words + quest_modals

for verb in modal_words:
    rules.extend([
        (
            f'MOD {verb}',
             [
                {'TAG': {'IN':['VB', 'MD']}, 'LOWER': verb},
                non_verbs,
                {'TAG': 'VB'},
            ]
        ),
        (
            f'?MOD {verb}',
             [
                {'TAG': {'IN':['VB', 'MD']}, 'LOWER': verb},
                non_verbs,
                {'TAG': 'VBP'},
            ]
        ),
    ])

for verb in quest_modals:
    rules.append(
         (
            f'MOD quest {verb}',
             [
                {'TAG': {'IN':['VB', 'MD']}, 'LOWER': verb, 'IS_TITLE': True},
                non_verbs,
                {'TAG': {'IN': ['VBP', 'VB']}, 'LOWER': {'NOT_IN': ['turn']}},
            ]
        )
)
