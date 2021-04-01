"""
Use Spacy Matcher object to match English tense constructions
"""

import json
import pickle
import collections
import textwrap
from pathlib import Path
import spacy
from spacy.tokens import Doc, Token, Span
from spacy.matcher import Matcher
from spacy.util import filter_spans  # filter overlaps; nice tip: https://stackoverflow.com/a/63303480/8351428
from spacy.gold import align # align different tokenizations: https://spacy.io/usage/linguistic-features#aligning-tokenization
from tf.core.timestamp import Timestamp

# custom modules
from tense_rules import rules as tense_rules
from imperatives import test_span_impv
from gbi_functions import id2ref

ts = Timestamp()

# load the parsing data
def load_pickle(path):
    with open(path, 'rb') as infile:
        return pickle.load(infile)
english_verbs = load_pickle(snakemake.input.eng_verbs)
parsed_verses = load_pickle(snakemake.input.parsedverses)
word_data = load_pickle(snakemake.input.word_data)
verse2words = load_pickle(snakemake.input.verse2words)

# A "span" is a Spacy object, which is here
# returned from the matcher object
def attach_span(spans):
    """Connect a token to its span explicitly."""
    for span in spans:
        for token in span:
            token._.tense_span = span

def bequeath_tense(span, affix='?'):
    """Pass on a tense category to conjunction children."""

    #if span._.tense_tag != tense:
    #    return None

    # prevent recursive string labels
    tense = span._.tense_tag
    if tense.startswith('c-'):
        return None

    for token in span:
        for child in token.conjuncts:
            if all([
                len(child._.tense_span or []) == 1,
                child.tag_ == 'VB',
                child.i > token.i, # prevent backwards inheritance, e.g. Ex 28:1 see parsed
            ]):
                child._.tense_span._.tense_tag = affix+tense

# NB: use this to apply more advanced rules to spans
def correct_span(span):
    """Apply corrective checks to a span."""
    #if span._.tense_tag == 'IMPV':
    #    if not test_span_impv(span):
    #        span._.tense_tag = 'VB not-impv' # remove the label
    pass

def match_spans(parsed_verses, matcher):
    """For every verse, apply custom matcher rules and
     isolate the set of relevant spans which match the tense 
    rules and map to verse reference tuples. The identified 
    spans can then later be matched with words from the verbs
    dictionaries.
    """
    verse2spans = collections.defaultdict(dict)
    for trans, ref_tuples in parsed_verses.items():
        for ref_tuple, spacy_doc in ref_tuples.items():
            matches = matcher(spacy_doc)
            
            # retrieve Spacy Span objects
            # and give them tense tags
            spans = []
            for m_id, start, end in matches:
                span = spacy_doc[start:end]
                span._.tense_tag = nlp.vocab.strings[m_id]
                correct_span(span)
                spans.append(span)
            
            filtered_spans = filter_spans(spans)  # filter out overlapping spans; keep longest
            attach_span(filtered_spans) # ensure tokens are mapped to their matched span
            for span in filtered_spans:
                bequeath_tense(span) 
            
            # save positive matches; unmatched verses will
            # be recognized later
            if filtered_spans:
                verse2spans[trans][ref_tuple] = filtered_spans
            else:
                continue

    return verse2spans

def trans_to_span(para_words, spans, 
                  verse_words, aligner):
    """Match given words with its tense span.
    
    A match is an overlap of known parallel words and 
    a span of matched words, based on overlapping GBI ids.
    Thus all Spacy token indicies are converted to GBI indicies
    and used to lookup the corresponding GBI ids for set comparison.
    
    Args:
        para_words: list of gbi word ids for a known parallel alignment 
        spans: list of Spacy Span objects from the Matcher, with tense_tag attributes
        verse_words: list of gbi ids within a verse; is indexed with remapped indices
            from the Span tokens, which have attributes `start` and `end` which
            correspond with their index in the Spacy doc (verse). Those indices are
            remapped to their GBI positions with the `aligner`.        
        aligner: remaps Spacy indicies to GBI indicies for tokens
    """
    for span in spans:
        start, end = aligner(span.start), aligner(span.end-1) 
        end += 1 # -1 above to avoid IndexError since end might be +1 longer than end for index slicing
        span_words = set(verse_words[start:end])
        if set(para_words) & span_words:
            return span

class InspectionDoc:

    def __init__(self):
        self.data = collections.defaultdict(lambda: collections.defaultdict(str))

    def make(self, esv_path, niv_path):
        # export inspection file
        for trans, path in [('esv', esv_path), ('niv', niv_path)]:
            doc = ''
            path = Path(path)
            for verse, message in self.data[trans].items():
                doc += '{} {}:{}'.format(*verse) + '\n'
                verse_text = str(parsed_verses[trans][verse])
                verse_text = '\n'.join(textwrap.wrap(verse_text, 80)) + '\n'
                doc += verse_text
                doc += message
                doc += '\n'
            path.write_text(doc)

def link_spans():
    """Match spans with verbs and assign data.

    The primary challenge is that there is no easy
    way to link the Spacy-parsed doc data with the
    word data we have for GBI. Remember that the full
    verse texts are compiled from the word lists and 
    then parsed as a full verse. The challenge is to
    use the indexing from before the word lists were
    compiled to match with the indexing of the Spacy
    doc object. That then needs to be cross-referenced
    with the Spacy Matcher object.
    """

    ts.indent(0, reset=True)
    ts.info('matching spans...')
            
    inspect = InspectionDoc()
    bhsa2eng = collections.defaultdict(dict)

    for trans, bhsa_nodes in english_verbs.items():
        
        for bhsa_node, para_words in bhsa_nodes.items():

            # get GBI-side data
            verse_ref = id2ref(para_words[0], 'translation')        
            para_text = ' '.join(word_data[trans][w]['text'] for w in para_words)
            verse_words = verse2words[trans][verse_ref]
            verse_tokens = [word_data[trans][w]['text'] for w in verse_words]
            verse_tokens = [t.replace(';','.') for t in verse_tokens]
            
            # get Spacy-side data
            verse_parsing = parsed_verses[trans][verse_ref]
            spacy_tokens = [str(t) for t in verse_parsing]
            
            # map Spacy tokens back to GBI tokens using indicies
            # Spacy tokenizes words with apostrophes differently (for e.g. `he'll` == `he` + `'ll`)
            # They can be re-aligned: https://spacy.io/usage/linguistic-features#aligning-tokenization
            cost, a2b, b2a, a2b_multi, b2a_multi = align(spacy_tokens, verse_tokens) # alignment of indicies here
            aligner = lambda i: a2b_multi.get(i, a2b[i]) # returns 1-to-1 or many-to-1 aligned index
            
            # try to retrieve span links with advanced tense tags
            verse_sents = list(verse_parsing.sents)
            spans = verse2spans[trans].get(verse_ref, [])
            span_match = trans_to_span(para_words, spans, verse_words, aligner) or '' # search for overlapping GBI id sets
            if span_match:
                tense_tag = span_match._.tense_tag
                sentence_i = verse_sents.index(span_match[-1].sent) 
            else:
                tense_tag = ''
                sentence_i = None
            
            # retrieve basic parsings
            raw_tokens = []
            for i, token in enumerate(verse_parsing):
                if verse_words[aligner(i)] in para_words:
                    raw_tokens.append(token)
                    
            vb_tokens = [t for t in raw_tokens if t.tag_.startswith('VB')]

           
                
            # save the data
            data = {
                'eng_ref': verse_ref,
                'words': para_text,
                'tags': '|'.join(t.tag_ for t in raw_tokens),
                'vb_tags': '|'.join(t.tag_ for t in vb_tokens),
                'tense': tense_tag,
                'tense_span': f'{span_match}',
                'sentence_i': sentence_i,
            }
            
            bhsa2eng[trans][bhsa_node] = data
                
            # add strings to inspection file
            if span_match and span_match._.tense_tag:
                inspect.data[trans][verse_ref] += f'\t\tMATCH: {bhsa_node}|{tense_tag}|{span_match}|{para_text}\n'
            else:
                inspect.data[trans][verse_ref] += f'\t\tMISS: {bhsa_node}|''|''|{para_text}\n'

                
    ts.info('done with matches')
    return (bhsa2eng, inspect)

# -- BUILD AND MATCH TENSE SPANS --

# load custom tense rules
nlp = spacy.load('en_core_web_sm')
Span.set_extension('tense_tag', default='', force=True)
Token.set_extension('tense_span', default=None, force=True)
matcher = Matcher(nlp.vocab)
for tag, ruleset in tense_rules:
    matcher.add(tag, None, ruleset)

# run the span matcher
verse2spans = match_spans(parsed_verses, matcher)

# match the spans with verbs
bhsa2eng, inspect_doc = link_spans()

# export
for trans, trans_data  in bhsa2eng.items():
    trans_file = snakemake.output[trans]
    with open(trans_file, 'w') as outfile:
        json.dump(trans_data, outfile, ensure_ascii=False, indent=2)

# export inspection document
inspect_doc.make(snakemake.output.inspect_esv, snakemake.output.inspect_niv)
