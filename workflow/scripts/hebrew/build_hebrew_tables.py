import re
import sys
import csv
import json
import pickle
import collections
from pathlib import Path
from tf.fabric import Fabric
from book_formats import get_book_maps, etcbc2sbl, etcbc2abbr
from verb_form import get_verbform, get_cl_verbform
from modify_domain import permissive_q
from synvar_carc import in_dep_calc as clause_relator
from modify_cltype import simplify_cl_type
from tag_args import clause_objects, get_loca_assocs, clause_locas, clause_time, clause_args

# NB that working directory when script is executed is 
# /workflow; because we have some utilities that we want
# to run from above directory, we need to append it to path
sys.path.append('scripts')
from build_tables import build_sample_tables

# fire up Text-Fabric with BHSA data
TF = Fabric(snakemake.input['tf_mods'], silent='deep')
features = """
sp pdp vs vt ps gn nu
lex language gloss voc_lex voc_lex_utf8
function number label 
typ code rela mother domain txt 
genre
sense
nhead
funct_assoc
"""
bhsa = TF.load(features, silent='deep')
F, E, T, L, Fs, = bhsa.F, bhsa.E, bhsa.T, bhsa.L, bhsa.Fs

# load GBI Hebrew data
with open(snakemake.input.bhsa2gbi, 'rb') as infile:
    bhsa2gbi = pickle.load(infile)

# preprocess data
bookmap = get_book_maps(bhsa)
loca_lexs = get_loca_assocs(bhsa)

def join_on(nodes, jchar='_', default=''):
    """Join words on a char and ensure they are pre/appended with that char.
    
    The pre/appending provides easy-to-match word boundaries.
    """
    joined_string = f'{jchar}'.join(nodes)
    if not joined_string:
        return default
    else:
        return f'{jchar}{joined_string}{jchar}'

def get_preceding_words(node, context='clause'):
    """Retrieves words from before a verb within a context"""
    context_node = L.u(node, context)[0]
    context_words = L.d(context_node, 'word')
    prec_words = context_words[:context_words.index(node)]
    return prec_words

def main_row(node):
    """Compile all relevant BHSA data for a given node."""

    # data on this clause itself
    book, chapter, verse = T.sectionFromNode(node)
    booksbl = etcbc2sbl[book]
    bookabbr = etcbc2abbr[book]
    ref_string = f'{book} {chapter}:{verse}'
    ref_sbl = f'{booksbl} {chapter}:{verse}'
    ref_abbr = f'{bookabbr} {chapter}:{verse}'
    verse_node = L.u(node, 'verse')[0]
    clause_atom = L.u(node, 'clause_atom')[0]
    clause = L.u(node, 'clause')[0]
    sent = L.u(node, 'sentence')[0]
    clause_type = F.typ.v(clause)
    preceding_words = get_preceding_words(node)
    prec_lexes = join_on((F.lex.v(w) for w in preceding_words), default='Ø') 
    prec_pos = join_on((F.pdp.v(w) for w in preceding_words), default='Ø')
    domain2 = permissive_q(clause, bhsa)
    cl_type_simp = simplify_cl_type(clause_atom, prec_lexes, bhsa)
    cl_args = clause_args(node, bhsa)
    has_q = ('Q' in cl_args) * 1 # look for question particles
#    cl_args = re.match('.*V', cl_args)[0] # NB ignore post-verbal arguments

    # collect preceding particles only
    particle_types = {'nega', 'advb', 'prep', 'conj', 'prde', 'prin', 'intj', 'inrg'}
    prec_particles = join_on(
        (F.lex.v(w) for w in preceding_words
            if F.pdp.v(w) in particle_types)
    , default='Ø') 
    
    null_string = ''
    
    row_data = {
            'bhsa_node': node,
            'ref': ref_string, 
            'book': book, 
            'book_super': bookmap['super'].get(book, book),
            'canon_part': bookmap['tripart'][book],
            'period': bookmap['period'].get(book, ''),
            'genre': F.genre.v(verse_node),
            'domain2': domain2,
            'text_full': F.g_word_utf8.v(node),
            'text_plain': F.g_cons_utf8.v(node),
            'lex': F.lex_utf8.v(node),
            'lex_etcbc': F.lex.v(node),
            'gloss': F.gloss.v(node),
            'verb_form': get_verbform(node, bhsa, bhsa2gbi),
            'stem': F.vs.v(node),
            'person': F.ps.v(node),
            'gender': F.gn.v(node),
            'number': F.nu.v(node),
            'valence': F.sense.v(node),
            'clause_atom': T.text(clause_atom),
            'clause': T.text(clause),
            'sentence': T.text(sent),
            'txt_type': F.txt.v(clause),
            'clause_type': clause_type,
            'cltype_simp': cl_type_simp,
            'clause_rela': clause_relator(clause, bhsa),
#            'cl_args': cl_args,
            'is_question': has_q,
            'prec_lexes': prec_lexes,
            'prec_pos': prec_pos,
            'prec_part': prec_particles,
            'ref_sbl': ref_sbl,
            'ref_abbr': ref_abbr,
    }

    # provide clause argument data
    # objects
    row_data.update(
        clause_objects(node, clause_atom, clause, bhsa)
    )
    # locatives
    row_data.update(
        clause_locas(node, loca_lexs, bhsa)
    )
    row_data.update(
        clause_time(node, bhsa)
    )
    # convert to boolean 0 or 1 to avoid indexing
    # pivot tables with booleans
    row_data['has_objc'] = 1 * row_data['has_objc']
    row_data['has_loca'] = 1 * row_data['has_loca']

    return row_data 

def nearby_clatom_data(clatom_lookup, starting_clatom):
    """Retrieve data on a nearby clause_atom, if it exists
    
    Args: 
        clatom_lookup: iterable of clause_atom nodes or empty
    Returns:
        dict of data on the first clause_atom in the lookup, if
        one was found, else an empty dict 
    """
    rel_dat = {
        'clause':'', 'cl_atom': '', 'clause_atom':'', 
        'rela': '', 'domain2': '', 'verbtype': '',
        'type': '', 'verb_ps': '', 'verb_lex': '', 
        'verbplain': '', 'intertext': ''
    }
    # retrive data on first clause in the lookup; if there is one
    if clatom_lookup:
        cl_atom = rel_dat['cl_atom'] = clatom_lookup[0]
        cl = L.u(cl_atom, 'clause')[0]
        verb = next((w for w in L.d(cl_atom, 'word') if F.pdp.v(w) == 'verb'), 0)
        rel_dat['verb_lex'] = F.lex.v(verb)
        rel_dat['verb_ps'] = F.ps.v(verb)
        rel_dat['type'] = F.typ.v(cl_atom)
        rel_dat['verbplain'] = F.g_cons_utf8.v(verb)
        rel_dat['verbtype'] = get_cl_verbform(cl_atom, bhsa, bhsa2gbi)
        rel_dat['domain2'] = permissive_q(cl, bhsa) # domain with permissive Q
        rel_dat['rela'] = clause_relator(cl, bhsa)
        rel_dat['clause_atom'] = T.text(cl_atom)
        rel_dat['clause'] = T.text(cl)

        # capture text in between starting node and this one
        if cl_atom - starting_clatom <= 3:
            if cl_atom < starting_clatom:
                interm_clatoms = list(range(cl_atom, starting_clatom))
            else:
                interm_clatoms = list(range(starting_clatom+1, cl_atom+1))
            for cl in interm_clatoms:
                rel_dat['intertext'] += T.text(cl)

    return rel_dat

def clrela_row(node):
    """Retrieve data on related clauses."""
     
    clause_atom = L.u(node, 'clause_atom')[0]

    # build data on the mother/daughter clause
    relas = {
        'mother': nearby_clatom_data(E.mother.f(clause_atom), clause_atom),
        'daught': nearby_clatom_data(E.mother.t(clause_atom), clause_atom)
    }

    row_data = {'bhsa_node': node}
    for relcl, rcdata in relas.items():
        row_data.update({
            f'{relcl}_{k}': rcdata[k] for k in rcdata            
        })

    return row_data

rowmakers = [main_row, clrela_row]

build_sample_tables(
    rowmakers,
    snakemake.input.sample,
    snakemake.output
)
