import csv
import json
import collections
from pathlib import Path
from tf.fabric import Fabric
from book_formats import get_book_maps
from verb_form import get_verbform, get_cl_verbform
from modify_domain import permissive_q
from synvar_carc import in_dep_calc as clause_relator
from modify_cltype import simplify_cl_type

# fire up Text-Fabric with BHSA data
TF = Fabric(snakemake.input['tf_mods'], silent='deep')
features = """
sp pdp vs vt ps gn nu
lex language gloss voc_lex voc_lex_utf8
function number label 
typ code rela mother domain txt 
genre
sense
"""
bhsa = TF.load(features, silent='deep')
F, E, T, L, Fs, = bhsa.F, bhsa.E, bhsa.T, bhsa.L, bhsa.Fs

bookmap = get_book_maps(bhsa)

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
    ref_string = f'{book} {chapter}:{verse}'
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

    # collect preceding particles only
    particle_types = {'advb', 'prep', 'conj', 'prde', 'prin', 'inj', 'inrg'}
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
            'verb_form': get_verbform(node, bhsa),
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
            'prec_lexes': prec_lexes,
            'prec_pos': prec_pos,
            'prec_part': prec_particles,
    }
    return row_data 

def nearby_clatom_data(clatom_lookup):
    """Retrieve data on a nearby clause_atom, if it exists
    
    Args: 
        clatom_lookup: iterable of clause_atom nodes or empty
    Returns:
        dict of data on the first clause_atom in the lookup, if
        one was found, else an empty dict 
    """
    rel_dat = {}
    # retrive data on first clause in the lookup; if there is one
    if clatom_lookup:
        cl_atom = rel_dat['cl_atom'] = clatom_lookup[0]
        cl = rel_dat['cl'] = L.u(cl_atom, 'clause')[0]
        rel_dat['typ'] = F.typ.v(cl_atom)
        rel_dat['verb_type'] = get_cl_verbform(cl_atom, bhsa)
        rel_dat['domain2'] = permissive_q(cl, bhsa) # domain with permissive Q
        rel_dat['rela'] = clause_relator(cl, bhsa)
        rel_dat['txt'] = T.text(cl_atom)
    return rel_dat

def clrela_row(node):
    """Retrieve data on related clauses."""
     
    clause_atom = L.u(node, 'clause_atom')[0]

    # build data on the mother/daughter clause
    relas = {
        'mother': nearby_clatom_data(E.mother.f(clause_atom)),
        'daught': nearby_clatom_data(E.mother.t(clause_atom))
    }

    row_data = {'bhsa_node': node}
    for relcl, rcdata in relas.items():
        row_data.update({
            f'{relcl}_clause': rcdata.get('txt', ''),
            f'{relcl}_type': rcdata.get('typ', ''),
            f'{relcl}_verbtype': rcdata.get('verb_type', ''),
            f'{relcl}_rela': rcdata.get('rela', ''),
            f'{relcl}_domain2': rcdata.get('domain2', ''),
    })

    return row_data

# tables to make with their matched
# row-building functions
table_builds = {
    'bhsa': main_row,
    'bhsa_clrela': clrela_row,
}

# iterate through all samples and build up the data
for samp_set in snakemake.input['samples']:

    set_file = Path(samp_set)    
    verb_form = set_file.stem
    nodes = sorted(json.loads(set_file.read_text()))
    table_data = collections.defaultdict(list)

    for node in nodes:
        for table, row_getter in table_builds.items(): 
            table_data[table].append(row_getter(node)) 

    # make the exports
    outdir = Path(snakemake.params.outdir)
    for table, rows in table_data.items():
        file = outdir.joinpath(f'{verb_form}/{table}.csv')
        with open(file, 'w') as outfile:
            header = rows[0].keys()
            writer = csv.DictWriter(outfile, header)
            writer.writeheader()
            writer.writerows(rows)
