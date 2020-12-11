"""
Apply some custom corrections to the qtl tables.
Many of these corrections are in place temporarily
to weed out tendentious taggings in the dataset.
"""

import pandas as pd
input = snakemake.input

print('applying qatal corrections...')

raw_tables = [
    pd.read_csv(input.bhsa, index_col='bhsa_node'),
    pd.read_csv(input.bhsa_clrela, index_col='bhsa_node'),
    pd.read_csv(input.eng, index_col='bhsa_node'),
    pd.read_csv(input.eng_text, index_col='bhsa_node'),
    pd.read_csv(input.lxx, index_col='bhsa_node'),
]

bhsa, bhsa_clrela, eng, eng_text, lxx = raw_tables

# Use Pandas operations to weed out unwanted samples
# the list of unwanted samples will then be used to filter
# out the samples from all of the tables

#  all imperative tags for qatal are fallacious
impv_cases = eng[
    (eng.niv_TAM == 'PRES..IMPV') 
    | (eng.esv_TAM == 'PRES..IMPV')
].index

# certain words in English have forms for which the parser cannot well-distinguish
# between the simple past and the present; we mark these cases manually using 
# the LXX as an extra check
bad_past_re = r'.*([Pp]ut|[Ss]et|[Cc]ut|[Ll]ay|[Cc]ast|[Ss]pread|[Ss]pit|[Rr]ead|rid)|darken'

bad_past_esv = bhsa[
    (eng_text.esv.str.match(bad_past_re, na=False))
    & (eng.esv_TAM == 'PAST..IND')
    & (lxx.lxx_tm.str.match('.*(present|future)'))
].index

bad_past_niv = bhsa[
    (eng_text.niv.str.match(bad_past_re, na=False))
    & (eng.niv_TAM == 'PAST..IND')
    & (lxx.lxx_tm.str.match('.*(present|future)'))
].index

# mark bad irrealis values
bad_irrealis = bhsa[
    (
        (~eng.niv_TAM.isin(['PRES..IMPV', 'FUT..IND', 'PRES..MOD', 'PRES..IND']) & (eng.esv_TAM.isin(['PRES..IMPV'])))
        | (eng.niv_TAM.isin(['PRES..IMPV']) & (~eng.esv_TAM.isin(['PRES..IMPV', 'FUT..IND', 'PRES..MOD', 'PRES..IND'])))
    )
    & (~lxx.lxx_tm.str.match('.*(future|impv)', na=False))
].index

not_safe = [bad_past_esv, bad_past_niv, bad_irrealis, impv_cases]
not_safe_nodes = set(n for iset in not_safe for n in iset)

print(len(not_safe_nodes), 'nodes listed as unsafe...')

# add new columns to BHSA
# note that we really only need the column 
# here, not on every table
bhsa['safe'] = ~bhsa.index.isin(not_safe_nodes)

# export the tables
for table, outpath in zip(raw_tables, snakemake.output):
    table.to_csv(outpath, index=True)
