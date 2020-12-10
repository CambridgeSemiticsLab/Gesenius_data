import sys
from pathlib import Path
repo = Path.
tools = Path.home().joinpath()
sys.path.append('')
from clause_relas import in_dep_calc as clause_relator

# load BHSA features with genre module
locations = [
    '~/github/etcbc/bhsa/tf/c', 
    '~/github/etcbc/genre_synvar/tf/c',
    '~/github/etcbc/valence/tf/c'
]
TF = Fabric(locations)
extra_features = '''
domain txt ps gn 
nu genre sense
mother sp
'''
features = tf_tools.standard_features + extra_features
api = TF.load(features)
bhsa = use('bhsa', api=api)
F, E, T, L, Fs, = bhsa.api.F, bhsa.api.E, bhsa.api.T, bhsa.api.L, bhsa.api.Fs
