"""
Use Spacy to translate English Bible translations which are aligned with
the Hebrew Bible data.
"""

import sys
import json
import collections
import re
import pandas as pd
from pathlib import Path

# Spacy tools
import spacy
from spacy.matcher import Matcher
from spacy.tokens import Doc, Token, Span
from spacy.util import filter_spans # filter overlaps; nice tip: https://stackoverflow.com/a/63303480/8351428
from spacy.gold import align # align different tokenizations: https://spacy.io/usage/linguistic-features#aligning-tokenization

# custom modules 
from gbi_functions import id2ref
