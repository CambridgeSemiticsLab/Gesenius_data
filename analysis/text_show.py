import numpy as np

def pretty_hebrew(val):
    """Render Hebrew in a dataframe."""
    return 'font-size:20px; font-family: Times New Roman; text-align: right; max-width: 500px'

def get_spread(array, n): 
    """Retrieve an even spread of indices an array/Series.
    
    https://stackoverflow.com/a/50685454/8351428
    
    Args:
        array: either numpy array or Pandas Series 
            (with multiple indices allowed)
        n: number of indices to return
    Returns:
        indexed array or series
    """
    end = len(array) - 1 
    spread = np.ceil(np.linspace(0, end, n)).astype(int)
    indices = np.unique(spread)
    try:
        return array[indices]
    except KeyError:
        return array.iloc[indices]

class TextShower:
    """Show stylized text examples from a Dataframe of data"""
    
    def __init__(self, default=['ref', 'sentence', 'text_full'], 
                 stylize=['sentence', 'text_full']):
        self.default = default
        self.stylize = stylize

    def show(self, df, extra=[], spread=0):
        """Display text from pandas dataframe in a readable way."""
        original_shape = df.shape
        df = df[self.default + extra]
        if spread > 0:
            spread_i = get_spread(df.index, spread)
            df = df.loc[spread_i]
        print(f'showing {df.shape[0]} of {original_shape[0]}')
        return df.style.applymap(pretty_hebrew, subset=self.stylize)
