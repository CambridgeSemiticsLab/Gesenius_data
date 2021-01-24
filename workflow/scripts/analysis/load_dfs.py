from pathlib import Path
import pandas as pd

class DfLoader:

    """Loads tables of data with custom selections.

    Args:   
        csv_dir: directory containing csv files
        csvs: str with kleene star or list of files to load
    """

    def __init__(self, csv_dir, csvs='*.csv'):

        csv_dir = Path(csv_dir)

        load_order = {
            'bhsa': 1, 'bhsa_clrela': 2,
            'eng': 3, 'eng_text': 4, 'lxx': 5,
        }
        sort_key = lambda f: load_order.get(f.stem, max(load_order.values())+1)
        # get iterable of table files
        if type(csvs) == str: 
            csvs = sorted(csv_dir.glob(csvs), key=sort_key)
        
        # populate list with pandas dataframes
        tables = [] 
        for file in csvs:
            tables.append(pd.read_csv(file, index_col='bhsa_node'))
            
        # concatenate into single df
        self.df = pd.concat(tables, 1)

        # cache dataframes that are already called
        self.cache = {}

    # Wrapper function for retrieving a df or its cache.
    # the ability to cache allows functions to depend on 
    # other dataframes in any order without wasting memory;
    # this wrapper function first attempts to retrieve the 
    # dataframe from the self.cache dict using the function's name;
    # it will run the function once if it is not found and store
    # the resulting dataframe in self.cache
    def get_df(get_func):

        def get_cache(self):

            # check if DF in cache
            name = get_func.__name__
            if name in self.cache:
                return self.cache[name]

            # build DF for the first time
            else:
                df = get_func(self)
                self.cache[name] = df
                return df

        return get_cache

    @get_df
    def eng_agree(self):
        """Get DF where english translations agree."""
        df = self.df
        df = df[df.eng_agree == 1]
        return df

    @get_df
    def eng_simp_agree(self):
        """Get DF where english translations agree."""
        df = self.df
        df = df[df.eng_simp_agree == 1]
        return df

    @get_df
    def esv(self):
        df = self.df
        df = df[df.esv_TAM.str.match('.*', na=False)]
        return df

    @get_df
    def niv(self):
        df = self.df
        df = df[df.niv_TAM.str.match('.*', na=False)]
        return df

    @get_df
    def eng_both(self):
        df = self.df
        df = df[
            (df.esv_TAM.str.match('.*', na=False))
            & (df.niv_TAM.str.match('.*', na=False))
        ]
        return df

    @get_df
    def eng_disagree(self):
        df = self.eng_both()
        df = df[df.eng_agree == 0]
        return df

    @get_df
    def eng_simp_disagree(self):
        df = self.eng_both()
        df = df[df.eng_simp_agree == 0]
        return df

    @get_df
    def df_safe_qtl(self):
        """Gives a dataframe with safe data for qtl.
    
        ! DEPRECATED !
        This is only kept for qtl data, which uses
        a different set of configurations.
        """
                 
        # select with the 'safe' column
        dfs = self.df[self.df.safe]

        # exclude uses of "did not" for now since these 
        # have unique semantic considerations
        # remove cases of 'did not' for now since these are semantically ambiguous
        dfs = dfs[
            (~dfs.esv_TAMspan.str.match('.*did not.*', na=False))
            & (~dfs.niv_TAMspan.str.match('.*did not.*', na=False))
        ].copy()

        return dfs
