import matplotlib.pyplot as plt
import seaborn as sns

def plot_bar_1D(data_1D, ax=None, title='', xlabel='', ylabel='', **plot_kwargs):
    """Make barplot from 1D data to show simple counts."""
    if not ax: 
        fig, ax = plt.subplots()
    if not {'color', 'cmap', 'edgecolor'} & set(plot_kwargs):
        plot_kwargs.update({'color': 'lightgrey', 'edgecolor':'black'})
    data_1D.plot(kind='bar', ax=ax, **plot_kwargs)
    ax.set_title(title)
    ax.grid(axis='y')
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def heatmap(df, **kwargs):
    """Draw seaborne heatmap with custom settings"""
    default_kwargs = { 
        'center': 1.3,
        'cmap': sns.diverging_palette(220, 10, as_cmap=True),
        'square': True,
        'linewidth': 0.5,
        'annot': True,
    }   
    default_kwargs.update(kwargs)
    sns.heatmap(df, **default_kwargs)

