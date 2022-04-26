import pandas as pd
from itertools import product
from typing import NamedTuple

url = 'https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html'


def make_dict(df):
    ''' Make scoring dictionary from given dataframe values '''

    perms = [''.join(p) for p in product(list(df.columns), repeat=2)]
    dct = {a_combo: 0 for a_combo in perms}
    for key in dct.keys():
        a1, a2 = key[0], key[1]
        dct[key] = round(abs(float(df[a1]) - float(df[a2])) * -1, 2)
    return dct


def clean_df(df):
    '''Cleans df and adds property name as a column value'''

    df.columns = df.iloc[0]
    df.drop([0, 2], inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    df['*'] = [0]
    return df


# Make clean dfs for different properties
df_hydropathy = clean_df(pd.read_html(url)[0])
df_hydropathy.rename(columns={'W (1)': 'W'}, inplace=True)

df_volume = clean_df(pd.read_html(url)[1])


class AminoAcid(NamedTuple):
    hydropthy_dict: dict = make_dict(df_hydropathy)
    volume_dict: dict = make_dict(df_volume)
