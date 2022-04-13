import pandas as pd

url = 'https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html'


def clean_df(df, property_name):
    '''Cleans df and adds property name as a column value'''

    df.columns = df.iloc[0]
    df.drop([0, 2], inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    df['Property'] = property_name
    df.insert(0, 'Property', df.pop('Property'))

    return df

# Make clean dfs for different properties
df_hydropathy = clean_df(pd.read_html(url)[0], 'hydropathy') #  Kyte-Doolitle scale
df_hydropathy.rename(columns={'W (1)': 'W'}, inplace=True)

df_volume = clean_df(pd.read_html(url)[1], 'volume') # Angstrom cubed

# Concatenate the dfs to get one df
df_property = pd.concat([df_hydropathy, df_volume])
print(df_property)



class AminoAcid():

    def __init__(self):

        #  Kyte-Doolitle scale
        #self.hydropathy
        pass