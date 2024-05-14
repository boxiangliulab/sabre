import pandas as pd
import sys


ID = sys.argv[1]

result_df = None

def pool(s):
    return list(s)[0]

def replace(s):
    result = []


for i in range(1, 23):
    df = pd.read_csv('./output/{}/singular_cells/singular_cell_linkage_chr{}.txt'.format(ID, i), sep='\t', dtype={'geno':'str'})
    if df.empty: continue
    df['chr'] = 'chr{}'.format(i)
    df = df.groupby(['chr', 'var', 'correct', 'global_oppo_support', 'barcode']).agg({'geno':pool, 'support':'sum', 'oppo_support':'sum'}).reset_index()
    df.columns = ['chr', 'var', 'correct', 'global_oppo_support', 'barcode', 'geno', 'support', 'oppo_support']
    agg_result = df.groupby(['chr','var','correct','global_oppo_support']).agg({'geno':pool,'support':'sum', 'barcode':'count', 'oppo_support':'sum'}).reset_index()
    agg_result.columns = ['chr','var', 'correct', 'global_oppo_support', 'geno', 'support', 'barcode','oppo_support']
    if result_df is None:
        result_df = agg_result
    else:
        result_df = pd.concat([result_df, agg_result], ignore_index=True)

result_df.to_csv('./output/{}/singular_cells/result_singular_cells.csv'.format(ID),index=False)
