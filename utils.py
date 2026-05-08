import numpy as np
import sympy as sm
import pandas as pd

# Import ModelData components
from ModelData import (
    ModelData, 
    SpeciesGroup, 
    species_groups, 
    model_diet_datas,
    get_seq2name as _get_seq2name,
    get_DC as _get_DC
)


# SpeciesGroup, ModelData and related functions are imported from ModelData.py


def mat_from_np(x, resolution=None):
    """Convert numpy array to sympy Matrix with optional rational resolution."""
    x = np.array(x)

    if resolution is not None:
        if x.ndim == 1:
            rationals_list = [sm.Rational(str(element)).limit_denominator(resolution) for element in x]
        elif x.ndim == 2:
            rationals_list = [[sm.Rational(str(element)).limit_denominator(resolution) for element in row] for row in x]
        else:
            raise Exception(f'array has too many dimensions: {x.ndim}')
    else:
        if x.ndim == 1:
            rationals_list = [sm.Rational(str(element)) for element in x]
        elif x.ndim == 2:
            rationals_list = [[sm.Rational(str(element)) for element in row] for row in x]
        else:
            raise Exception(f'array has too many dimensions: {x.ndim}')
    return sm.Matrix(rationals_list)


def move_scattered_identity(df):
    """
    Identifies rows and columns that form an identity matrix,
    even if they are scattered, and moves them to the top-left.
    """
    identity_pairs = [] # List of (row_label, col_label)
    
    # Iterate through rows to find unit vectors
    for row_label, row in df.iterrows():
        # Check if row has exactly one '1' and the rest are '0'
        if np.sum(row == 1) == 1 and np.sum(row == 0) == (len(row) - 1):
            col_idx = np.where(row == 1)[0][0]
            col_label = df.columns[col_idx]
            identity_pairs.append((row_label, col_label))

    if not identity_pairs:
        return df, df.index, df.columns

    # Separate the pairs into row and column lists
    id_rows = [p[0] for p in identity_pairs]
    id_cols = [p[1] for p in identity_pairs]

    # Create the new order (Identity components first)
    new_row_order = id_rows + [r for r in df.index if r not in id_rows]
    new_col_order = id_cols + [c for c in df.columns if c not in id_cols]

    # Apply to both DataFrames
    return df.reindex(index=new_row_order, columns=new_col_order), new_row_order, new_col_order

def remove_cycles(Z_input, new=False):
    if new:
        return remove_cycles_new(Z_input)
    else:
        return remove_cycles_old(Z_input)

def remove_cycles_new(Z_input):
    # Work on a copy to avoid mutating the original matrix
    Z = Z_input.copy()
    
    # 1. Identify input type and convert to numpy for the logic
    is_pandas = isinstance(Z_input, pd.DataFrame)
    is_sm_matrix = isinstance(Z, sm.Matrix)
    if is_pandas:
        Z = Z_input.to_numpy(dtype=float, copy=True)
    elif is_sm_matrix:
        Z = sm.matrix2numpy(Z)
    else:
        Z = np.array(Z_input, dtype=float, copy=True)

    from remove_cycles_fix import remove_cycles as rc
    Z, _ = rc(Z)

    # 3. Re-wrap in DataFrame if necessary
    if is_pandas:
        return pd.DataFrame(Z, index=Z_input.index, columns=Z_input.columns)
    if is_sm_matrix:
        return mat_from_np(Z)
    
    return Z

def remove_cycles_old(Z_input):
    """
    Removes cycles from a flow matrix Z using the Ulanowicz method.
    Z[i, j] represents flow from prey i to consumer j.
    """
    # Work on a copy to avoid mutating the original matrix
    Z = Z_input.copy()
    
    # 1. Identify input type and convert to numpy for the logic
    is_pandas = isinstance(Z_input, pd.DataFrame)
    is_sm_matrix = isinstance(Z, sm.Matrix)
    if is_pandas:
        Z = Z_input.to_numpy(dtype=float, copy=True)
    elif is_sm_matrix:
        Z = sm.matrix2numpy(Z)
    else:
        Z = np.array(Z_input, dtype=float, copy=True)
    
    def find_cycle(matrix):
        """Finds a single cycle using DFS."""
        num_nodes = matrix.shape[0]
        visited = [False] * num_nodes
        stack = []
        parent = [-1] * num_nodes

        def dfs(u):
            visited[u] = True
            stack.append(u)
            for v in range(num_nodes):
                if matrix[u, v] > 0:
                    if v in stack: # Cycle detected
                        cycle_path = stack[stack.index(v):]
                        return cycle_path
                    if not visited[v]:
                        res = dfs(v)
                        if res: return res
            stack.pop()
            return None

        for i in range(num_nodes):
            if not visited[i]:
                cycle = dfs(i)
                if cycle: return cycle
        return None

    while True:
        cycle = find_cycle(Z)
        if not cycle:
            break  # No more cycles found
        
        # 1. Find the "Weakest Link" (Minimum flow in the cycle)
        # Create pairs of (from_node, to_node) for the cycle
        links = []
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            links.append((u, v))
            
        min_flow = min(Z[u, v] for u, v in links)
        
        # 2. Subtract the bottleneck flow from all links in the cycle
        for u, v in links:
            Z[u, v] -= min_flow
            # Clean up tiny floating point errors
            if Z[u, v] < 1e-20:
                Z[u, v] = 0

    # 3. Re-wrap in DataFrame if necessary
    if is_pandas:
        return pd.DataFrame(Z, index=Z_input.index, columns=Z_input.columns)
    if is_sm_matrix:
        return mat_from_np(Z)

    return Z


def get_model_groups_data(model_number):
    def to_df_row(c: SpeciesGroup):
        dct = c.__dict__.copy()
        dct.pop("taxons_included")
        return pd.DataFrame([dct])

    df_rows = [to_df_row(species_groups[i]) for i in range(len(species_groups)) if species_groups[i].model_number == model_number]
    df = pd.concat(df_rows, ignore_index=True)
    df = df.replace("-9999", np.nan).replace(-9999, np.nan)

    df = df.rename(columns={  # change names
        'export': 'catch',
        'prop_unassimilated_food': 'gs',
        'gross_efficiency': 'ge'
        })
    
    df['p'] = df['pb'] * df['biomass']  # production
    df['q'] = df['qb'] * df['biomass']  # consumption
    df['M0'] = df['p'] * (1-df['ee'])  # other mortality
    df['net_migration'] = df['emigration'] - df['immigration']  # net migration
    
      # production*EE = catch + predation + biomass_accum + net_migration:
    df['predation'] = df['p'] * df['ee'] - (df['catch'] + df['biomass_accum'] + df['net_migration'])
    df['egestion'] = df['q'] * df['gs']  # gs

    df['flow_to_det'] = df['egestion'] + df['M0']

    cols_to_return = ['group_name', 'trophic_info', 'tl', 'ge', 'ee', 'catch', 
        'biomass', 'pb', 'qb', 'p', 'q', 'predation', 'M0', 'gs', 'egestion', 'respiration', 'biomass_accum', 'emigration', 'immigration', 'net_migration',
         'flow_to_det', 'detritus_import',        
        ]

    return df.set_index('group_seq')[cols_to_return]


def SPPR_2015(model_number):

    # general data:
    seq2name = _get_seq2name(model_number, model_diet_datas)
    groups_data = get_model_groups_data(model_number).fillna(0)

    # DC and detritus_fate matrices:
    DC, det_fate = _get_DC(model_number, model_diet_datas)

    Z = DC.mul(groups_data['q'], axis='index')
    DET_seq = list(groups_data[groups_data['trophic_info'] == 'DET'].index.values)
    if len(DET_seq) == 1:
        Z.loc[DET_seq[0], :] = groups_data['flow_to_det']
    elif len(DET_seq) == 2:  # det_date acts as a switch between DET groups
        Z.loc[DET_seq[0], :] = groups_data['flow_to_det'] * (DC * (1-det_fate)).sum(axis=1)
        Z.loc[DET_seq[1], :] = groups_data['flow_to_det'] * (DC * (det_fate)).sum(axis=1)
    else:
        raise Exception('too many detritus groups(?)')

    production = groups_data['p'].copy()

    # combine PP to single row:
    PP_seq_list = sorted(groups_data.index[groups_data['trophic_info'] == 'PP'].values)
    PP_seq = PP_seq_list[-1]
    seq_to_drop = PP_seq_list[:-1]
    if len(PP_seq_list) > 1:

        production.loc[PP_seq] += production.loc[seq_to_drop].sum()
        production = production.drop(index=seq_to_drop)

        Z.loc[:, PP_seq] += Z.loc[:, seq_to_drop].sum(axis=1)
        Z = Z.drop(index=seq_to_drop, columns=seq_to_drop)

    # combine part of DET that is PP into PP row:
    DET_seq = groups_data[groups_data['trophic_info'] == 'DET'].index.values[0]  # assuming there is only one DET
    Z_without_DET = Z.copy()
    percent_of_det_that_is_PP = Z_without_DET.loc[DET_seq, PP_seq] / Z_without_DET.loc[DET_seq, :].sum()  # 98%
    Z_without_DET.loc[:, PP_seq] += percent_of_det_that_is_PP * Z_without_DET.loc[:, DET_seq]
    Z_without_DET = Z_without_DET.drop(index=DET_seq, columns=DET_seq)

    # production of living compartments:
    new_index = Z_without_DET.index
    P = production[new_index].copy()
    ee = groups_data.loc[new_index, "ee"]
    non_PP = groups_data['trophic_info'] == 'Regular'
    P[non_PP] = P[non_PP].mul(ee[non_PP])   # P*EE = (export + predation + growth + net_migration), without M0
    # EE = (export + predation + growth + net_migration) / (export + predation + growth + net_migration + M0)
    P = P.sort_index(ascending=False)

    # production-normalized transaction matrix:
    A = Z_without_DET.T.sort_index(ascending=False).sort_index(axis=1, ascending=False) / P

    # production requirement matrix:
    seq_to_drop = P.index[P == 0]
    A = A.drop(columns=seq_to_drop, index=seq_to_drop)
    L = pd.DataFrame(np.linalg.inv((np.identity(A.shape[0]) - A)), index=A.index, columns=A.columns)

    SPPR = L.loc[PP_seq, :]

    return SPPR, seq2name, DC, Z, P, A, L, groups_data, PP_seq, DET_seq, new_index


# Data loading is now handled in ModelData.py