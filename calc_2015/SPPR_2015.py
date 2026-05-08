from pathlib import Path
current_dir = Path(__file__).parent

import numpy as np
import pandas as pd
import json
from dataclasses import dataclass

def load_json_dict(filename):
    # Open a file in read mode ('r')
    with open(f'{current_dir}/{filename}.json', 'r') as json_file:
        # Use json.load() to read the data and convert it back to a dictionary
        return json.load(json_file)
    
def get_model_metadata(model_id):
    return model_metadatas[str(model_id)]

def get_model_data(model_id):
    return model_datas[str(model_id)]

def get_model_diet_data(model_id):
    return model_diet_datas[str(model_id)]

@dataclass
class SpeciesGroup:
    group_name: str
    group_seq: int
    model_number: int
    model_name: str
    model_country: str
    model_year: str
    lme: int
    tl: float
    biomass: float
    pb: float
    qb: float
    ee: float
    M0b: float
    flow_to_det: float
    gross_efficiency: float
    prop_unassimilated_food: float
    respiration: float
    biomass_accum: float
    biomass_accum_rate: float
    immigration: float
    emigration: float
    export: float
    trophic_info: str
    detritus_import: float
    taxons_included: list[dict] 

    @classmethod
    def from_dict(cls, data: dict):
        """Factory method to create an instance from a dictionary."""
        return cls(**data)

    def __str__(self):
        return f"[{self.model_number}] '{self.model_name}' {self.model_country} ({self.model_year}) -> {self.group_seq}: {self.group_name} | tl: {self.tl}, taxons included: {[self.taxons_included[i]['taxon_name'] + ' (' + self.taxons_included[i]['AphiaID'] + ')' for i in range(len(self.taxons_included))]}"
    
    def to_df_row(self):
        dct = self.__dict__.copy()
        dct.pop("taxons_included")
        return pd.DataFrame(dct)

def decoder(dct):
    return SpeciesGroup.from_dict(dct)

def read_json_SpeciesGroup_list(filename):
    # Open the JSON file and load the data
    with open(f'{current_dir}/{filename}.json', 'r') as f:
        # Use json.load and specify the object_hook
        return json.load(f, object_hook=decoder)

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

def get_seq2name(model_number):
    model = get_model_diet_data(model_number)
    if not isinstance(model, dict):
        return None
    
    groups = model.get('group', {})
    if len(groups) == 0:
        return None
    
    return {int(g["group_seq"]): g["group_name"] for g in groups}

def get_DC(model_number):
    model = get_model_diet_data(model_number)
    if not isinstance(model, dict):
        return None
    
    groups = model.get('group', {})
    if len(groups) == 0:
        return None
        
    DC_dict = {}
    detritus_fate_dict = {}

    for g in groups:
        diet_descr = g.get('diet_descr', {})
        diet_descr = diet_descr if diet_descr else {}
        diet = diet_descr.get('diet', None)
        if not diet:
            DC_dict[int(g['group_seq'])] = {int(g['group_seq']): None}
            detritus_fate_dict[int(g['group_seq'])] = {int(g['group_seq']): None}
        else:
            diet = diet if isinstance(diet, list) else [diet]
            DC_dict[int(g['group_seq'])] = {int(d['prey_seq']): float(d['proportion']) for d in diet}
            detritus_fate_dict[int(g['group_seq'])] = {int(d['prey_seq']): float(d['detritus_fate']) for d in diet}
    
    # add missing columns:
    DC = pd.DataFrame.from_dict(DC_dict, orient='index').fillna(0)
    for g in groups:
        if int(g['group_seq']) not in DC.columns:
            DC[int(g['group_seq'])] = 0

    DC = DC.sort_index(ascending=False).sort_index(axis=1)

    detritus_fate = pd.DataFrame.from_dict(detritus_fate_dict, orient='index').fillna(0)
    detritus_fate = detritus_fate.sort_index(ascending=False).sort_index(axis=1)

    return DC, detritus_fate


def SPPR_2015(model_number):

    # general data:
    seq2name = get_seq2name(model_number)
    groups_data = get_model_groups_data(model_number).fillna(0)

    # DC and detritus_fate matrices:
    DC, det_fate = get_DC(model_number)

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


# Load Data:
model_metadatas = load_json_dict('model_metadatas')
model_datas = load_json_dict('model_datas')
model_diet_datas = load_json_dict('model_diet_datas')
species_groups = read_json_SpeciesGroup_list("SpeciesGroups")