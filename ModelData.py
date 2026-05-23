from pathlib import Path
current_dir = Path(__file__).parent

import numpy as np
import pandas as pd
import json
from dataclasses import dataclass


def load_json_dict(filename):
    """Load a JSON file and return as dictionary."""
    with open(f'{current_dir}/{filename}.json', 'r') as json_file:
        return json.load(json_file)


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
    diet_import: float
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


def read_json_SpeciesGroup_list(filename):
    """Load JSON file and deserialize to list of SpeciesGroup objects."""
    def decoder(dct):
        return SpeciesGroup.from_dict(dct)

    with open(f'{current_dir}/{filename}.json', 'r') as f:
        return json.load(f, object_hook=decoder)


def get_seq2name(model_number, model_diet_datas):
    """Get mapping from group sequence number to group name."""
    model = model_diet_datas[str(model_number)]
    if not isinstance(model, dict):
        return None
    
    groups = model.get('group', {})
    if len(groups) == 0:
        return None
    
    seq2name = {int(g["group_seq"]): g["group_name"] for g in groups}
    import_seq = max(seq2name) + 1
    seq2name[import_seq] = 'diet_import'

    return seq2name


def get_DC(model_number, model_diet_datas):
    """Get diet composition (DC) and detritus fate matrices."""
    model = model_diet_datas[str(model_number)]
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

    DC = DC.sort_index().sort_index(axis=1)

    detritus_fate = pd.DataFrame.from_dict(detritus_fate_dict, orient='index').fillna(0)
    detritus_fate = detritus_fate.sort_index().sort_index(axis=1)

    return DC, detritus_fate


class ModelData:

    def __init__(self, model_number: int):
        def _groups_data(model_number):
            def to_df_row(c: SpeciesGroup):
                dct = c.__dict__.copy()
                dct.pop("taxons_included")
                return pd.DataFrame([dct])
            
            df_rows = []
            taxons = {}
            for i in range(len(species_groups)):
                if species_groups[i].model_number == model_number:
                    model_index = i
                    row = to_df_row(species_groups[i])
                    df_rows.append(row)
                    taxons[species_groups[i].group_seq] = species_groups[i].taxons_included

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
            df['egestion'] = df['q'] * df['gs']

            df['flow_to_det'] = df['egestion'] + df['M0']

            cols_to_return = ['group_name', 'trophic_info', 'tl', 'ge', 'ee', 'catch', 
                'biomass', 'pb', 'qb', 'p', 'q', 'predation', 'M0', 'gs', 'egestion', 'respiration', 'biomass_accum', 'emigration', 'immigration', 'net_migration',
                'flow_to_det', 'detritus_import'    
            ]

            groups_data = df.set_index('group_seq')[cols_to_return]
            groups_data = groups_data.sort_index(ascending=False)
            diet_import = df.set_index('group_seq')['diet_import'].sort_index(ascending=False)

            return groups_data, taxons, model_index, diet_import

        self.groups_data: pd.DataFrame
        self.groups_taxons: dict[int, list[dict]]  # group_seq: taxons_list
        
        self.groups_data, self.groups_taxons, model_index, diet_import = _groups_data(model_number)
        model_species_group_example = species_groups[model_index]

        self.model_number: int = model_number
        self.model_name: str = model_species_group_example.model_name
        self.model_country: str = model_species_group_example.model_country
        self.model_year: str = model_species_group_example.model_year
        self.lme: int = model_species_group_example.lme
        
        self.seq2name = get_seq2name(model_number, model_diet_datas)
        self.name2seq = {v: k for k, v in self.seq2name.items()}
        import_seq = max(self.seq2name)

        self.groups_data.loc[import_seq] = np.nan
        self.groups_data.loc[import_seq, 'group_name'] = 'diet_import'
        self.groups_data.loc[import_seq, 'trophic_info'] = 'Import'
        self.groups_data.loc[import_seq, 'tl'] = 1.0
        self.groups_data.loc[import_seq, 'respiration'] = 0.0
        diet_import.loc[import_seq] = 0

        # DC and Detritus fate:
        DC, det_fate = get_DC(model_number, model_diet_datas)
        DC.loc[import_seq] = 0  # add import data to DC
        DC[import_seq] = diet_import  # add import data to DC
        self.DC = DC.sort_index(ascending=False).sort_index(axis=1, ascending=False)
        self.det_fate = det_fate.sort_index(ascending=False).sort_index(axis=1, ascending=False)


# Load Data:
model_metadatas = load_json_dict('real_models/model_metadatas')
model_datas = load_json_dict('real_models/model_datas')
model_diet_datas = load_json_dict('real_models/model_diet_datas')
species_groups = read_json_SpeciesGroup_list("real_models/SpeciesGroups")


def get_model_metadata(model_id):
    return model_metadatas[str(model_id)]


def get_model_data(model_id):
    return model_datas[str(model_id)]


def get_model_diet_data(model_id):
    return model_diet_datas[str(model_id)]
