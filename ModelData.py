from pathlib import Path
current_dir = Path(__file__).parent

import numpy as np
import pandas as pd
import json
import re
from dataclasses import dataclass


def load_json_dict(filename):
    """Load a JSON file and return as dictionary."""
    with open(f'{current_dir}/{filename}.json', 'r') as json_file:
        return json.load(json_file)


# Old SpeciesGroup for backward compatibility with legacy data
@dataclass
class SpeciesGroupLegacy:
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
        return pd.DataFrame([dct])


# New SpeciesGroup for JSON API
@dataclass
class SpeciesGroup:
    group_name: str
    group_seq: int
    biomass: float
    pb: float
    qb: float
    ee: float
    M0b: float
    tl: float
    flow_to_det: float
    ge: float
    gs: float
    respiration: float
    biomass_accum: float
    biomass_accum_rate: float
    immigration: float
    emigration: float
    export: float
    trophic_info: str
    detritus_import: float
    diet_import: float
    taxon_descr: str

    @classmethod
    def from_dict(cls, data: dict):
        """Factory method to create an instance from a dictionary."""
        return cls(**data)

    def to_df_row(self):
        dct = self.__dict__.copy()
        return pd.DataFrame([dct])


def read_json_SpeciesGroup_list(filename):
    """Load JSON file and deserialize to list of SpeciesGroupLegacy objects."""
    def decoder(dct):
        return SpeciesGroupLegacy.from_dict(dct)

    with open(f'{current_dir}/{filename}.json', 'r') as f:
        return json.load(f, object_hook=decoder)


class ModelData:
    """ModelData class supporting both old API (model_number) and new API (json_filepath)."""

    def __init__(self, model_input):
        """Initialize ModelData from either a model_number (int) or json_filepath (str).
        
        Args:
            model_input: Either an int (model_number) for legacy API or a str (json_filepath) for new API
        """
        if isinstance(model_input, int):
            # Old API: load from pre-existing JSON files by model number
            self._init_from_model_number(model_input)
        elif isinstance(model_input, str):
            # New API: load from a JSON filepath
            self._init_from_json_filepath(model_input)
        else:
            raise TypeError(f"model_input must be int or str, got {type(model_input)}")
    
    def _init_from_json_filepath(self, json_filepath):
        """Initialize ModelData from a JSON filepath.
        
        Filename format: {model_number}_{model_name}_{model_years}.json
        Example: 227_Iceland_(1950).json
        """
        # Parse filepath
        filepath = Path(json_filepath)
        filename_no_ext = filepath.stem  # e.g., "227_Iceland_(1950)"
        
        # Extract model_number, model_name, and model_year from filename
        self.model_number, self.model_name, self.model_year = self._parse_filename(filename_no_ext)
        
        # Country is not available from the JSON, set to NaN
        self.model_country = np.nan
        
        # Load JSON data
        self.groups_data: pd.DataFrame
        self.groups_taxons: dict[int, list[dict]]  # group_seq: taxons_list
        
        self.data_json = ModelData.load_json_dict_static(str(json_filepath))
        species_groups = ModelData.get_species_groups(self.data_json)
        self.groups_data, diet_import = ModelData.get_groups_df(species_groups)

        self.lme = 13  # Default, could be extracted from metadata if needed
        
        self.seq2name = ModelData.get_seq2name(self.data_json)
        self.name2seq = {v: k for k, v in self.seq2name.items()}
        import_seq = max(self.seq2name)

        self.groups_data.loc[import_seq] = np.nan
        self.groups_data.loc[import_seq, 'group_name'] = 'diet_import'
        self.groups_data.loc[import_seq, 'trophic_info'] = 'Import'
        self.groups_data.loc[import_seq, 'tl'] = 1.0
        self.groups_data.loc[import_seq, 'respiration'] = 0.0
        diet_import.loc[import_seq] = 0

        # DC and Detritus fate:
        DC, det_fate = ModelData.get_DC(self.data_json)
        DC.loc[import_seq] = 0  # add import data to DC
        DC[import_seq] = diet_import  # add import data to DC
        self.DC = DC.sort_index(ascending=False).sort_index(axis=1, ascending=False)
        
        det_fate.loc[import_seq] = 0  # add diet_import to det_fate
        det_fate[import_seq] = 0  # add diet_import to det_fate
        det_fate.loc[self.groups_data['trophic_info']=='Import', self.groups_data['trophic_info']=='DET'] = 1
        self.det_fate = det_fate.sort_index(ascending=False).sort_index(axis=1, ascending=False)
    
    def _init_from_model_number(self, model_number: int):
        """Initialize ModelData from a model number (legacy API)."""
        def _groups_data(model_number):
            def to_df_row(c):
                dct = c.__dict__.copy()
                if "taxons_included" in dct:
                    dct.pop("taxons_included")
                if "model_number" in dct:
                    dct.pop("model_number")
                if "model_name" in dct:
                    dct.pop("model_name")
                if "model_country" in dct:
                    dct.pop("model_country")
                if "model_year" in dct:
                    dct.pop("model_year")
                if "lme" in dct:
                    dct.pop("lme")
                return pd.DataFrame([dct])
            
            df_rows = []
            taxons = {}
            model_index = None
            for i in range(len(species_groups)):
                if hasattr(species_groups[i], 'model_number') and species_groups[i].model_number == model_number:
                    model_index = i
                    row = to_df_row(species_groups[i])
                    df_rows.append(row)
                    if hasattr(species_groups[i], 'taxons_included'):
                        taxons[species_groups[i].group_seq] = species_groups[i].taxons_included

            if not df_rows:
                raise ValueError(f"No species groups found for model {model_number}")
                
            df = pd.concat(df_rows, ignore_index=True)
            df = df.replace("-9999", np.nan).replace(-9999, np.nan)

            df = df.rename(columns={  # change names
                'export': 'catch',
                'prop_unassimilated_food': 'gs',
                'gross_efficiency': 'ge'
                }, errors='ignore')
            
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
        
        det_fate.loc[import_seq] = 0  # add diet_import to det_fate
        det_fate[import_seq] = 0  # add diet_import to det_fate
        det_fate.loc[self.groups_data['trophic_info']=='Import', self.groups_data['trophic_info']=='DET'] = 1
        self.det_fate = det_fate.sort_index(ascending=False).sort_index(axis=1, ascending=False)

    @staticmethod
    def _parse_filename(filename_no_ext):
        """Parse filename to extract model_number, model_name, and model_year.
        
        Example input: "227_Iceland_(1950)"
        Returns: (227, "Iceland", "1950")
        """
        # Extract year from parentheses at the end
        year_match = re.search(r'\(([^)]+)\)$', filename_no_ext)
        if not year_match:
            raise ValueError(f"Could not extract year from filename: {filename_no_ext}")
        
        model_year = year_match.group(1)
        
        # Remove the year part to get the rest
        rest = filename_no_ext[:year_match.start()].rstrip('_')
        
        # Split by underscore: first is model_number, rest is model_name
        parts = rest.split('_', 1)  # Split on first underscore only
        if len(parts) != 2:
            raise ValueError(f"Could not parse model_number and model_name from: {rest}")
        
        model_number = int(parts[0])
        model_name = parts[1]
        
        return model_number, model_name, model_year
    
    @classmethod
    def load_json_dict_static(cls, filepath):
        """Load a JSON file and return as dictionary."""
        with open(filepath, 'r') as json_file:
            return json.load(json_file)
    
    @classmethod
    def get_seq2name(cls, model):
        """Get mapping from group sequence number to group name."""
        if not isinstance(model, dict):
            return None
        
        groups = model.get('group', {})
        if len(groups) == 0:
            return None
        
        seq2name = {int(g["group_seq"]): g["group_name"] for g in groups}
        import_seq = max(seq2name) + 1
        seq2name[import_seq] = 'diet_import'

        return seq2name

    @classmethod
    def get_DC(cls, model_json_dict):
        """Get diet composition (DC) and detritus fate matrices."""
        if not isinstance(model_json_dict, dict):
            return None
        else:
            model = model_json_dict
        
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
                DC_dict[int(g['group_seq'])] = {int(float(d['prey_seq'])): float(d['proportion']) for d in diet}
                detritus_fate_dict[int(g['group_seq'])] = {int(float(d['prey_seq'])): float(d['detritus_fate']) for d in diet}
        
        # add missing columns:
        DC = pd.DataFrame.from_dict(DC_dict, orient='index').fillna(0)
        for g in groups:
            if int(g['group_seq']) not in DC.columns:
                DC[int(g['group_seq'])] = 0

        DC = DC.sort_index().sort_index(axis=1)

        detritus_fate = pd.DataFrame.from_dict(detritus_fate_dict, orient='index').fillna(0)
        detritus_fate = detritus_fate.sort_index().sort_index(axis=1)

        return DC, detritus_fate
    
    @classmethod
    def get_species_groups(cls, model_json):
        """Extract species groups from JSON model data."""
        groups = model_json.get('group', {})
        if len(groups) == 0:
            return []
        species_groups = []
        for j in range(len(groups)):
            group_diet_data = groups[j]
            group_info = {
                "group_name": group_diet_data["group_name"],
                "group_seq": int(group_diet_data["group_seq"]),
                "biomass": float(group_diet_data["biomass"]) if group_diet_data.get("biomass", np.nan) != "-9999" else np.nan,
                "pb": float(group_diet_data["pb"]) if group_diet_data.get("pb", np.nan) != "-9999" else np.nan,
                "qb": float(group_diet_data["qb"]) if group_diet_data.get("qb", np.nan) != "-9999" else np.nan,
                "ee": float(group_diet_data["ee"]) if group_diet_data.get("ee", np.nan) != "-9999" else np.nan,
                "M0b": float(group_diet_data["other_mort"]) if group_diet_data.get("other_mort") != "-9999" else np.nan,
                "ge": float(group_diet_data.get("ge", np.nan)) if group_diet_data.get("ge", np.nan) != "-9999" else np.nan,
                "gs": float(group_diet_data.get("gs", np.nan)) if group_diet_data.get("gs", np.nan) != "-9999" else np.nan,
                "respiration": float(group_diet_data.get("respiration", np.nan)) if group_diet_data.get("respiration", np.nan) != "-9999" else np.nan,
                "immigration": float(group_diet_data.get("immigration", np.nan)) if group_diet_data.get("immigration", np.nan) != "-9999" else np.nan,
                "emigration": float(group_diet_data.get("emigration", np.nan)) if group_diet_data.get("emigration", np.nan) != "-9999" else np.nan,
                "biomass_accum": float(group_diet_data.get("biomass_accum", np.nan)) if group_diet_data.get("biomass_accum", np.nan) != "-9999" else np.nan,
                "biomass_accum_rate": float(group_diet_data.get("biomass_accum_rate", np.nan)) if group_diet_data.get("biomass_accum_rate", np.nan) != "-9999" else np.nan,
                "export": float(group_diet_data.get("export", np.nan)) if group_diet_data.get("export", np.nan) != "-9999" else np.nan,
                "trophic_info": (lambda x: ["Regular", "PP", "DET"][int(x)])(float(group_diet_data["pp"])),
                "detritus_import": float(group_diet_data.get("detritus_import", 0)) if group_diet_data["pp"] == "2" else 0,
                "diet_import": float(group_diet_data.get("diet_imp", np.nan)) if group_diet_data.get("diet_imp", np.nan) != "-9999" else np.nan,
                "taxon_descr": group_diet_data.get("taxon_descr", None),
                "tl": np.nan,
            }
            
            # flow to det = flow from unassimilated food + flow from other mortality
            group_info['M0b'] = group_info['pb'] * (1-group_info['ee'])
            flow_from_food = group_info["biomass"] * group_info["qb"] * group_info["gs"]
            flow_from_bodies = group_info["biomass"] * group_info["M0b"]
            group_info['flow_to_det'] = flow_from_food + flow_from_bodies

            species_groups.append(SpeciesGroup.from_dict(group_info))

        return species_groups

    @classmethod
    def get_groups_df(cls, species_groups):
        """Convert species groups to DataFrame."""
        df_rows = [species_groups[i].to_df_row() for i in range(len(species_groups))]
        df = pd.concat(df_rows, ignore_index=True)
        df = df.replace("-9999", np.nan).replace(-9999, np.nan)

        df = df.rename(columns={  # change names
            'export': 'catch',
            })
        
        df['p'] = df['pb'] * df['biomass']  # production
        df['q'] = df['qb'] * df['biomass']  # consumption
        df['M0'] = df['p'] * (1-df['ee'])  # other mortality
        df['net_migration'] = df['emigration'] - df['immigration']  # net migration
        
        # production*EE = catch + predation + biomass_accum + net_migration:
        df['predation'] = df['p'] * df['ee'] - (df['catch'] + df['biomass_accum'] + df['net_migration'])
        df['egestion'] = df['q'] * df['gs']

        df['flow_to_det'] = df['egestion'] + df['M0']

        cols_to_return = ['group_name', 'trophic_info', 'taxon_descr', 'tl', 'ge', 'ee', 'catch', 
            'biomass', 'pb', 'qb', 'p', 'q', 'predation', 'M0', 'gs', 'egestion', 'respiration', 'biomass_accum', 'emigration', 'immigration', 'net_migration',
            'flow_to_det', 'detritus_import'
        ]

        groups_data = df.set_index('group_seq')[cols_to_return]
        groups_data = groups_data.sort_index(ascending=False)
        diet_import = df.set_index('group_seq')['diet_import'].sort_index(ascending=False)

        return groups_data, diet_import


# Legacy functions for old API
def get_seq2name(model_number, model_diet_datas):
    """Get mapping from group sequence number to group name (legacy function)."""
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
    """Get diet composition (DC) and detritus fate matrices (legacy function)."""
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


# Load Data (for legacy API):
try:
    model_metadatas = load_json_dict('real_models/model_metadatas')
    model_datas = load_json_dict('real_models/model_datas')
    model_diet_datas = load_json_dict('real_models/model_diet_datas')
    species_groups = read_json_SpeciesGroup_list("real_models/SpeciesGroups")
except Exception as e:
    # If loading fails, these will be empty and legacy API will fail gracefully
    model_metadatas = {}
    model_datas = {}
    model_diet_datas = {}
    species_groups = []


def get_model_metadata(model_id):
    return model_metadatas[str(model_id)]


def get_model_data(model_id):
    return model_datas[str(model_id)]


def get_model_diet_data(model_id):
    return model_diet_datas[str(model_id)]
