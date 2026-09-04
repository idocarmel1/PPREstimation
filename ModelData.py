from __future__ import annotations

from pathlib import Path
current_dir = Path(__file__).parent

import numpy as np
import pandas as pd
import json
import re
import warnings
from dataclasses import dataclass
from typing import Any, Optional
import warnings


def load_json_dict(filename: str) -> dict:
    """Load a JSON file (relative to this module's directory) and return it as a dict.

    The ``.json`` extension is appended automatically, and the path is resolved
    relative to the directory containing this module. This is the loader used for
    the bundled legacy data files (e.g. ``real_models/model_datas``).

    Args:
        filename (str): Path/name of the JSON file WITHOUT the ``.json`` suffix,
            relative to this module's directory (``current_dir``). For example,
            ``'real_models/model_metadatas'`` loads ``<current_dir>/real_models/model_metadatas.json``.

    Returns:
        dict: The parsed JSON content. (The top-level structure is whatever the
            file contains; for the bundled files it is a dict keyed by model id.)
    """
    # Resolve the file relative to this module's directory and append the extension,
    # so callers never have to know where the package lives on disk.
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
    def from_dict(cls, data: dict) -> "SpeciesGroupLegacy":
        """Build a ``SpeciesGroupLegacy`` from a flat dict of field values.

        Factory used when deserializing legacy JSON records. The dict keys are
        expected to match the dataclass field names exactly.

        Args:
            data (dict): Mapping whose keys are the dataclass field names
                (``group_name``, ``group_seq``, ``biomass``, ... ``taxons_included``).

        Returns:
            SpeciesGroupLegacy: A new instance populated from ``data``.
        """
        # Unpack the dict directly into the dataclass constructor; keys must match fields.
        return cls(**data)

    def __str__(self) -> str:
        """Return a human-readable one-line summary of this legacy species group.

        Includes the model identity, the group sequence/name, its trophic level,
        and the list of taxa (name + AphiaID) folded into the group. Intended for
        logging and quick inspection.

        Returns:
            str: A formatted description string.
        """
        return f"[{self.model_number}] '{self.model_name}' {self.model_country} ({self.model_year}) -> {self.group_seq}: {self.group_name} | tl: {self.tl}, taxons included: {[self.taxons_included[i]['taxon_name'] + ' (' + self.taxons_included[i]['AphiaID'] + ')' for i in range(len(self.taxons_included))]}"

    def to_df_row(self) -> pd.DataFrame:
        """Convert this legacy group into a single-row DataFrame of its scalar fields.

        The nested ``taxons_included`` list is dropped because it cannot be stored
        as a flat tabular value; everything else becomes one column.

        Returns:
            pd.DataFrame: A one-row DataFrame whose columns are the dataclass
                fields except ``taxons_included``.
        """
        # Copy so we do not mutate the dataclass's own __dict__.
        dct = self.__dict__.copy()
        # Drop the nested taxon list which is not a flat column value.
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
    def from_dict(cls, data: dict) -> "SpeciesGroup":
        """Build a ``SpeciesGroup`` from a flat dict of field values.

        Factory used by :meth:`ModelData.get_species_groups` after it has parsed
        and cleaned the raw JSON for one group. The dict keys must match the
        dataclass field names exactly.

        Args:
            data (dict): Mapping whose keys are the dataclass field names
                (``group_name``, ``group_seq``, ``biomass``, ... ``taxon_descr``).

        Returns:
            SpeciesGroup: A new instance populated from ``data``.
        """
        # Unpack the dict directly into the dataclass constructor; keys must match fields.
        return cls(**data)

    def to_df_row(self) -> pd.DataFrame:
        """Convert this group into a single-row DataFrame of all its fields.

        Unlike the legacy variant there is no nested taxon list to strip, so every
        field becomes a column.

        Returns:
            pd.DataFrame: A one-row DataFrame whose columns are all the dataclass fields.
        """
        # Copy so we do not mutate the dataclass's own __dict__, then wrap in a DataFrame.
        dct = self.__dict__.copy()
        return pd.DataFrame([dct])


def read_json_SpeciesGroup_list(filename: str) -> Any:
    """Load a JSON file and deserialize every object into a ``SpeciesGroupLegacy``.

    Uses ``json``'s ``object_hook`` so that EVERY JSON object encountered during
    parsing is fed through :meth:`SpeciesGroupLegacy.from_dict`. This is used to
    load the bundled legacy ``SpeciesGroups`` file as a flat list of dataclass
    instances. Note: because the hook fires on every nested object, the input is
    expected to be a list of records whose only objects are the group records.

    Args:
        filename (str): Path/name of the JSON file WITHOUT the ``.json`` suffix,
            relative to this module's directory.

    Returns:
        Any: Whatever the JSON top level decodes to after the hook runs — for the
            bundled data this is a ``list[SpeciesGroupLegacy]``.
    """
    # object_hook is invoked for each decoded JSON object; turn each into a dataclass.
    def decoder(dct: dict) -> SpeciesGroupLegacy:
        return SpeciesGroupLegacy.from_dict(dct)

    with open(f'{current_dir}/{filename}.json', 'r') as f:
        return json.load(f, object_hook=decoder)


class ModelData:
    """Container for one Ecopath model's group table, diet composition and detritus fate.

    Supports two construction APIs selected by the type of the constructor argument:
      * Legacy API: pass an ``int`` model number; data is pulled from the bundled
        module-level legacy structures (``species_groups``, ``model_diet_datas``).
      * New API: pass a ``str`` path to a per-model JSON file; data is parsed
        directly from that file.

    After construction the key attributes are:
      * ``groups_data`` (pd.DataFrame): per-group ecological parameters indexed by group_seq.
      * ``DC`` (pd.DataFrame): diet composition matrix (predator columns, prey rows).
      * ``det_fate`` (pd.DataFrame): detritus fate matrix.
      * ``seq2name`` / ``name2seq`` (dict): group sequence <-> name lookups.
      * ``model_number``, ``model_name``, ``model_country``, ``model_year``, ``lme``: metadata.
    """

    def __init__(self, model_input: int | str) -> None:
        """Initialize ModelData from either a model_number (int) or json_filepath (str).

        Dispatches to the correct loader based on the runtime type of the argument.

        Args:
            model_input (int | str): Selects the construction path:
                * ``int`` — a model number; loads from the bundled legacy data via
                  :meth:`_init_from_model_number`.
                * ``str`` — a path to a per-model JSON file; loads via
                  :meth:`_init_from_json_filepath`.

        Returns:
            None: The instance is populated in place.

        Raises:
            TypeError: If ``model_input`` is neither an ``int`` nor a ``str``.
        """
        if isinstance(model_input, int):
            # Old API: load from pre-existing JSON files by model number
            self._init_from_model_number(model_input)
        elif isinstance(model_input, str):
            # New API: load from a JSON filepath
            self._init_from_json_filepath(model_input)
        else:
            raise TypeError(f"model_input must be int or str, got {type(model_input)}")

    def _init_from_json_filepath(self, json_filepath: str) -> None:
        """Initialize ModelData from a per-model JSON file (new API).

        Parses model metadata out of the filename, loads the JSON, builds the group
        parameter table and the diet/detritus matrices, then injects a synthetic
        "diet_import" group (an extra sequence index that represents imported diet
        as if it were an external prey source) into all three structures so that
        downstream PPR computations can treat imports uniformly.

        Filename format: ``{first_number}_{model_number}_{model_name}_{model_years}.json``
        Example: ``227_227_Iceland_(1950).json`` or ``13_10013_Humboldt_Current_(1980).json``

        Args:
            json_filepath (str): Path to the per-model JSON file. Its stem must
                follow the ``{first_number}_{model_number}_{name}_({year})`` convention.

        Returns:
            None: Populates ``self`` attributes in place (``model_number``,
                ``model_name``, ``model_country``, ``model_year``, ``data_json``,
                ``groups_data``, ``lme``, ``seq2name``, ``name2seq``, ``DC``, ``det_fate``).
        """
        # Parse filepath
        filepath = Path(json_filepath)
        filename_no_ext = filepath.stem  # e.g., "227_Iceland_(1950)"

        # Extract model_number, model_name, and model_year from filename
        self.model_number, self.model_name, self.model_year = self._parse_filename(filename_no_ext)

        # Country is not available from the JSON, set to NaN
        self.model_country = np.nan

        # Load JSON data
        # (Type-only annotations below document the intended attribute types.)
        self.groups_data: pd.DataFrame
        self.groups_taxons: dict[int, list[dict]]  # group_seq: taxons_list

        # Read the raw model dict, extract the cleaned per-group records, then
        # fold them into a parameter DataFrame plus the per-group diet_import series.
        self.data_json = ModelData.load_json_dict_static(str(json_filepath))
        species_groups = ModelData.get_species_groups(self.data_json)
        self.groups_data, diet_import = ModelData.get_groups_df(species_groups)

        self.lme = 13  # Default, could be extracted from metadata if needed

        # Build the seq<->name lookups; get_seq2name already appends a 'diet_import' entry,
        # so the maximum sequence index is the synthetic import group's index.
        self.seq2name = ModelData.get_seq2name(self.data_json)
        self.name2seq = {v: k for k, v in self.seq2name.items()}
        import_seq = max(self.seq2name)

        # Add a new row for the synthetic import group and fill the few fields it needs;
        # everything else stays NaN because imports have no biomass/production of their own.
        self.groups_data.loc[import_seq] = np.nan
        self.groups_data.loc[import_seq, 'group_name'] = 'diet_import'
        self.groups_data.loc[import_seq, 'trophic_info'] = 'Import'
        self.groups_data.loc[import_seq, 'tl'] = 1.0
        self.groups_data.loc[import_seq, 'respiration'] = 0.0
        # The import group consumes nothing, so its own diet_import contribution is 0.
        diet_import.loc[import_seq] = 0

        # DC and Detritus fate:
        DC, det_fate = ModelData.get_DC(self.data_json)
        DC.loc[import_seq] = 0  # add import data to DC  (import group eats nothing -> zero row)
        DC[import_seq] = diet_import  # add import data to DC  (import as a prey column = per-group diet_import)
        # Sort both axes descending so the matrix orientation matches downstream expectations.
        self.DC = DC.sort_index(ascending=False).sort_index(axis=1, ascending=False)

        det_fate.loc[import_seq] = 0  # add diet_import to det_fate (zero row for the import group)
        det_seq = self.groups_data.index[self.groups_data['trophic_info'] == 'DET']
        # Normalize orientation, force the detritus self-identity, and default a fully-missing
        # (degenerate) single-DET matrix to the closed system. See _finalize_det_fate.
        self.det_fate = ModelData._finalize_det_fate(det_fate, self.groups_data)

        # fix export of detritus groups to be 0:
        self.groups_data.loc[det_seq, 'catch'] = 0

    def _init_from_model_number(self, model_number: int) -> None:
        """Initialize ModelData from a model number using the bundled legacy data.

        Scans the module-level ``species_groups`` list for all groups belonging to
        ``model_number``, derives flow quantities (production, consumption, predation,
        egestion, flow to detritus, etc.) from the raw Ecopath parameters, and builds
        the group table, diet composition and detritus fate matrices. As in the new
        API, a synthetic "diet_import" group is appended to every structure.

        Args:
            model_number (int): The model identifier to load. Must match the
                ``model_number`` attribute on at least one entry of the bundled
                ``species_groups`` list.

        Returns:
            None: Populates ``self`` attributes in place (``groups_data``,
                ``groups_taxons``, ``model_number``, ``model_name``,
                ``model_country``, ``model_year``, ``lme``, ``seq2name``,
                ``name2seq``, ``DC``, ``det_fate``).

        Raises:
            ValueError: If no species group matches ``model_number``.
        """
        def _groups_data(model_number: int) -> tuple[pd.DataFrame, dict[int, list[dict]], int, pd.Series]:
            """Build the legacy group table, taxon map, example index and diet_import series.

            Inner helper that filters the bundled ``species_groups`` for this model,
            renames legacy columns to the canonical names, and computes the derived
            flow columns.

            Args:
                model_number (int): Model id to filter on.

            Returns:
                tuple[pd.DataFrame, dict[int, list[dict]], int, pd.Series]:
                    * groups_data — per-group parameters indexed by group_seq (descending).
                    * taxons — mapping group_seq -> list of taxon dicts.
                    * model_index — index into ``species_groups`` of one matching row
                      (used later to read shared model metadata).
                    * diet_import — per-group imported-diet proportions indexed by group_seq.

            Raises:
                ValueError: If no group matches ``model_number``.
            """
            def to_df_row(c: SpeciesGroupLegacy) -> pd.DataFrame:
                """Turn one legacy group into a one-row DataFrame, dropping model-level fields.

                The model-identity and taxon fields are removed because they are
                constant across the model (read separately) or non-tabular.

                Args:
                    c (SpeciesGroupLegacy): The group to convert.

                Returns:
                    pd.DataFrame: One-row DataFrame of the group's per-group fields.
                """
                dct = c.__dict__.copy()
                # Strip out non-tabular / model-level fields so only per-group columns remain.
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
            # Collect every group row that belongs to the requested model number.
            for i in range(len(species_groups)):
                if hasattr(species_groups[i], 'model_number') and species_groups[i].model_number == model_number:
                    model_index = i
                    row = to_df_row(species_groups[i])
                    df_rows.append(row)
                    # Keep the taxon list keyed by group sequence for later lookup.
                    if hasattr(species_groups[i], 'taxons_included'):
                        taxons[species_groups[i].group_seq] = species_groups[i].taxons_included

            if not df_rows:
                raise ValueError(f"No species groups found for model {model_number}")

            # Stack all the one-row frames into the full group table.
            df = pd.concat(df_rows, ignore_index=True)
            # The sentinel -9999 (string or numeric) means "missing"; convert to NaN.
            df = df.replace("-9999", np.nan).replace(-9999, np.nan)

            df = df.rename(columns={  # change names
                'export': 'catch',
                'prop_unassimilated_food': 'gs',
                'gross_efficiency': 'ge'
                }, errors='ignore')

            # Derived absolute flows from the Ecopath rate parameters:
            df['p'] = df['pb'] * df['biomass']  # production
            df['q'] = df['qb'] * df['biomass']  # consumption
            df['M0'] = df['p'] * (1-df['ee'])  # other mortality
            df['net_migration'] = df['emigration'] - df['immigration']  # net migration

            # production*EE = catch + predation + biomass_accum + net_migration:
            # solve that mass balance for predation (the unknown sink).
            df['predation'] = df['p'] * df['ee'] - (df['catch'] + df['biomass_accum'] + df['net_migration'])
            df['egestion'] = df['q'] * df['gs']

            # Flow to detritus = unassimilated food (egestion) + non-predation mortality (M0).
            df['flow_to_det'] = df['egestion'] + df['M0']

            cols_to_return = ['group_name', 'trophic_info', 'tl', 'ge', 'ee', 'catch',
                'biomass', 'pb', 'qb', 'p', 'q', 'predation', 'M0', 'gs', 'egestion', 'respiration', 'biomass_accum', 'emigration', 'immigration', 'net_migration',
                'flow_to_det', 'detritus_import'
            ]

            # Index by group sequence and keep only the canonical column subset, sorted descending.
            groups_data = df.set_index('group_seq')[cols_to_return]
            groups_data = groups_data.sort_index(ascending=False)
            # diet_import is carried separately as a Series (one value per group).
            diet_import = df.set_index('group_seq')['diet_import'].sort_index(ascending=False)

            return groups_data, taxons, model_index, diet_import

        # (Type-only annotations below document the intended attribute types.)
        self.groups_data: pd.DataFrame
        self.groups_taxons: dict[int, list[dict]]  # group_seq: taxons_list

        # Build the core table plus the example-row index and the diet_import series.
        self.groups_data, self.groups_taxons, model_index, diet_import = _groups_data(model_number)
        # Any matching row carries the shared model-level metadata; use the first one.
        model_species_group_example = species_groups[model_index]

        self.model_number: int = model_number
        self.model_name: str = model_species_group_example.model_name
        self.model_country: str = model_species_group_example.model_country
        self.model_year: str = model_species_group_example.model_year
        self.lme: int = model_species_group_example.lme

        # Build seq<->name lookups (legacy helper already appends the 'diet_import' entry).
        self.seq2name = get_seq2name(model_number, model_diet_datas)
        self.name2seq = {v: k for k, v in self.seq2name.items()}
        import_seq = max(self.seq2name)

        # Append the synthetic import group row, mirroring the new-API path.
        self.groups_data.loc[import_seq] = np.nan
        self.groups_data.loc[import_seq, 'group_name'] = 'diet_import'
        self.groups_data.loc[import_seq, 'trophic_info'] = 'Import'
        self.groups_data.loc[import_seq, 'tl'] = 1.0
        self.groups_data.loc[import_seq, 'respiration'] = 0.0
        diet_import.loc[import_seq] = 0

        # DC and Detritus fate:
        DC, det_fate = get_DC(model_number, model_diet_datas)
        DC.loc[import_seq] = 0  # add import data to DC  (import group eats nothing -> zero row)
        DC[import_seq] = diet_import  # add import data to DC  (import as a prey column)
        self.DC = DC.sort_index(ascending=False).sort_index(axis=1, ascending=False)

        det_fate.loc[import_seq] = 0  # add diet_import to det_fate (zero row)
        det_fate[import_seq] = 0  # add diet_import to det_fate (zero column)
        # Normalize orientation, force the detritus self-identity, and default a fully-missing
        # (degenerate) single-DET matrix to the closed system. See _finalize_det_fate.
        self.det_fate = ModelData._finalize_det_fate(det_fate, self.groups_data)

    @staticmethod
    def _finalize_det_fate(det_fate: pd.DataFrame, groups_data: pd.DataFrame, tol: float = 1e-9) -> pd.DataFrame:
        """Normalize a raw detritus-fate matrix into a clean groups x detritus-pool matrix.

        Takes the per-prey ``detritus_fate`` matrix assembled during loading and turns it into
        the canonical detritus-fate matrix the rest of the pipeline expects: rows indexed by
        every group, columns restricted to the detritus (DET) groups, with entry ``[g, d]`` the
        fraction of group ``g``'s ``flow_to_det`` (= M0 + egestion) routed to detritus pool ``d``.

        The steps are:

        1. Restrict the columns to the DET groups and 0-fill, dropping any spurious non-detritus
           columns that survived the raw per-prey construction.
        2. Force the detritus<->detritus sub-block to the identity, so each detritus pool maps to
           itself by default.
        3. Zero the Import group's row (the import pseudo-group produces no detritus; its
           ``flow_to_det`` is 0, so this has no numerical effect but keeps the matrix clean).
        4. Detect a *fully degenerate* matrix -- one where no living (non-import, non-detritus)
           group carries any routing at all, which happens when the source data has no detritus
           fate (e.g. reconstructed models whose raw ``detritus_fate`` is all zero). For a
           single-DET model this is defaulted to the closed system (every living group routes
           1.0 to the sole pool); a multi-DET degenerate matrix is left untouched so that
           :meth:`validate_det_fate` raises (the per-pool split cannot be inferred).

        Partially populated matrices (some rows summing to <1) are left exactly as given: a row
        summing to less than 1 is treated as legitimate export out of the system and is only
        flagged, not altered, by :meth:`validate_det_fate`.

        Args:
            det_fate (pd.DataFrame): the raw detritus-fate matrix (rows = groups, columns =
                whatever prey carried a ``detritus_fate`` value).
            groups_data (pd.DataFrame): the per-group table, used for ``trophic_info`` to
                identify the DET and Import groups.
            tol (float): magnitude below which a row sum is treated as zero when testing for the
                fully-degenerate case. Defaults to 1e-9.

        Returns:
            pd.DataFrame: the cleaned detritus-fate matrix (groups x DET pools), sorted by
            descending seq on both axes.
        """
        det_seq = list(groups_data.index[groups_data['trophic_info'] == 'DET'])
        import_seq = list(groups_data.index[groups_data['trophic_info'] == 'Import'])

        # Restrict to detritus-pool columns over the full group index, 0-filling everything else.
        det_fate = det_fate.reindex(index=groups_data.index, columns=det_seq).fillna(0.0)

        # Force the detritus<->detritus identity block.
        det_fate.loc[det_seq, :] = 0.0
        if len(det_seq) > 0:
            det_fate.loc[det_seq, det_seq] = np.eye(len(det_seq))

        # The import pseudo-group routes nothing to detritus.
        det_fate.loc[import_seq, :] = 0.0

        # Degenerate (no routing data at all) -> default single-DET to the closed system.
        living = [s for s in groups_data.index if s not in set(det_seq) and s not in set(import_seq)]
        if living and (det_fate.loc[living].sum(axis=1) <= tol).all():
            if len(det_seq) == 1:
                det_fate.loc[living, det_seq[0]] = 1.0
            # multi-DET degenerate: leave as-is -> validate_det_fate raises.

        return det_fate.sort_index(ascending=False).sort_index(axis=1, ascending=False)

    @staticmethod
    def validate_det_fate(det_fate: pd.DataFrame, groups_data: pd.DataFrame, tol: float = 1e-3) -> None:
        """Validate a finalized detritus-fate matrix, raising on inconsistencies.

        Checks that the detritus-fate matrix is usable for an open-system detritus model, where
        each group routes its ``flow_to_det`` (= M0 + egestion) to the detritus pools and any
        shortfall (row sum < 1) is treated as export out of the system:

        - **raise** ``ValueError`` if any non-import row sums to more than ``1 + tol`` (an
          impossible over-allocation -- a group cannot send more than 100% of its dead matter to
          detritus);
        - **raise** ``ValueError`` if the matrix is fully degenerate (no living group routes any
          flow to detritus), which for a multi-DET model means the per-pool split is unknown and
          cannot be defaulted (single-DET degeneracy is already repaired by
          :meth:`_finalize_det_fate`);
        - **warn** (``RuntimeWarning``, listing the offending groups) if any non-import row sums
          to less than ``1 - tol`` -- this is allowed (the remainder is exported) but flagged so
          the export is never silent.

        Args:
            det_fate (pd.DataFrame): a finalized detritus-fate matrix (groups x DET pools), e.g.
                from :meth:`_finalize_det_fate`.
            groups_data (pd.DataFrame): the per-group table, used for ``trophic_info`` and group
                names in the messages.
            tol (float): tolerance for the row-sum comparisons. Defaults to 1e-3.

        Returns:
            None.

        Raises:
            ValueError: on over-allocated rows or a fully-degenerate multi-DET matrix.
        """
        det_seq = set(groups_data.index[groups_data['trophic_info'] == 'DET'])
        import_seq = set(groups_data.index[groups_data['trophic_info'] == 'Import'])
        names = groups_data['group_name'] if 'group_name' in groups_data.columns else pd.Series(dtype=object)

        non_import = [s for s in det_fate.index if s not in import_seq]
        rowsums = det_fate.loc[non_import].sum(axis=1)

        def _label(seq_list):
            return ", ".join(f"{s} ({names.get(s, '?')})" for s in seq_list)

        # Over-allocation: physically impossible.
        over = [s for s in non_import if rowsums[s] > 1 + tol]
        if over:
            raise ValueError(
                "det_fate has rows summing to more than 1 (over-allocated detritus routing) for "
                f"group(s): {_label(over)}. Each group can route at most 100% of its flow_to_det."
            )

        # Fully degenerate: no living group routes anything (multi-DET missing data).
        living = [s for s in non_import if s not in det_seq]
        if living and (rowsums.loc[living] <= tol).all():
            raise ValueError(
                "det_fate carries no detritus routing for any living group. For a multi-DET model "
                "the per-pool split cannot be inferred -- supply DetritusFate data for "
                f"model with detritus groups {sorted(det_seq)}."
            )

        # Partial rows: legitimate export, but flag.
        under = [s for s in non_import if rowsums[s] < 1 - tol]
        if under:
            warnings.warn(
                f"{len(under)} group(s) route less than 100% of their flow_to_det to detritus; "
                f"the remainder is treated as export out of the system: {_label(under)}.",
                RuntimeWarning,
            )

    @staticmethod
    def validate_DC(DC: pd.DataFrame, groups_data: pd.DataFrame, normalize: bool = True, tol: float = 1e-3) -> None:
        """Validate that every consumer's diet composition sums to 1, raising otherwise.

        Each feeding (Regular) group's diet -- including its imported-diet column -- must sum to
        1 by the Ecopath definition of a diet composition. Non-feeders (primary producers,
        detritus and the import pseudo-group) have all-zero diet rows and are exempt.

        Args:
            DC (pd.DataFrame): the diet-composition matrix (rows = predators, columns = prey,
                including the imported-diet column).
            groups_data (pd.DataFrame): the per-group table, used for ``trophic_info`` (to pick
                the consumer rows) and group names in the message.
            normalize (bool): wether to normalize the DC rows to 1.
            tol (float): allowed absolute deviation of a consumer row sum from 1; chosen loose
                enough (1e-3) to absorb the rounding present in published diet matrices.
                Defaults to 1e-3.

        Returns:
            None.

        Raises:
            ValueError: if any consumer (Regular) row sums to a value more than ``tol`` away
                from 1.
        """
        is_regular = groups_data['trophic_info'] == 'Regular'
        consumer_seq = [s for s in DC.index if s in set(groups_data.index[is_regular])]
        names = groups_data['group_name'] if 'group_name' in groups_data.columns else pd.Series(dtype=object)

        rowsums = DC.loc[consumer_seq].sum(axis=1)
        bad = [s for s in consumer_seq if abs(rowsums[s] - 1.0) > tol]
        if bad and not normalize:
            detail = ", ".join(f"{s} ({names.get(s, '?')})={rowsums[s]:.4f}" for s in bad)
            raise ValueError(
                f"{len(bad)} consumer group(s) have a diet composition (including diet_import) "
                f"that does not sum to 1 (tol={tol}): {detail}."
            )
        elif bad:
            detail = ", ".join(f"{s} ({names.get(s, '?')})={rowsums[s]:.4f}" for s in bad)
            # Issue a standard user warning
            warnings.warn(f"{len(bad)} consumer group(s) have a diet composition (including diet_import) "
                            f"that does not sum to 1 (tol={tol}): {detail}.", RuntimeWarning)
        if normalize:
            DC = DC.div(DC.sum(axis=1), axis=0).fillna(0)
        return DC

    @staticmethod
    def _parse_filename(filename_no_ext: str) -> tuple[int, str, str]:
        """Parse a model filename stem into (model_number, model_name, model_year).

        Expects the convention ``{first_number}_{model_number}_{model_name}_({model_year})``,
        i.e. two leading numeric tokens. The year is taken from the parentheses at the
        end; the first token is a grouping/source id and is discarded; the *second*
        token is the model id; everything after it is the name.

        The second token is used as the model number so a single parser handles both
        file families uniformly:
          * Multi-model source files carry a distinct model id in the second token,
            e.g. ``"13_10013_Humboldt_Current_(1980)"`` -> ``(10013, "Humboldt_Current", "1980")``.
          * Single-model files repeat their number, so the second token reproduces the
            original id, e.g. ``"227_227_Iceland_(1950)"`` -> ``(227, "Iceland", "1950")``.

        Args:
            filename_no_ext (str): The filename stem (no directory, no ``.json``).

        Returns:
            tuple[int, str, str]: ``(model_number, model_name, model_year)`` where
                ``model_number`` is the second numeric token and the other two are strings.

        Raises:
            ValueError: If the trailing ``(year)`` is missing, or the number/name
                portion cannot be split into the two numbers plus a name.
        """
        # Extract year from parentheses at the end
        year_match = re.search(r'\(([^)]+)\)$', filename_no_ext)
        if not year_match:
            raise ValueError(f"Could not extract year from filename: {filename_no_ext}")

        model_year = year_match.group(1)

        # Remove the year part to get the rest (and trim a trailing underscore separator).
        rest = filename_no_ext[:year_match.start()].rstrip('_')

        # Split off the first two underscore-separated tokens (the two numbers); the
        # remainder is the model name. Format: {first_number}_{model_number}_{name}.
        parts = rest.split('_', 2)  # Split on the first two underscores only
        if len(parts) != 3:
            raise ValueError(f"Could not parse two numbers and a model_name from: {rest}")

        # parts[0] is the grouping/source id (discarded); parts[1] is the model number.
        model_number = int(parts[1])
        model_name = parts[2]

        return model_number, model_name, model_year

    @classmethod
    def load_json_dict_static(cls, filepath: str) -> dict:
        """Load a JSON file from an absolute/explicit path and return it as a dict.

        Unlike the module-level :func:`load_json_dict`, this takes a full path
        (extension included) and does not prepend the module directory.

        Args:
            filepath (str): Full path to the ``.json`` file to read.

        Returns:
            dict: The parsed JSON content (for model files, the model dict).
        """
        with open(filepath, 'r') as json_file:
            return json.load(json_file)

    @classmethod
    def get_seq2name(cls, model: dict) -> Optional[dict]:
        """Build a {group_seq: group_name} mapping from a model dict, plus a diet_import entry.

        Reads the ``group`` list out of the model JSON and maps each group's integer
        sequence to its name, then appends one extra synthetic entry (max seq + 1)
        named ``'diet_import'`` to represent imported diet.

        Args:
            model (dict): The model JSON dict. Must contain a ``'group'`` list whose
                items have ``'group_seq'`` and ``'group_name'`` keys.

        Returns:
            Optional[dict]: A ``{int: str}`` mapping including the synthetic
                ``diet_import`` entry, or ``None`` if ``model`` is not a dict or has
                no groups.
        """
        # Guard: only dicts can carry the expected structure.
        if not isinstance(model, dict):
            return None

        groups = model.get('group', {})
        if len(groups) == 0:
            return None

        # Map each group's integer sequence to its name.
        seq2name = {int(g["group_seq"]): g["group_name"] for g in groups}
        # Reserve the next index above the maximum for the synthetic import group.
        import_seq = max(seq2name) + 1
        seq2name[import_seq] = 'diet_import'

        return seq2name

    @classmethod
    def get_DC(cls, model_json_dict: dict) -> Optional[tuple[pd.DataFrame, pd.DataFrame]]:
        """Build the diet composition (DC) and detritus fate matrices from a model dict.

        For each group, reads its diet description and produces two square matrices
        indexed by group sequence: ``DC`` (proportion of each prey in a predator's
        diet) and ``detritus_fate`` (the detritus-routing fraction associated with
        each diet entry). Groups with no diet get a single placeholder cell that
        becomes 0 after the ``fillna``.

        Args:
            model_json_dict (dict): The model JSON dict, expected to contain a
                ``'group'`` list where each group may carry a
                ``diet_descr.diet`` list of ``{prey_seq, proportion, detritus_fate}``.

        Returns:
            Optional[tuple[pd.DataFrame, pd.DataFrame]]: ``(DC, detritus_fate)`` —
                two square DataFrames indexed/columned by group sequence and sorted
                ascending on both axes. Returns ``None`` if the input is not a dict
                or has no groups.
        """
        # Accept only dicts; otherwise there is nothing to parse.
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
            # Pull the diet list for this group, tolerating missing/empty diet_descr.
            diet_descr = g.get('diet_descr', {})
            diet_descr = diet_descr if diet_descr else {}
            diet = diet_descr.get('diet', None)
            if not diet:
                # No diet: seed a single self-referential placeholder cell (becomes 0 after fillna).
                DC_dict[int(g['group_seq'])] = {int(g['group_seq']): None}
                detritus_fate_dict[int(g['group_seq'])] = {int(g['group_seq']): None}
            else:
                # The diet may be a single dict or a list of dicts; normalize to a list.
                diet = diet if isinstance(diet, list) else [diet]
                # Map prey_seq -> proportion and prey_seq -> detritus_fate for this predator.
                # (prey_seq can arrive as a float-like string, hence int(float(...))).
                DC_dict[int(g['group_seq'])] = {int(float(d['prey_seq'])): float(d['proportion'].replace('-9999', '0')) for d in diet}
                detritus_fate_dict[int(g['group_seq'])] = {int(float(d['prey_seq'])): float(d['detritus_fate'].replace('-9999', '0')) for d in diet}

        # add missing columns:
        # Build the matrix from {predator: {prey: value}} with predators as the row index.
        DC = pd.DataFrame.from_dict(DC_dict, orient='index').fillna(0)
        # Ensure every group also exists as a column (prey that nobody ate would be absent).
        for g in groups:
            if int(g['group_seq']) not in DC.columns:
                DC[int(g['group_seq'])] = 0

        # Sort both axes so rows and columns share the same ascending group order.
        DC = DC.sort_index().sort_index(axis=1)

        detritus_fate = pd.DataFrame.from_dict(detritus_fate_dict, orient='index').fillna(0)
        detritus_fate = detritus_fate.sort_index().sort_index(axis=1)

        return DC, detritus_fate

    @classmethod
    def get_species_groups(cls, model_json: dict) -> list["SpeciesGroup"]:
        """Parse and clean the raw group records of a model into ``SpeciesGroup`` objects.

        Iterates the JSON ``group`` list, coercing each numeric field to float while
        mapping the ``-9999`` missing sentinel to NaN, decoding the Ecopath ``pp``
        flag into a trophic category, and computing the flow-to-detritus quantity
        from egestion plus non-predation mortality.

        Args:
            model_json (dict): The model JSON dict; expected to contain a ``'group'``
                list of raw group records.

        Returns:
            list[SpeciesGroup]: One cleaned dataclass per group (empty list if the
                model has no groups).

        Notes:
            The ``pp`` flag is decoded as: ``0`` -> ``"Regular"`` (consumer),
            ``1`` -> ``"PP"`` (primary producer), ``2`` -> ``"DET"`` (detritus).
            ``detritus_import`` is only read for detritus groups (``pp == "2"``),
            otherwise forced to 0.
        """
        groups = model_json.get('group', {})
        if len(groups) == 0:
            return []
        species_groups = []
        for j in range(len(groups)):
            group_diet_data = groups[j]
            # Coerce each field to float, treating the "-9999" sentinel as missing (NaN).
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
                # Decode Ecopath pp flag: 0->Regular consumer, 1->primary producer, 2->detritus.
                "trophic_info": (lambda x: ["Regular", "PP", "DET"][int(x)])(float(group_diet_data["pp"])),
                # detritus_import only meaningful for detritus groups (pp == "2"); else 0.
                "detritus_import": float(group_diet_data.get("detritus_import", 0)) if group_diet_data["pp"] == "2" else 0,
                "diet_import": float(group_diet_data.get("diet_imp", np.nan)) if group_diet_data.get("diet_imp", np.nan) != "-9999" else np.nan,
                "taxon_descr": group_diet_data.get("taxon_descr", None),
                "tl": np.nan,
            }

            # flow to det = flow from unassimilated food + flow from other mortality
            # Recompute M0b from pb and ee so it is consistent with the other fields.
            group_info['M0b'] = group_info['pb'] * (1-group_info['ee'])
            # Unassimilated food flow = biomass * consumption rate * unassimilated fraction.
            flow_from_food = group_info["biomass"] * group_info["qb"] * group_info["gs"]
            # Body/other-mortality flow = biomass * non-predation mortality rate.
            flow_from_bodies = group_info["biomass"] * group_info["M0b"]
            group_info['flow_to_det'] = flow_from_food + flow_from_bodies

            species_groups.append(SpeciesGroup.from_dict(group_info))

        return species_groups

    @classmethod
    def get_groups_df(cls, species_groups: list["SpeciesGroup"]) -> tuple[pd.DataFrame, pd.Series]:
        """Assemble a per-group parameter table and diet_import series from SpeciesGroup objects.

        Concatenates the one-row frames of each group, normalizes the ``-9999``
        sentinel to NaN, renames ``export`` to ``catch``, and derives the absolute
        flow columns (production, consumption, other mortality, net migration,
        predation, egestion, flow to detritus). Finally indexes by group sequence.

        Args:
            species_groups (list[SpeciesGroup]): The cleaned group objects, typically
                produced by :meth:`get_species_groups`.

        Returns:
            tuple[pd.DataFrame, pd.Series]:
                * groups_data — DataFrame of the canonical column subset, indexed by
                  ``group_seq`` and sorted descending.
                * diet_import — Series of per-group imported-diet proportions, indexed
                  by ``group_seq`` and sorted descending.
        """
        # Turn each SpeciesGroup into a one-row frame and stack them.
        df_rows = [species_groups[i].to_df_row() for i in range(len(species_groups))]
        df = pd.concat(df_rows, ignore_index=True)
        # Normalize the -9999 missing sentinel (string or numeric) to NaN.
        df = df.replace("-9999", np.nan).replace(-9999, np.nan)

        df = df.rename(columns={  # change names
            'export': 'catch',
            })

        # Derived absolute flows from the Ecopath rate parameters:
        df['p'] = df['pb'] * df['biomass']  # production
        df['q'] = df['qb'] * df['biomass']  # consumption
        df['M0'] = df['p'] * (1-df['ee'])  # other mortality
        df['net_migration'] = df['emigration'] - df['immigration']  # net migration

        # production*EE = catch + predation + biomass_accum + net_migration:
        # solve the mass-balance for predation (the unknown sink).
        df['predation'] = df['p'] * df['ee'] - (df['catch'] + df['biomass_accum'] + df['net_migration'])
        df['egestion'] = df['q'] * df['gs']

        # Flow to detritus = unassimilated food (egestion) + non-predation mortality (M0).
        df['flow_to_det'] = df['egestion'] + df['M0']

        cols_to_return = ['group_name', 'trophic_info', 'taxon_descr', 'tl', 'ge', 'ee', 'catch',
            'biomass', 'pb', 'qb', 'p', 'q', 'predation', 'M0', 'gs', 'egestion', 'respiration', 'biomass_accum', 'emigration', 'immigration', 'net_migration',
            'flow_to_det', 'detritus_import'
        ]

        # Index by group sequence, keep the canonical columns, sort descending.
        groups_data = df.set_index('group_seq')[cols_to_return]
        groups_data = groups_data.sort_index(ascending=False)
        # diet_import carried separately as a Series.
        diet_import = df.set_index('group_seq')['diet_import'].sort_index(ascending=False)

        return groups_data, diet_import


# Legacy functions for old API
def get_seq2name(model_number: int, model_diet_datas: dict) -> Optional[dict]:
    """Build a {group_seq: group_name} mapping for a model from the bundled legacy data.

    Legacy counterpart of :meth:`ModelData.get_seq2name`: it first looks the model
    up by id inside ``model_diet_datas`` and then maps each group's integer
    sequence to its name, appending a synthetic ``'diet_import'`` entry.

    Args:
        model_number (int): The model id; used (as a string key) to index
            ``model_diet_datas``.
        model_diet_datas (dict): The bundled mapping ``{str(model_id): model_dict}``.

    Returns:
        Optional[dict]: A ``{int: str}`` mapping including the synthetic
            ``diet_import`` entry, or ``None`` if the looked-up model is not a dict
            or has no groups.
    """
    # Look the model up by its string id in the bundled legacy diet data.
    model = model_diet_datas[str(model_number)]
    if not isinstance(model, dict):
        return None

    groups = model.get('group', {})
    if len(groups) == 0:
        return None

    # Map each group's integer sequence to its name, then add the synthetic import group.
    seq2name = {int(g["group_seq"]): g["group_name"] for g in groups}
    import_seq = max(seq2name) + 1
    seq2name[import_seq] = 'diet_import'

    return seq2name


def get_DC(model_number: int, model_diet_datas: dict) -> Optional[tuple[pd.DataFrame, pd.DataFrame]]:
    """Build the diet composition (DC) and detritus fate matrices from the bundled legacy data.

    Legacy counterpart of :meth:`ModelData.get_DC`: looks the model up by id in
    ``model_diet_datas`` and then constructs the two square matrices indexed by
    group sequence exactly as the class method does (here ``prey_seq`` is parsed
    directly as an int rather than via float).

    Args:
        model_number (int): The model id; used (as a string key) to index
            ``model_diet_datas``.
        model_diet_datas (dict): The bundled mapping ``{str(model_id): model_dict}``.

    Returns:
        Optional[tuple[pd.DataFrame, pd.DataFrame]]: ``(DC, detritus_fate)`` — two
            square DataFrames indexed/columned by group sequence and sorted ascending
            on both axes. Returns ``None`` if the looked-up model is not a dict or has
            no groups.
    """
    # Look the model up by its string id in the bundled legacy diet data.
    model = model_diet_datas[str(model_number)]
    if not isinstance(model, dict):
        return None

    groups = model.get('group', {})
    if len(groups) == 0:
        return None

    DC_dict = {}
    detritus_fate_dict = {}

    for g in groups:
        # Pull the diet list for this group, tolerating missing/empty diet_descr.
        diet_descr = g.get('diet_descr', {})
        diet_descr = diet_descr if diet_descr else {}
        diet = diet_descr.get('diet', None)
        if not diet:
            # No diet: seed a single self-referential placeholder cell (becomes 0 after fillna).
            DC_dict[int(g['group_seq'])] = {int(g['group_seq']): None}
            detritus_fate_dict[int(g['group_seq'])] = {int(g['group_seq']): None}
        else:
            # Normalize a single diet dict to a list, then map prey_seq -> proportion / detritus_fate.
            diet = diet if isinstance(diet, list) else [diet]
            DC_dict[int(g['group_seq'])] = {int(d['prey_seq']): float(d['proportion']) for d in diet}
            detritus_fate_dict[int(g['group_seq'])] = {int(d['prey_seq']): float(d['detritus_fate']) for d in diet}

    # add missing columns:
    # Build the matrix with predators as the row index, then ensure every group is also a column.
    DC = pd.DataFrame.from_dict(DC_dict, orient='index').fillna(0)
    for g in groups:
        if int(g['group_seq']) not in DC.columns:
            DC[int(g['group_seq'])] = 0

    # Sort both axes so rows and columns share the same ascending group order.
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


def get_model_metadata(model_id: int | str) -> Any:
    """Return the bundled metadata record for a model id.

    Thin accessor over the module-level ``model_metadatas`` mapping.

    Args:
        model_id (int | str): The model id; coerced to ``str`` for the dict lookup.

    Returns:
        Any: The metadata value stored for that model (typically a dict).

    Raises:
        KeyError: If no metadata exists for ``model_id``.
    """
    return model_metadatas[str(model_id)]


def get_model_data(model_id: int | str) -> Any:
    """Return the bundled core data record for a model id.

    Thin accessor over the module-level ``model_datas`` mapping.

    Args:
        model_id (int | str): The model id; coerced to ``str`` for the dict lookup.

    Returns:
        Any: The data value stored for that model (typically a dict).

    Raises:
        KeyError: If no data exists for ``model_id``.
    """
    return model_datas[str(model_id)]


def get_model_diet_data(model_id: int | str) -> Any:
    """Return the bundled diet data record for a model id.

    Thin accessor over the module-level ``model_diet_datas`` mapping.

    Args:
        model_id (int | str): The model id; coerced to ``str`` for the dict lookup.

    Returns:
        Any: The diet data value stored for that model (typically a dict).

    Raises:
        KeyError: If no diet data exists for ``model_id``.
    """
    return model_diet_datas[str(model_id)]
