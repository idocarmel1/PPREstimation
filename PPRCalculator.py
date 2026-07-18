from __future__ import annotations

import numpy as np
import pandas as pd
import sympy as sm
from tqdm.notebook import tqdm
from scipy.stats import gamma
import igraph as ig
from scipy.optimize import minimize
import warnings
from typing import Any, Optional

from ModelData import ModelData
from utils import mat_from_np, remove_cycles, move_scattered_identity
from copy import deepcopy

class PPRCalculator:

    # constructors:
    def __init__(self, model_number: int | str, underdetermined: bool = False, zero_catch: bool = True, zero_biomass_accum: bool = True, default_gs: bool = True, weight_flow: float = 1.0, weight_guess: float = 1.0) -> None:
        """Primary constructor: build the calculator directly from a model identifier.

        Loads a ModelData from the given model number / json filepath, delegates to
        from_modeldata to run the full property-filling pipeline, then copies the
        resulting instance's attributes onto self so PPRCalculator(path) returns a
        ready-to-use object.

        Args:
            model_number (int | str): model number or json filepath understood by ModelData.
            underdetermined (bool): if True, missing mass-balance variables are filled by
                the Linear Inverse Model (apply_lim) instead of the deterministic Ecopath
                defaults alone. Defaults to False.
            zero_catch (bool): if True, treat missing catch values as 0. Defaults to True.
            zero_biomass_accum (bool): if True, treat missing biomass-accumulation values
                as 0. Defaults to True.
            default_gs (bool): if True, assign the default GS=0.2 to regular groups lacking
                a gs value. Defaults to True.
            weight_flow (float): weight of the minimum-total-flow (parsimony) penalty in the
                LIM objective. Defaults to 1.0.
            weight_guess (float): weight of the deviation-from-biological-guess penalty in the
                LIM objective. Defaults to 1.0.

        Returns:
            None. Populates self in place.
        """
        # load model:
        instance = self.from_modeldata(
            ModelData(model_number),
            underdetermined=underdetermined,
            zero_catch=zero_catch, 
            zero_biomass_accum=zero_biomass_accum, 
            default_gs=default_gs, 
            weight_flow=weight_flow, 
            weight_guess=weight_guess
        )
        self.__dict__.update(instance.__dict__)
         
    @classmethod
    def from_dict(cls, data_dict: dict, underdetermined: bool = False, zero_catch: bool = True, zero_biomass_accum: bool = True, default_gs: bool = True, weight_flow: float = 1.0, weight_guess: float = 1.0) -> "PPRCalculator":
        """Alternative constructor: rebuild an instance from a dict of pre-existing attributes.

        Used to reconstruct a calculator from a previously serialized state (e.g. a toy model)
        rather than from a ModelData / file. The dict is copied onto a blank instance and the
        full pipeline (apply_ecopath_defaults -> optional apply_lim -> _fill_properties ->
        _sort) is re-run so the result is fully consistent.

        Args:
            data_dict (dict): mapping of attribute name -> value to seed the instance with;
                must at least contain a valid _groups_df.
            underdetermined (bool): if True, fill missing flows via apply_lim. Defaults to False.
            zero_catch (bool): treat missing catch as 0. Defaults to True.
            zero_biomass_accum (bool): treat missing biomass accumulation as 0. Defaults to True.
            default_gs (bool): assign default GS=0.2 to regular groups lacking gs. Defaults to True.
            weight_flow (float): LIM minimum-flow penalty weight. Defaults to 1.0.
            weight_guess (float): LIM guess-deviation penalty weight. Defaults to 1.0.

        Returns:
            PPRCalculator: a fully built, sorted instance.
        """
        instance = cls.__new__(cls)
        instance.__dict__.update(data_dict)

        groups_df = deepcopy(instance._groups_df)
        instance.is_underdetermined = underdetermined

        # define all properties:
        groups_df = instance.apply_ecopath_defaults(
            df=groups_df,
            DC=instance.get_DC(DET_as_PP=True),
            det_fate=getattr(instance, "_det_fate", None),
            zero_catch=zero_catch, zero_biomass_accum=zero_biomass_accum, default_gs=default_gs
        )

        if underdetermined:
            groups_df = instance.apply_lim(groups_df, weight_flow=weight_flow, weight_guess=weight_guess)
        
        # fill properties now that groups_df is full:
        instance._fill_properties(groups_df)

        # sort according to seq:
        instance = instance._sort()
        
        return instance

    @classmethod
    def from_modeldata(cls, modeldata: ModelData, underdetermined: bool = False, zero_catch: bool = True, zero_biomass_accum: bool = True, default_gs: bool = True, weight_flow: float = 1.0, weight_guess: float = 1.0) -> "PPRCalculator":
        """Core constructor used by __init__: build the calculator from a loaded ModelData.

        Copies the groups table, diet-composition (DC) matrix, detritus-fate matrix and the
        name<->seq dicts out of the ModelData, verifies at least one detritus group exists,
        then runs the full pipeline: apply_ecopath_defaults -> (optional) apply_lim ->
        _fill_properties -> _sort. Groups are kept sorted by descending seq throughout.

        Args:
            modeldata (ModelData): the loaded model whose groups_data, DC, det_fate and
                seq/name dicts are copied into the new instance.
            underdetermined (bool): if True, fill missing flows via apply_lim. Defaults to False.
            zero_catch (bool): treat missing catch as 0. Defaults to True.
            zero_biomass_accum (bool): treat missing biomass accumulation as 0. Defaults to True.
            default_gs (bool): assign default GS=0.2 to regular groups lacking gs. Defaults to True.
            weight_flow (float): LIM minimum-flow penalty weight. Defaults to 1.0.
            weight_guess (float): LIM guess-deviation penalty weight. Defaults to 1.0.

        Returns:
            PPRCalculator: a fully built, sorted instance.

        Raises:
            Exception: if the model contains no DET (detritus) group.
        """
        instance = cls.__new__(cls)
        instance._model = modeldata
        instance.is_underdetermined = underdetermined

        groups_df = instance._model.groups_data.sort_index(ascending=False)

        instance._groups_df = groups_df
        instance.n_groups = groups_df.shape[0]

        # DC and Detritus fate:
        instance._DC = instance._model.DC.sort_index(ascending=False).sort_index(axis=1, ascending=False).copy()
        instance._det_fate = instance._model.det_fate.sort_index(ascending=False).sort_index(axis=1, ascending=False).copy()

        # useful dicts:
        instance.seq2name = instance._model.seq2name.copy()
        instance.name2seq = instance._model.name2seq.copy()

        # verify at least one DET group exists:
        DET_seq = instance.get_DET_seq()
        if len(DET_seq) == 0:
            raise Exception('no DET group found')

        # Validate the input matrices before building the model: a diet composition must sum to 1
        # per consumer (raises otherwise), and detritus routing must be consistent -- a fully
        # degenerate multi-DET matrix raises, partial rows (legitimate export) only warn.
        ModelData.validate_DC(instance._DC, instance._groups_df)
        ModelData.validate_det_fate(instance._det_fate, instance._groups_df)

        # define all properties:
        groups_df = instance.apply_ecopath_defaults(
            df=groups_df,
            DC=instance.get_DC(DET_as_PP=True),
            det_fate=getattr(instance, "_det_fate", None),
            zero_catch=zero_catch, zero_biomass_accum=zero_biomass_accum, default_gs=default_gs
        )

        if underdetermined:
            groups_df = instance.apply_lim(groups_df, weight_flow=weight_flow, weight_guess=weight_guess)
        
        # fill properties now that groups_df is full:
        instance._fill_properties(groups_df)

        # sort:
        instance._sort()

        return instance

    def _sort(self) -> "PPRCalculator":
        """Normalize the ordering of every Series/DataFrame attribute on the instance.

        Iterates over all attributes and re-sorts each pandas Series by descending index and
        each DataFrame by descending index and columns, so all vectors and matrices share a
        consistent group ordering. _groups_df is deliberately skipped so it keeps its own
        ordering.

        Returns:
            PPRCalculator: self, after in-place re-sorting of the relevant attributes.
        """
        for name, value in vars(self).items():
            if isinstance(value, pd.Series):
                setattr(self, name, value.sort_index(ascending=False))
            elif isinstance(value, pd.DataFrame) and name != '_groups_df':
                setattr(self, name, value.sort_index(ascending=False).sort_index(ascending=False, axis=1))
        return self

    def _fill_properties(self, groups_df: pd.DataFrame) -> None:
        """Unpack the fully-filled groups table into the individual named vectors.

        Splits the completed groups DataFrame into the per-group Series the rest of the class
        relies on (production p, consumption q, catch, predation, growth/biomass_accum,
        immigration/emigration/net_migration, natural mortality M0, respiration, egestion, and
        the ratios EE/GE/GS/TL). Missing columns (detritus_import, tl) are synthesized first.
        Finally runs an initial balance check and stores a mass-balanced copy of the model.

        Args:
            groups_df (pd.DataFrame): the completed per-group parameter table (post defaults/LIM).

        Returns:
            None. Sets the vector attributes, self.is_balanced and self.balanced_model in place.
        """
        self._groups_df = groups_df.copy()

        # add missing columns:
        if 'detritus_import' not in groups_df.columns:
            groups_df['detritus_import'] = 0
        if 'tl' not in groups_df.columns:
            groups_df['tl'] = self.get_TL(break_cycles=False, DET_as_PP=True)

        # important vectors:
        self.p = groups_df['p'].fillna(0).copy()
        self.q = groups_df['q'].fillna(0).copy()
        self.catch = groups_df['catch'].fillna(0).copy()
        self.predation = groups_df['predation'].fillna(0).copy()
        self.growth = groups_df['biomass_accum'].fillna(0).copy()
        self.immigration = groups_df['immigration'].fillna(0).copy()
        self.emigration = groups_df['emigration'].fillna(0).copy()
        self.net_migration = (groups_df['net_migration'] - groups_df['detritus_import']).fillna(0).copy()
        self.M0 = groups_df['M0'].fillna(0).copy()
        self.respiration = groups_df['respiration'].fillna(0).copy()
        self.egestion = groups_df['egestion'].fillna(0).copy()
        # Per-group flow_to_det (= M0 + egestion) exported out of the system (det_fate row < 1).
        self.det_export = (groups_df['det_export'] if 'det_export' in groups_df.columns
                           else pd.Series(0.0, index=groups_df.index)).fillna(0).copy()
        self.EE = groups_df['ee'].fillna(1).copy()
        self.GE = groups_df['ge'].fillna(1).copy()
        self.GS = groups_df['gs'].fillna(0).copy()
        self.TL = groups_df['tl'].fillna(1).copy()

        # balance:
        self.n_balance_runs = 0
        self.is_balanced, _, _ = self.is_model_balanced()
        self.balanced_model = self.balance_model(change_production=False)
    
    @classmethod
    def apply_ecopath_defaults(cls, df: pd.DataFrame, DC: pd.DataFrame, det_fate: Optional[pd.DataFrame] = None, zero_catch: bool = False, zero_biomass_accum: bool = False, default_gs: bool = False) -> pd.DataFrame:
        """Apply standard Ecopath defaults and sync the mass-balance flows with the ratios.

        Deterministically completes the per-group parameter table using the Ecopath
        assumptions: given p, q, ee and catch (and p, q absent for detritus), it derives
        M0 = p*(1-ee), egestion = q*gs (with gs defaulting to 0.2 for regular groups and 0
        for PP/DET/Import), treats missing biomass_accum and net_migration as 0, fills a single
        missing cell of the consumption equation (q = p + respiration + egestion) and of the
        production equation (p = M0 + catch + predation + net_migration + biomass_accum), and
        sets predation as the column sum of the flow matrix Z. Detritus q (= inflow) is built
        from each group's flow_to_det = M0 + egestion routed through det_fate; the unrouted
        remainder (det_fate row < 1) is recorded per group in 'det_export' as flow leaving the
        system.

        Args:
            df (pd.DataFrame): per-group parameter table to complete (modified copy returned).
            DC (pd.DataFrame): diet-composition matrix used to build the flow matrix Z and
                hence predation.
            det_fate (Optional[pd.DataFrame]): matrix giving, per group, the fraction of its
                flow_to_det reaching each detritus column. If None, all detritus inflow is
                pooled into a single total. Defaults to None.
            zero_catch (bool): if True, missing catch values are filled with 0. Defaults to False.
            zero_biomass_accum (bool): if True, missing biomass_accum values are filled with 0.
                Defaults to False.
            default_gs (bool): if True, regular groups lacking gs get the default 0.2.
                Defaults to False.

        Returns:
            pd.DataFrame: the completed table with all flow vectors and the ratios
            ee/ge/gs (and flow_to_det) filled in, sorted by descending seq.
        """
        df = df.sort_index(ascending=False)
        is_regular = df["trophic_info"] == "Regular"
        is_import = df["trophic_info"] == "Import"
        is_det = df["trophic_info"] == "DET"
        is_pp = df["trophic_info"] == "PP"

        # original_nas = df.isna()

        # cols with constant value of 0 or 1:
        cols_where_default_is_zero = ['detritus_import', 'immigration', 'emigration', 'net_migration']
        if zero_catch:
            cols_where_default_is_zero.append('catch')
        if zero_biomass_accum:
            cols_where_default_is_zero.append('biomass_accum')

        fill_dict = {col: 0 for col in cols_where_default_is_zero}
        df = df.fillna(fill_dict)
        df.loc[:, 'biomass'] = df['biomass'].fillna(1)
        # set biomass_accum for det to nan if it was nan before:
        df.loc[is_det, 'biomass_accum'] = np.nan
        # set biomass_accum to 0 for Import groups and and catch to 0 for all non-regular groups:
        df.loc[is_import, 'biomass_accum'] = df.loc[is_import, 'biomass_accum'].fillna(0)
        df.loc[~is_regular, 'catch'] = df.loc[~is_regular, 'catch'].fillna(0)

        # get Z which is correct for all regular groups:
        Z = DC.mul(df['q'].fillna(0), axis='index')
        Z.loc[is_det | is_pp | is_import, :] = 0
        Z = Z.sort_index(ascending=False).sort_index(ascending=False, axis=1)

        # set default values for non-regular groups:
        df.loc[is_pp | is_import | is_det, ['egestion', 'respiration']] = 0
        df.loc[is_import | is_det, 'M0'] = 0
        df.loc[is_import | is_det, 'ee'] = 1

        # fillna for gs and egestion:
        df.loc[~is_regular, 'gs'] = df.loc[~is_regular, 'gs'].fillna(0)
        if default_gs:
            df.loc[is_regular, 'gs'] = df.loc[is_regular, 'gs'].fillna(0.2)

        # Sync gs -> egestion
        mask_sync = df['egestion'].isna() & df['q'].notna() & df['gs'].notna()
        df.loc[mask_sync, 'egestion'] = df['q'] * df['gs']

        # Sync ee -> M0
        mask_sync = df['M0'].isna() & df["p"].notna() & df["ee"].notna()
        df.loc[mask_sync, 'M0'] = df["p"] * (1 - df["ee"])

        # q for detritus:
        df['flow_to_det'] = df['M0'] + df['egestion']
        if det_fate is not None:
            detritus_inflows = det_fate.mul(df['flow_to_det'], axis='index').sum(axis=0)
            df.loc[is_det, 'q'] = df.loc[is_det, 'q'].fillna(detritus_inflows)
        else:
            df.loc[is_det, 'q'] = df.loc[is_det, 'q'].fillna(df['M0'].sum() + df['egestion'].sum())  # set q of detritus as M0.sum() + egestion.sum()

        # predation:
        df.loc[:, 'predation'] = Z.sum(axis=0)

        # for non-regular groups, egestion and respiration are nan so q = p:
        df.loc[is_det, 'p'] = df['q']
        df.loc[is_pp | is_import, 'q'] = df['p']
        
        # pb and qb:
        df.loc[:, 'pb'] = (df['p'] / df['biomass']).fillna(0)
        df.loc[:, 'qb'] = (df['q'] / df['biomass']).fillna(0)

        def _solve_linear_equation(df, eq_cols, signs):
            # Fill a single missing term of a signed linear balance: sum(signs * cols) == 0.
            # Identify rows where exactly ONE variable is missing (otherwise it's unsolvable this way)
            solvable_mask = df[eq_cols].isna().sum(axis=1) == 1

            # Calculate the net sum of all known terms (Pandas ignores NaNs by default in .sum)
            balance_sum = (df[eq_cols] * signs).sum(axis=1, skipna=True)

            # Fill the missing cells
            for col in eq_cols:
                # Target rows where THIS column is the missing one, and the equation is solvable
                target_cells = df[col].isna() & solvable_mask
                
                # The missing value is the negative balance_sum divided by the column's sign
                df.loc[target_cells, col] = -balance_sum[target_cells] / signs[col]
            
            return df
        
        # consumption equation:
        eq_cols = ['q', 'p', 'respiration', 'egestion']
        signs = pd.Series({
            'q': 1, 
            'p': -1, 
            'respiration': -1, 
            'egestion': -1, 
        })
        df = _solve_linear_equation(df, eq_cols, signs)
        
        # production equation:
        eq_cols = ['p', 'M0', 'catch', 'predation', 'net_migration', 'biomass_accum']
        signs = pd.Series({
            'p': 1, 
            'M0': -1, 
            'catch': -1, 
            'predation': -1, 
            'net_migration': -1, 
            'biomass_accum': -1
        })
        df = _solve_linear_equation(df, eq_cols, signs)

        # finalize ratios:
        df.loc[is_det, 'p'] = df['q']
        df.loc[is_pp | is_import, 'q'] = df['p']
        df.loc[is_regular, 'ee'] = 1 - (df['M0'] / df['p'])
        df.loc[is_regular, 'gs'] = (df['egestion'] / df['q'])
        df['ge'] = (df['p'] / df['q'])
        df['flow_to_det'] = df['flow_to_det'].fillna(df['M0'] + df['egestion'])
        det_idx = df.index[is_det]
        # Open-system detritus inflow: each detritus pool receives only the det_fate-routed share
        # of every group's flow_to_det (single- and multi-DET handled identically). With det_fate
        # rows summing to 1 this equals the full flow (closed system); rows summing to <1 mean the
        # shortfall left the system as export (see det_export below).
        if det_fate is not None:
            for det_j in det_idx:
                if det_j in det_fate.columns:
                    df.loc[det_j, 'q'] = (df['flow_to_det'] * det_fate[det_j].reindex(df.index).fillna(0)).sum()
                else:
                    df.loc[det_j, 'q'] = df['flow_to_det'].sum()
        else:
            df.loc[is_det, 'q'] = df['flow_to_det'].sum()
        df.loc[is_det, 'p'] = df.loc[is_det, 'q']
        df.loc[is_det, 'biomass_accum'] = df.loc[is_det, 'p'] - (df.loc[is_det, 'predation'] + df.loc[is_det, 'net_migration'])

        # det_export: the portion of each group's flow_to_det (= M0 + egestion) that is NOT routed
        # to any detritus pool, i.e. exported out of the modeled system. Zero in the closed system
        # (det_fate rows sum to 1) and when det_fate is unavailable; positive for groups whose
        # det_fate row sums to <1 (e.g. seabirds/mammals/migratory groups dying outside the box).
        if det_fate is not None:
            det_fate_rowsum = det_fate.sum(axis=1).reindex(df.index).fillna(1.0)
            df['det_export'] = ((1.0 - det_fate_rowsum).clip(lower=0.0) * df['flow_to_det']).fillna(0.0)
        else:
            df['det_export'] = 0.0

        return df

    def apply_lim(self, df: pd.DataFrame, weight_flow: float = 1.0, weight_guess: float = 0.0) -> pd.DataFrame:
        """Fill missing mass-balance variables via a per-group Linear Inverse Model (SLSQP).

        For each group with one or more unknowns among {q, p, respiration, egestion, M0,
        biomass_accum}, builds smart biological initial guesses, dynamically registers only
        the constraints that involve a free variable (consumption balance, production balance,
        EE bounds [0, 0.95], GS bounds [0.10, 0.35]), and minimizes a weighted sum of a
        parsimony (minimum total flow) penalty and a deviation-from-guess penalty. On failure
        the original guesses are kept and a message is printed. Detritus rows are rebalanced
        from flow_to_det / det_fate afterwards.

        Args:
            df (pd.DataFrame): per-group parameter table with possible NaNs to be solved.
            weight_flow (float): weight on the minimum-total-flow (parsimony) penalty.
                Defaults to 1.0.
            weight_guess (float): weight on the squared deviation from the biological guesses.
                Defaults to 0.0.

        Returns:
            pd.DataFrame: the table with solved flows and recomputed ee/gs/ge/flow_to_det and
            rebalanced detritus rows.
        """
        df_out = df.copy()

        def _finalize_outputs(df_final):
            # EE
            df_final.loc[:, 'ee'] = np.where(
                df_final.loc[:, 'p'] != 0, 
                1 - (df_final.loc[:, 'M0'] / df_final.loc[:, 'p']), 
                0
            )
            # GS
            df_final.loc[:, 'gs'] = np.where(
                df_final.loc[:, 'q'] != 0, 
                df_final.loc[:, 'egestion'] / df_final.loc[:, 'q'], 
                0
            )
            # GE
            df_final['ge'] = np.where(df_final['q'] != 0, df_final['p'] / df_final['q'], 0)
            
            # Flow to Detritus
            df_final['flow_to_det'] = df_final['M0'] + df_final['egestion']

            # rebalance DET rows: each pool receives only the det_fate-routed share of every
            # group's flow_to_det (single- and multi-DET identical), consistent with the open
            # system used in apply_ecopath_defaults. det_export is the unrouted remainder.
            det_seqs = self.get_DET_seq()
            det_fate = getattr(self, '_det_fate', None)
            if det_fate is not None:
                for det_j in det_seqs:
                    if det_j in det_fate.columns:
                        df_final.loc[det_j, 'q'] = (df_final['flow_to_det'] * det_fate[det_j].reindex(df_final.index).fillna(0)).sum()
                    else:
                        df_final.loc[det_j, 'q'] = df_final['flow_to_det'].sum()
                det_fate_rowsum = det_fate.sum(axis=1).reindex(df_final.index).fillna(1.0)
                df_final['det_export'] = ((1.0 - det_fate_rowsum).clip(lower=0.0) * df_final['flow_to_det']).fillna(0.0)
            else:
                df_final.loc[det_seqs, 'q'] = df_final['flow_to_det'].sum()
                df_final['det_export'] = 0.0
            df_final.loc[det_seqs, 'p'] = df_final.loc[det_seqs, 'q']
            df_final.loc[det_seqs, 'biomass_accum'] = df_final.loc[det_seqs, 'p'] - (df.loc[det_seqs, 'predation'] + df.loc[det_seqs, 'net_migration'])

            return df_final

        cols_to_solve = ['q', 'p', 'respiration', 'egestion', 'M0', 'biomass_accum']
        records = df_out.to_dict('index')
        
        for idx, row in records.items():
            missing_cols = [c for c in cols_to_solve if pd.isna(row[c])]
            if not missing_cols:
                continue
                
            c_val = row['catch'] if pd.notna(row['catch']) else 0
            pr_val = row['predation'] if pd.notna(row['predation']) else 0
            nm_val = row['net_migration'] if pd.notna(row['net_migration']) else 0

            # --- 1. Coordinated Guess Variables ---
            q_guess = row['q'] if pd.notna(row['q']) else (row['p'] * 10 if pd.notna(row['p']) else 1.0)
            p_guess = row['p'] if pd.notna(row['p']) else (q_guess * 0.1)
            
            # Smart adaptation: If both P and Q are known, don't overshoot the available pool
            if pd.notna(row['q']) and pd.notna(row['p']):
                rem_pool = max(0.0, row['q'] - row['p'])
                resp_guess = row['respiration'] if pd.notna(row['respiration']) else (rem_pool * (0.7 / 0.9))
                egest_guess = row['egestion'] if pd.notna(row['egestion']) else (rem_pool * (0.2 / 0.9))
            else:
                resp_guess = row['respiration'] if pd.notna(row['respiration']) else (q_guess * 0.7)
                egest_guess = row['egestion'] if pd.notna(row['egestion']) else (q_guess * 0.2)
                
            m0_guess = row['M0'] if pd.notna(row['M0']) else (q_guess * 0.1)
            ba_guess = row['biomass_accum'] if pd.notna(row['biomass_accum']) else (p_guess - m0_guess - c_val - pr_val - nm_val)

            # Assemble variables and bounds
            x0 = []
            bounds = []
            for col in missing_cols:
                if col == 'q': x0.append(q_guess); bounds.append((0, None))
                elif col == 'p': x0.append(p_guess); bounds.append((0, None))
                elif col == 'respiration': x0.append(resp_guess); bounds.append((0, None))
                elif col == 'egestion': x0.append(egest_guess); bounds.append((0, None))
                elif col == 'M0': x0.append(m0_guess); bounds.append((0, None))
                elif col == 'biomass_accum': x0.append(ba_guess); bounds.append((None, None))

            x0_arr = np.array(x0)
            scale = np.where(np.abs(x0_arr) > 0.01, np.abs(x0_arr), 1.0)
            
            # --- 2. High-Speed Mapping Indices ---
            idx_q = missing_cols.index('q') if 'q' in missing_cols else -1
            idx_p = missing_cols.index('p') if 'p' in missing_cols else -1
            idx_r = missing_cols.index('respiration') if 'respiration' in missing_cols else -1
            idx_e = missing_cols.index('egestion') if 'egestion' in missing_cols else -1
            idx_m0 = missing_cols.index('M0') if 'M0' in missing_cols else -1
            idx_ba = missing_cols.index('biomass_accum') if 'biomass_accum' in missing_cols else -1

            # --- 3. Objective Function ---
            def objective(x):
                # Penalty A: Flow Penalty (Parsimony / Minimum Total Flow)
                # Minimizes the total magnitude of the actively solved flows
                flow_pen = np.sum(np.abs(x) / scale)
                
                # Penalty B: Guess Penalty (Biological Anchoring)
                # Minimizes squared deviation from your smart biological guesses (x0)
                guess_pen = np.sum(((x - x0_arr) / scale) ** 2)
                
                return (weight_flow * flow_pen) + (weight_guess * guess_pen)

            # --- 4. Dynamic Constraint Registration ---
            constraints = []
            
            # Consumption Balance (Only if Q, P, R, or U are being optimized)
            if any(c in missing_cols for c in ['q', 'p', 'respiration', 'egestion']):
                def eq_consumption(x):
                    q_curr = x[idx_q] if idx_q >= 0 else q_guess
                    p_curr = x[idx_p] if idx_p >= 0 else p_guess
                    r_curr = x[idx_r] if idx_r >= 0 else resp_guess
                    e_curr = x[idx_e] if idx_e >= 0 else egest_guess
                    return q_curr - (p_curr + r_curr + e_curr)
                constraints.append({'type': 'eq', 'fun': eq_consumption})
                
            # Production Balance (Only if P, M0, or BA are being optimized)
            if any(c in missing_cols for c in ['p', 'M0', 'biomass_accum']):
                def eq_production(x):
                    p_curr = x[idx_p] if idx_p >= 0 else p_guess
                    m0_curr = x[idx_m0] if idx_m0 >= 0 else m0_guess
                    ba_curr = x[idx_ba] if idx_ba >= 0 else ba_guess
                    return p_curr - (m0_curr + c_val + pr_val + nm_val + ba_curr)
                constraints.append({'type': 'eq', 'fun': eq_production})

            # EE bounds (Only if P or M0 can actually change)
            if 'p' in missing_cols or 'M0' in missing_cols:
                # Upper Bound: EE <= 0.95 (Equivalent to: M0 - 0.05 * P >= 0)
                def ineq_ee_upper(x):
                    p_curr = x[idx_p] if idx_p >= 0 else p_guess
                    m0_curr = x[idx_m0] if idx_m0 >= 0 else m0_guess
                    return m0_curr - (0.05 * p_curr)
                constraints.append({'type': 'ineq', 'fun': ineq_ee_upper})
                
                # Lower Bound: EE >= 0.0 (Equivalent to: P - M0 >= 0)
                def ineq_ee_lower(x):
                    p_curr = x[idx_p] if idx_p >= 0 else p_guess
                    m0_curr = x[idx_m0] if idx_m0 >= 0 else m0_guess
                    return p_curr - m0_curr
                constraints.append({'type': 'ineq', 'fun': ineq_ee_lower})
                
            # GS Upper limit (Only if Q or Egestion can actually change)
            if 'q' in missing_cols or 'egestion' in missing_cols:
                # Enforce GS <= 0.35
                def ineq_gs_upper(x):
                    q_curr = x[idx_q] if idx_q >= 0 else q_guess
                    e_curr = x[idx_e] if idx_e >= 0 else egest_guess
                    return 0.35 * q_curr - e_curr
                constraints.append({'type': 'ineq', 'fun': ineq_gs_upper})

                # Enforce GS >= 0.10
                def ineq_gs_lower(x):
                    q_curr = x[idx_q] if idx_q >= 0 else q_guess
                    e_curr = x[idx_e] if idx_e >= 0 else egest_guess
                    return e_curr - (0.10 * q_curr)
                constraints.append({'type': 'ineq', 'fun': ineq_gs_lower})

            # --- 5. Run Optimizer ---
            res = minimize(
                objective, 
                x0_arr, 
                method='SLSQP', 
                bounds=bounds, 
                constraints=constraints,
                tol=1e-5,
                options={'maxiter': 1500, 'ftol': 1e-5}
            )

            if not res.success:
                model = self.get_model()
                if model is not None:
                    print(f"Model {model.model_number} - LIM failed for group '{idx}': {res.message}")
                else:
                    print(f"Toy Model - LIM failed for group '{idx}': {res.message}")            
            else:
                for i, col in enumerate(missing_cols):
                    df_out.loc[idx, col] = res.x[i]                

        return _finalize_outputs(df_out)

    # balancing checks:
    def is_model_balanced(self) -> tuple[bool, pd.Series, pd.Series]:
        """Check the two Ecopath mass-balance identities hold (within tolerance).

        Recomputes production = catch + predation + growth + net_migration + M0 and
        consumption = production + egestion + respiration from the current vectors and compares
        them (np.isclose) against the stored p and q.

        Returns:
            tuple[bool, pd.Series, pd.Series]: (model_is_balanced, production, consumption),
            where model_is_balanced is True iff both identities hold, and the two Series are the
            recomputed production and consumption per group so callers can inspect the flows.
        """
        # production = catch + predation + growth + net_migration + M0
        production = self.catch + self.predation + self.growth + self.net_migration + self.M0
        p_is_balanced = all(np.isclose(production, self.p))

        # consumpotion = production + egestion + respiration
        consumption = production + self.egestion + self.respiration
        q_is_balanced =  all(np.isclose(consumption, self.q, atol=1e-6))

        model_is_balanced = p_is_balanced and q_is_balanced

        return model_is_balanced, production, consumption
    
    def balance_model(self, change_production: bool = False) -> "PPRCalculator":
        """Iteratively force the model onto exact mass balance and return a balanced copy.

        Works on a deepcopy so self is untouched. Each iteration recomputes growth from the
        consumption identity, then either net_migration (change_production=False) or production
        (change_production=True) from the production identity, and re-derives predation (column
        sum of Z) and M0 = p*(1-EE). It repeats until predation and M0 stop changing, i.e. the
        flows are self-consistent.

        Args:
            change_production (bool): if False (default) the residual is absorbed into growth
                and net_migration, keeping production p fixed; if True the residual is absorbed
                into production p instead.

        Returns:
            PPRCalculator: a deepcopied, mass-balanced instance (n_balance_runs incremented).
        """
        # TODO: change this function so I can decide which subset of parameters stays constant
        balanced_self = deepcopy(self)
        close_enough = False
        while not close_enough:
            balanced_self.growth = (
                balanced_self.q - (balanced_self.egestion + balanced_self.respiration + 
                (balanced_self.catch + balanced_self.net_migration + balanced_self.M0 + balanced_self.predation))
            )
            if not change_production:
                balanced_self.net_migration = balanced_self.p - (balanced_self.catch + balanced_self.M0 + balanced_self.predation + balanced_self.growth)
            else:
                balanced_self.p = balanced_self.catch + balanced_self.net_migration + balanced_self.M0 + balanced_self.predation + balanced_self.growth
            predation = balanced_self.get_Z(DET_as_PP=True).sum(axis=0)
            M0 = balanced_self.p * (1 - balanced_self.EE)
            close_enough = (
                all(np.isclose(predation, balanced_self.predation))
                and all(np.isclose(M0, balanced_self.M0))
            )
            balanced_self.predation = predation
            balanced_self.M0 = balanced_self.p * (1 - balanced_self.EE)
            balanced_self.n_balance_runs += 1
        
        return balanced_self
    
    def is_sppr_balanced(self, sppr: pd.DataFrame | pd.Series, diet_import_equations: Optional[tuple] = None) -> tuple[bool, float, float]:
        """Check an SPPR result is globally self-consistent (inflow == outflow).

        The primary-production inflow into the system must equal the production-required
        outflow implied by the exports (catch + growth + net_migration) weighted by each
        group's SPPR.

        Two regimes:
          - diet_import_equations is None (the '_old' helper): PP and Import groups form the
            production basis and inflow is their total production p.
          - diet_import_equations given (the '_new' helper): only PP forms the within-system
            basis, and the symbolic diet-import equations are solved so imported diet
            contributes its own DIET_SPPR-weighted inflow on top.

        Args:
            sppr (pd.DataFrame | pd.Series): an SPPR result from one of the SPPR_* methods;
                DataFrames are summed across columns to a per-group total.
            diet_import_equations (Optional[tuple]): (equations, variables) sympy system from a
                symbolic SPPR call; if provided the diet-import-aware ('_new') check is used.
                Defaults to None.

        Returns:
            tuple[bool, float, float]: (is_balanced, inflow, outflow), where is_balanced is
            np.isclose(inflow, outflow).
        """
        def _helper_old(sppr):
            sppr = PPRCalculator.rename_results(sppr, self.name2seq)

            if isinstance(sppr, pd.DataFrame):
                sppr = sppr.sum(axis=1)  # fix when sppr includes sppr_import: then we should sum only within the region

            basis_seq = list(self.get_PP_seq()) + list(self.get_Import_seq())
            inflow = self.p[basis_seq].sum()
            outflow = (self.catch + self.growth + self.net_migration).mul(sppr, fill_value=1).sum()

            is_balanced = bool(np.isclose(inflow, outflow))
            return is_balanced, inflow, outflow

        def _helper_new(sppr, diet_import_equations):
            
            equations, variabls = diet_import_equations

            sppr = PPRCalculator.rename_results(sppr, self.name2seq)
            basis_seq = list(self.get_PP_seq())
            inflow_from_within = self.p[basis_seq].sum()

            sol = sm.linsolve(equations, variabls)
            sol_tuple = list(sol)[0]
            sol_dict = dict(zip(variabls, sol_tuple))
            PP_dict = {sm.Symbol(f'SPPR_{self.seq2name[s]}'.replace(' ', '_')): 1 for s in self.get_PP_seq()}
            d = {k: float(v.subs(PP_dict)) for k, v in sol_dict.items()}
            r_dict = {k: sm.Symbol(f'DIET_SPPR_{self.seq2name[k]}'.replace(' ', '_')) for k in self.GE.index}
            a = self.get_Z()[self.get_Import_seq()].squeeze().rename(index=r_dict)
            inflow_from_diet_import = a.mul(pd.Series(d), fill_value=0).sum()

            inflow = inflow_from_within + inflow_from_diet_import
            outflow = (self.catch + self.growth + self.net_migration).mul(sppr.sum(axis=1), fill_value=1).sum()

            is_balanced = bool(np.isclose(inflow, outflow))
            return is_balanced, inflow, outflow
        
        if diet_import_equations is None:
            return _helper_old(sppr)
        else:
            return _helper_new(sppr, diet_import_equations)
    # getters:
    def get_model(self) -> Optional[ModelData]:
        """Return the underlying ModelData, or None for toy / from_dict instances.

        Returns:
            Optional[ModelData]: the backing ModelData if this instance was built from a
            file/model number, else None.
        """
        if hasattr(self, "_model"):
            return self._model
        else:
            return None

    def get_groups_df(self) -> pd.DataFrame:
        """Return a defensive (descending-seq sorted) copy of the per-group parameter table.

        Returns:
            pd.DataFrame: a copy of _groups_df sorted by descending seq.
        """
        return self._groups_df.copy().sort_index(ascending=False)

    def get_DC(self, DET_as_PP: bool = True, normalize: bool = False) -> pd.DataFrame:
        """Return the diet-composition (DC) matrix, optionally redefining detritus rows.

        Args:
            DET_as_PP (bool, optional): if True, the DET row is the stored DC (detritus treated
                as a basal source, row set to 1). If False, the DET row is rebuilt from the flow
                matrix as (M0 + egestion) shares, i.e. Z/Z.sum, with the non-DET rows rescaled to
                preserve their original row sums. Defaults to True.
            normalize (bool, optional): if True, every DC row is renormalized to sum to 1;
                otherwise rows keep their original sum. Defaults to False.

        Returns:
            pd.DataFrame: the (n x n) DC matrix, sorted by descending index and columns.
        """
        if DET_as_PP:
            DC = self._DC.copy()
        else:
            Z = self.get_Z(DET_as_PP=False)
            DC = Z.div(Z.sum(axis=1), axis=0).fillna(0)
            # keep the sum correct:
            non_DET_seq = [i for i in self._DC.index if i not in self.get_DET_seq()]
            DC.loc[non_DET_seq] = DC.loc[non_DET_seq].multiply(self._DC.sum(axis=1).loc[non_DET_seq], axis=0)

        if normalize:
            DC = DC.div(DC.sum(axis=1), axis=0).fillna(0)
        return DC.sort_index(ascending=False).sort_index(ascending=False, axis=1)

    def get_det_fate(self):
        return self._det_fate.copy()
    
    def get_Z(self, DET_as_PP: bool = False) -> pd.DataFrame:
        """Return the flow matrix Z = DC * q (consumption-weighted diet), with DET rows redefined.

        Z[i, j] is the flow from prey j into predator i. For non-detritus rows it is the diet
        composition scaled by the group's consumption q. The DET rows are special:

        Args:
            DET_as_PP (bool, optional): if False (default), each DET row is the group's
                flow_to_det = (M0 + egestion) split across prey by det_fate (or assigned whole
                when no fate column exists). If True, DET rows are set to 0 (detritus treated as
                a basal source contributing no outgoing diet flow). Defaults to False.

        Returns:
            pd.DataFrame: the (n x n) flow matrix Z, sorted by descending index and columns.
        """
        Z = self._DC.mul(self._groups_df['q'].fillna(0), axis='index')
        DET_seq = self.get_DET_seq()
        if not DET_as_PP:
            flow_to_det = (self.M0 + self.egestion).fillna(0)
            det_fate = getattr(self, '_det_fate', None)
            # Route each DET row's flow_to_det through det_fate for single- and multi-DET alike
            # (open system); with det_fate rows summing to 1 this reduces to the full flow.
            if det_fate is not None:
                for det_j in DET_seq:
                    if det_j in det_fate.columns:
                        fracs = det_fate[det_j].reindex(flow_to_det.index).fillna(0)
                        Z.loc[det_j, :] = flow_to_det * fracs
                    else:
                        Z.loc[det_j, :] = flow_to_det
            else:
                for det_j in DET_seq:
                    Z.loc[det_j, :] = flow_to_det
        else:
            for det_j in DET_seq:
                Z.loc[det_j, :] = 0
        return Z.sort_index(ascending=False).sort_index(ascending=False, axis=1)
    
    def get_DET_seq(self) -> list:
        """Return the sorted seq IDs of all detritus (DET) groups.

        Returns:
            list: ascending-sorted sequence numbers of the DET groups.
        """
        return sorted(self._groups_df.index[self._groups_df['trophic_info'] == 'DET'].values)

    def get_PP_seq(self) -> list:
        """Return the sorted seq IDs of all primary-producer (PP) groups.

        Returns:
            list: ascending-sorted sequence numbers of the PP groups.
        """
        return sorted(self._groups_df.index[self._groups_df['trophic_info'] == 'PP'].values)

    def get_Regular_seq(self) -> list:
        """Return the sorted seq IDs of all regular (consumer) groups.

        Returns:
            list: ascending-sorted sequence numbers of the Regular groups.
        """
        return sorted(self._groups_df.index[self._groups_df['trophic_info'] == 'Regular'].values)

    def get_Import_seq(self) -> list:
        """Return the sorted seq IDs of all imported-diet (Import) groups.

        Import groups represent an external production source feeding the modeled system.

        Returns:
            list: ascending-sorted sequence numbers of the Import groups.
        """
        return sorted(self._groups_df.index[self._groups_df['trophic_info'] == 'Import'].values)

    def get_TE(self, TE_option: str, DET_values: float = 1, as_matrix: bool = True, global_TE: str | float = 'mean') -> pd.DataFrame | pd.Series:
        """Build the per-group transfer-efficiency (TE) vector or matrix.

        Args:
            TE_option (str): selects how each group's transfer efficiency is computed; one of:
                'GE' -> gross efficiency p/q;
                'TE' -> gross efficiency times the ecotrophic fraction, (p/q)*(1 - M0/p);
                'With Egestion' -> (p/q)*(q/(q - egestion)), i.e. efficiency on assimilated intake;
                'global' -> a single scalar TE broadcast to all groups (see global_TE).
            DET_values (float, optional): TE assigned to the detritus rows (transfer from
                detritus to others). Defaults to 1.
            as_matrix (bool, optional): if True return an (n x n) DataFrame (the vector
                broadcast across columns); if False return the length-n Series. Defaults to True.
            global_TE (str | float, optional): used only when TE_option == 'global'. If 'mean'
                (default) the global value is the catch-weighted mean of the per-group 'TE'
                efficiency (biomass-weighted when total catch is 0); otherwise it is used as the
                literal global TE value.

        Returns:
            pd.DataFrame | pd.Series: the TE matrix if as_matrix else the TE Series, sorted by
            descending index (and columns for the matrix).

        Raises:
            Exception: if TE_option is not one of the four supported strings.
        """
        TE_options = ['GE', 'TE', 'With Egestion', 'global']
        if TE_option == 'GE':
            te = (self.p / self.q).fillna(1)
        elif TE_option == 'TE':
            te = ((self.p / self.q).fillna(1)) * (1 - self.M0 / self.p).fillna(1)
        elif TE_option == 'With Egestion':
            # te = ((self.p + self.egestion) / self.q).fillna(1)
            te = (self.p / self.q).fillna(1) * (self.q / (self.q - self.egestion)).fillna(1)
        elif TE_option == 'global':
            if global_TE == 'mean':
                weights=self.catch
                if weights.sum() == 0:
                    weights = self.get_groups_df()['biomass'].fillna(0)
                te = ((self.p / self.q).fillna(1)) * (1 - self.M0 / self.p).fillna(1)
                te = np.average(te, weights=weights)
            else:
                te = global_TE
            te = pd.Series(np.ones(self.n_groups) * te, index=self.GE.index)
        else:
            raise Exception(f"TE_option should be one of {TE_options}")

        te.loc[self.get_DET_seq()] = DET_values
        # te.loc[self.get_DET_seq()] = self.EE[self.get_DET_seq()]

        if as_matrix:
            te = pd.DataFrame([te] * len(te), index=te.index, columns=te.index).T
            te = te.sort_index(axis=1, ascending=False)
        
        return te.sort_index(ascending=False)
    
    def get_TL(self, break_cycles: bool, DET_as_PP: bool, TE_option: str = 'With Egestion') -> pd.Series:
        """Compute the trophic level of every group via the standard linear-algebra definition.

        Solves TL = (I - DC)^-1 . 1, i.e. each group's TL is one plus the diet-weighted mean TL
        of its prey. Detritus rows of DC are first redefined per TE_option so basal sources sit
        at TL 1.

        Args:
            break_cycles (bool): if True, remove cycles from the flow matrix (Ulanowicz
                algorithm) and rebuild DC before inverting, so the inversion is well-behaved.
            DET_as_PP (bool): passed through to get_Z / get_DC; whether detritus is treated as
                a basal production source.
            TE_option (str, optional): how detritus rows of DC are redefined before inversion:
                'TE' zeroes detritus feeding regular groups; 'GE' rebuilds the detritus diet
                from M0 / det_fate; other options leave DC as-is. Defaults to 'With Egestion'.

        Returns:
            pd.Series: per-group trophic level, indexed by seq, sorted descending.
        """
        DET_seq = self.get_DET_seq()
        Regular_seq = self.get_Regular_seq()
        
        Z = self.get_Z(DET_as_PP=DET_as_PP)
        DC = self.get_DC(DET_as_PP=DET_as_PP)

        # breake cycles:
        if break_cycles:
            Z = remove_cycles(Z, new=False)
            DC = Z.div(Z.sum(axis=1), axis=0).fillna(0)
        
        # change DC according to TE_option:
        flow2det = (self.M0 + self.egestion).sum()
        if TE_option == 'TE':
            DC.loc[DET_seq, Regular_seq] = 0
        elif TE_option == 'GE':
            det_fate = getattr(self, '_det_fate', None)
            for det_j in DET_seq:
                q_det_j = Z.loc[det_j, :].sum() if not DET_as_PP else flow2det / max(len(DET_seq), 1)
                if det_fate is not None and det_j in det_fate.columns and q_det_j > 0:
                    fracs = det_fate[det_j].reindex(self.M0.index).fillna(0)
                    DC.loc[det_j, :] = ((self.M0 * fracs) / q_det_j).fillna(0)
                else:
                    DC.loc[det_j, :] = ((self.M0) / flow2det).fillna(0)

        # calculate TL:
        B = np.ones(DC.shape[0])
        TL = np.linalg.inv((np.identity(DC.shape[0]) - DC)) @ B
        TL = pd.Series(TL, index=DC.index).sort_index(ascending=False)
        return TL

    def get_PPR(self, sppr: pd.DataFrame | pd.Series, only_inner: bool = False, only_pp: bool = False) -> pd.DataFrame | pd.Series:
        """Convert a per-group SPPR into total primary production required (PPR) by the catch.

        Weights each group's SPPR (primary production required per unit production) by its
        catch and sums: PPR = catch . SPPR. The SPPR input is relabeled to seq, reindexed onto
        the catch index, and infinities are zeroed first.

        Args:
            sppr (pd.DataFrame | pd.Series): an SPPR result from one of the SPPR_* methods.
            only_inner (bool, optional): if True, drop the Import columns so only
                within-system production is counted. Defaults to False.
            only_pp (bool, optional): if True, drop the Import and Detritus columns so only
                within-system primary production is counted. Defaults to False.

        Returns:
            pd.DataFrame | pd.Series: total PPR per basal source; a 1-row DataFrame if sppr is a
            DataFrame, else a Series.
        """
        sppr = sppr.copy().fillna(0) # sppr is an output of an SPPR calculating method from this class.
        sppr = PPRCalculator.rename_results(sppr, self.name2seq)
        sppr = sppr.reindex(self.catch.index, fill_value=0)
        sppr = sppr.replace(np.inf, 0)

        # only_pp is a stronger option than only_inner
        if only_pp: only_inner=False

        if only_inner and (set(self.get_Import_seq()).issubset(set(sppr.columns))):
            sppr = sppr.drop(columns=self.get_Import_seq())
        
        if only_pp and ((set(self.get_Import_seq()) | set(self.get_DET_seq())).issubset(set(sppr.columns))):
            sppr = sppr.drop(columns=self.get_Import_seq())
            sppr = sppr.drop(columns=self.get_DET_seq())

        if isinstance(sppr, pd.DataFrame):
            return self.catch.dot(sppr).to_frame().T
        else:
            return self.catch.dot(sppr)
    
    def get_NPP(self, only_inner: bool = True) -> float:
        """Return the net primary production (NPP) of the system.

        NPP is the total production p of the primary-producer (PP) groups.

        Args:
            only_inner (bool, optional): must be True; only the within-system case is
                implemented. Defaults to True.

        Returns:
            float: summed production of the PP groups.

        Raises:
            Exception: if only_inner is False (not implemented).
        """
        if only_inner:
            return self.p[self.get_PP_seq()].sum()
        else:
            raise Exception('not implemented yet')

    def get_PPR2NPP_ratio(self, sppr: pd.DataFrame | pd.Series, only_pp: bool = False) -> float:
        """Return the fraction of available NPP appropriated by the catch (PPR / NPP).

        Args:
            sppr (pd.DataFrame | pd.Series): an SPPR result from one of the SPPR_* methods.
            only_pp (bool, optional): if True, drop the Import and Detritus columns so only
                within-system primary production is counted. Defaults to False.

        Returns:
            float: total within-system PPR divided by total NPP.
        """
        return self.get_PPR(sppr, only_inner=True, only_pp=only_pp).sum(axis=1).sum() / self.get_NPP(only_inner=True)
        
    ##########################################################################################
    ############################## SPPR calculating methods ##################################
    ##########################################################################################
    def SPPR_1986(self) -> pd.DataFrame:
        """Pauly & Christensen (1986)-style SPPR using a single catch-weighted trophic level.

        Computes one catch-weighted mean trophic level for the whole catch, then SPPR =
        TE^(1-TL) with a fixed transfer efficiency TE = 0.1. The same value is assigned to
        every group (constant for a given model). Returns all zeros if there is no catch.

        Returns:
            pd.DataFrame: a single 'sppr' column indexed by group seq.
        """
        if all(self.catch == 0):
            return pd.DataFrame(0, index=self.GE.index, columns=['sppr'])
        TE = 0.1
        TL = np.average(self.get_TL(break_cycles=True, DET_as_PP=True), weights=self.catch)  # break cylces in DC, force DET to be of TL=1 before TL calculation.
        SPPR = TE ** (1 - TL)
        SPPR = pd.DataFrame(SPPR, index=self.GE.index, columns=['sppr'])
        return SPPR
    
    def SPPR_1995(self, global_TE: str | float = 0.1) -> pd.DataFrame:
        """Christensen & Pauly (1995)-style per-group SPPR = TE^(1-TL).

        Uses each group's own continuous trophic level together with a single global transfer
        efficiency.

        Args:
            global_TE (str | float, optional): the global TE; either the literal float or
                'mean' (catch- / biomass-weighted mean efficiency, see get_TE). Defaults to 0.1.

        Returns:
            pd.DataFrame: a single 'sppr' column indexed by group seq.
        """
        TE = self.get_TE(TE_option='global', global_TE=global_TE, as_matrix=False)
        TL = self.get_TL(break_cycles=True, DET_as_PP=True)  # break cylces in DC and force DET to be of TL=1.
        SPPR = TE ** (1 - TL)
        return SPPR.to_frame(name='sppr')
    
    def SPPR_1995_TL_fix(self, global_TE: str | float = 0.1) -> pd.DataFrame:
        """SPPR_1995 variant that linearly interpolates between bracketing integer trophic levels.

        Instead of raising 1/TE to a non-integer power directly, it blends the two integer
        levels bracketing each group's fractional TL:
        SPPR = (1-frac)*(1/TE)^(TLint-1) + frac*(1/TE)^TLint, which avoids the discontinuity of
        the direct fractional exponent.

        Args:
            global_TE (str | float, optional): the global TE; literal float or 'mean'
                (see get_TE). Defaults to 0.1.

        Returns:
            pd.DataFrame: a single 'sppr' column indexed by group seq.
        """
        TE = self.get_TE(TE_option='global', global_TE=global_TE, as_matrix=False)
        TL = self.get_TL(break_cycles=True, DET_as_PP=True)  # break cylces in DC and force DET to be of TL=1.
        TL_fraction = TL % 1
        TL_int = TL.astype(int)
        sppr = (1-TL_fraction) * (1/TE)**(TL_int-1) + TL_fraction * (1/TE)**(TL_int)
        sppr = sppr.fillna(1)
        return sppr.to_frame(name='sppr')
        
    def SPPR_EwE(self, TE_option: str, use_EE: bool = True, return_paths: bool = True, silent: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
        """Path-enumeration SPPR in the style of Ecopath with Ecosim (EwE) flow-network analysis.

        Enumerates every simple path from each group down to a basal terminal node (a node with
        an all-zero DC row) and sums the product of the edge weights along each path.
        A = DC/TE is the per-edge production-required weight, so for each (group, basal-source)
        pair SPPR is the sum over all paths of the product of A along the path. When use_EE, the
        rows are scaled by ecotrophic efficiency.

        Args:
            TE_option (str): transfer-efficiency mode for building A; one of 'GE', 'TE',
                'With Egestion', 'global' (see get_TE).
            use_EE (bool, optional): if True, scale each row by the group's ecotrophic
                efficiency EE. Defaults to True.
            return_paths (bool, optional): if True use the slow implementation that also returns
                the explicit per-(source, sink) path lists; if False use the fast vectorized
                implementation that returns an empty dict for paths. Defaults to True.
            silent (bool, optional): if True, suppress the tqdm progress bars. Defaults to True.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, dict]: (SPPR, A, paths_dict), where SPPR is the
            (groups x basal-terminals) production-required matrix, A is the per-edge weight
            matrix DC/TE, and paths_dict maps source -> sink -> list of node-name paths (empty
            dict when return_paths is False).
        """
        def _get_paths_with_safety_valve(g, start_node, terminal_indices, max_paths=1_000_000):
            """
            Increments depth until either all paths are found
            or the 10 million path limit is breached.
            """
            final_paths = []

            # Define steps of depth.
            # Most food webs are 'thin' enough that steps of 5 are fine.
            # We go up to 45 as you requested earlier.
            n_paths_prev = -1
            for depth in range(1, len(self.EE), 2):
                try:
                    # Get paths at current depth limit
                    current_paths = g.get_all_simple_paths(
                        start_node,
                        to=terminal_indices,
                        maxlen=depth
                    )

                    # Update our collection
                    final_paths = current_paths

                    # Check if we've hit the "explosion" threshold
                    n_paths = len(final_paths)
                    # print(start_node, depth, n_paths)
                    # print(current_paths)
                    # print("=" * 50)

                    if n_paths >= max_paths:
                        # We stop here because the next depth step
                        # will likely hang or exceed memory
                        break

                    # OPTIONAL: If the number of paths didn't increase from the
                    # last depth, it means we've found ALL possible paths already.
                    # (Checking this requires keeping the previous count)
                    if n_paths_prev == n_paths:
                        break

                    n_paths_prev = n_paths

                except Exception as e:
                    # If igraph itself runs out of memory or hits an internal error
                    print(f"Depth {depth} too complex, returning paths from depth {depth - 3}: {e}")
                    break

            return final_paths
        def _slow_EwE_with_paths(TE_option, use_EE, silent):
            DC = self.get_DC(DET_as_PP=True)
            TE = self.get_TE(TE_option=TE_option, as_matrix=True, DET_values=1)
            A = (DC / TE).fillna(0)

            # 1. Cache nodes as a tuple for hyper-fast, memory-efficient string referencing
            nodes_tuple = tuple(DC.index)
            num_nodes = len(nodes_tuple)

            # 2. Convert matrix to a pure Python list-of-lists.
            # Indexing A_list[i][j] is MUCH faster than A_vals[i, j] in a tight loop.
            A_list = A.values.tolist()
            DC_vals = DC.values

            # 3. Identify terminals
            terminal_mask = (DC_vals == 0).all(axis=1)
            terminal_indices = np.where(terminal_mask)[0].tolist()
            terminal_nodes = [nodes_tuple[i] for i in terminal_indices]
            num_terminals = len(terminal_indices)

            # Map terminal indices to their column index in the final 2D array
            term_idx_to_col = {t_idx: col for col, t_idx in enumerate(terminal_indices)}

            # 4. Pre-allocate the SPPR matrix and Output Dictionary
            sppr_matrix = np.zeros((num_nodes, num_terminals), dtype=np.float64)
            paths_dict = {n: {s: [] for s in terminal_nodes} for n in nodes_tuple}

            # 5. Create igraph
            adj_matrix = (DC_vals != 0).astype(int)
            g = ig.Graph.Adjacency(adj_matrix.tolist(), mode="directed")

            # 6. Iterate and calculate
            for start_idx, start_node in tqdm(enumerate(nodes_tuple), desc="Outer Loop", total=len(nodes_tuple), disable=silent):
                # Edge case: start node is already a terminal
                if start_idx in terminal_indices:
                    paths_dict[start_node][start_node].append([start_node])
                    sppr_matrix[start_idx, term_idx_to_col[start_idx]] = 1.0
                    continue

                # Fetch all simple paths at C-speed (returns list of integer lists)
                # paths = g.get_all_simple_paths(start_idx, to=terminal_indices)
                paths = _get_paths_with_safety_valve(g, start_idx, terminal_indices)

                for path in tqdm(paths, disable=silent, desc="Inner Loop", leave=False):
                    sink_idx = path[-1]
                    sink_node = nodes_tuple[sink_idx]

                    # Production required along this path = product of edge weights A[u][v].
                    prod = 1.0
                    for i in range(len(path) - 1):
                        prod *= A_list[path[i]][path[i + 1]]

                    # Accumulate into the (source, basal-sink) cell -- summing over all paths.
                    # Add directly to the pre-allocated matrix using mapped indices
                    sppr_matrix[start_idx, term_idx_to_col[sink_idx]] += prod

                    # Construct the string path using our cached tuple
                    paths_dict[start_node][sink_node].append([nodes_tuple[v] for v in path])

            # 7. Wrap the matrix in a DataFrame instantly
            SPPR = pd.DataFrame(
                sppr_matrix,
                index=DC.index,
                columns=terminal_nodes
            )

            if use_EE:
                SPPR = SPPR.mul(self.EE, axis='index')

            return SPPR, A, paths_dict
        def _fast_EwE_no_paths(TE_option, use_EE, silent):
            """
            High-performance SPPR calculation using vectorized edge lookups
            and segmented products.
            """
            # 1. Data Preparation
            DC = self.get_DC(DET_as_PP=True)
            TE = self.get_TE(TE_option=TE_option, as_matrix=True, DET_values=1)

            # Pre-clean to avoid NaNs/Infs
            TE = TE.replace(0, np.nan).fillna(1.0)
            A = (DC / TE).replace([np.inf, -np.inf], 0).fillna(0)

            A_vals = A.values
            nodes = DC.index.tolist()
            num_nodes = len(nodes)

            # 2. Identify Terminal Nodes (Sinks)
            terminal_mask = (DC.values == 0).all(axis=1)
            terminal_indices = np.where(terminal_mask)[0]
            terminal_nodes = [nodes[i] for i in terminal_indices]
            term_idx_to_col = {t_idx: col for col, t_idx in enumerate(terminal_indices)}

            # 3. Graph Traversal (C-Backend)
            # Generate the graph structure for igraph
            adj_matrix = (DC.values != 0).astype(int)
            g = ig.Graph.Adjacency(adj_matrix.tolist(), mode="directed")

            # Fetch all simple paths as integer index lists
            # This is the fastest way to traverse millions of paths
            all_paths = []
            for start_idx in tqdm(range(num_nodes), disable=silent, desc="fetching paths"):
                # paths = g.get_all_simple_paths(start_idx, to=terminal_indices)
                paths = _get_paths_with_safety_valve(g, start_idx, terminal_indices)
                all_paths.extend(paths)

            if not all_paths:
                return pd.DataFrame(0.0, index=DC.index, columns=terminal_nodes), A

            # 4. Vectorized "Segmented" Math
            # We flatten all paths into edges to calculate products in bulk
            edge_starts = []
            edge_ends = []
            lengths = []
            metadata = []  # (start_node_idx, sink_node_idx)

            for p in tqdm(all_paths, disable=silent, desc="traversing paths"):
                metadata.append((p[0], p[-1]))
                if len(p) > 1:
                    edge_starts.extend(p[:-1])
                    edge_ends.extend(p[1:])
                    lengths.append(len(p) - 1)
                else:
                    lengths.append(0)

            # Convert to arrays for BLAS-speed operations
            edge_starts = np.array(edge_starts)
            edge_ends = np.array(edge_ends)
            lengths = np.array(lengths)

            # Bulk-fetch all edge weights from A at once
            all_edge_values = A_vals[edge_starts, edge_ends]

            # Calculate products using reduceat (C-level segmented product).
            # reduceat multiplies each contiguous segment of the flattened edge array, so each
            # path's edge weights collapse to a single product in one vectorized call.
            path_products = np.ones(len(all_paths))
            if len(all_edge_values) > 0:
                # Calculate the starting index of each path in the flattened edge array
                indices = np.zeros(len(lengths), dtype=int)
                indices[1:] = np.cumsum(lengths)[:-1]

                # Only process paths that actually have edges (length > 0)
                valid_mask = lengths > 0
                path_products[valid_mask] = np.multiply.reduceat(all_edge_values, indices[valid_mask])

            # 5. Aggregate into final SPPR Matrix
            sppr_matrix = np.zeros((num_nodes, len(terminal_indices)))

            # Map metadata to coordinate indices
            meta_arr = np.array(metadata)
            rows = meta_arr[:, 0]
            cols = np.array([term_idx_to_col[tid] for tid in meta_arr[:, 1]])

            # Vectorized 'Scatter-Add': sum every path's product into its (source, sink) cell.
            np.add.at(sppr_matrix, (rows, cols), path_products)

            # 6. Final DataFrame formatting
            SPPR = pd.DataFrame(sppr_matrix, index=DC.index, columns=terminal_nodes)

            if use_EE:
                SPPR = SPPR.mul(self.EE, axis='index')

            return SPPR, A, {}

        if return_paths:
            return _slow_EwE_with_paths(TE_option=TE_option, use_EE=use_EE, silent=silent)
        else:
            return _fast_EwE_no_paths(TE_option=TE_option, use_EE=use_EE, silent=silent)

    def SPPR_EwE_Ulanowicz(self, TE_option: str, global_TE: str | float = 'mean', use_EE: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Matrix (nullspace) reformulation of the EwE path-summation SPPR.

        Instead of enumerating paths, build the per-edge weight matrix A = DC/TE with cycles
        removed, replace each basal (producer) row by an identity row, and find the steady-state
        SPPR as the nullspace of L = A - I (so A.x = x). The nullspace basis is RREF-normalized
        so each resulting column is anchored to exactly one basal source. When use_EE the rows
        are scaled by ecotrophic efficiency.

        Args:
            TE_option (str): transfer-efficiency mode for building A; one of 'GE', 'TE',
                'With Egestion', 'global' (the 'global' case requires global_TE).
            global_TE (str | float, optional): the global TE used when TE_option == 'global';
                'mean' (catch-weighted, or biomass-weighted when total catch is 0) or a literal
                float. Defaults to 'mean'.
            use_EE (bool, optional): if True, scale each row by ecotrophic efficiency EE.
                Defaults to True.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: (SPPR, A, L), the per-group SPPR
            (groups x basal sources), the per-edge weight matrix A, and L = A - I.

        Raises:
            ValueError: if L has an empty nullspace (no steady-state solution / disconnected).
        """
        DC = self.get_DC(DET_as_PP=True)
        DCNoCyc = remove_cycles(DC.copy(), new=False)  # should work on Z instead?
        # zero values of DC where DCNoCyc is 0:
        DC[DCNoCyc == 0] = 0
        # DC = DC.div(DC.sum(axis=1), axis=0).fillna(0)

        # according to 2015's article, EwE uses TE = GE*EE. in to EwE user guide, they use just GE.
        TE = self.get_TE(TE_option=TE_option, global_TE=global_TE, as_matrix=True,  DET_values=1)
        
        A = (DC / TE).fillna(0).values
        A[TE.values == 0] = 0
        for i in range(len(DC)):# Replace producer rows with identity rows
            if DC.values[i, :].sum() == 0:
                A[i, :] = sm.zeros(1, len(DC))
                A[i, i] = 1
        A = pd.DataFrame(A, index=DC.index, columns=DC.columns)
        A, new_index, new_columns = move_scattered_identity(A)
        A = mat_from_np(A)

        # calculate L matrix
        L = A - sm.eye(A.rows)
        
        # solve for SPPR
        ns = L.nullspace()

        if len(ns) == 0:
            raise ValueError("No steady-state solution found. Check matrix connectivity.")
        
        # Stack the nullspace vectors as rows and RREF them so each basis row is pivoted on a
        # single basal source (a leading 1), giving one SPPR column per basal source.
        ns = [np.array(s).T.flatten() for s in ns]
        M = sm.Matrix(ns)
        rref_matrix, _ = M.rref()
        # Drop all-zero rows and convert the symbolic RREF rows to float; transpose so columns
        # are the basal sources and rows are the groups.
        basis = [np.array(rref_matrix.row(i).evalf()).astype(float).flatten()
                    for i in range(rref_matrix.rows)
                    if not rref_matrix.row(i).is_zero]
        basis = np.array(basis, dtype=float).T

        # turn back to DataFrames:
        SPPR = pd.DataFrame(basis, index=new_index)
        # Name each column by the group whose row holds the pivot 1 (its anchoring basal source).
        SPPR = SPPR.rename(columns=lambda c: SPPR.index.values[SPPR[c] == 1][0])
        if use_EE:
            SPPR = SPPR.mul(self.EE, axis='index')
        A = pd.DataFrame(np.array(A.tolist(), dtype=float), index=new_index, columns=new_columns)
        L = pd.DataFrame(np.array(L.tolist(), dtype=float), index=new_index, columns=new_columns)

        return SPPR, A, L

    def SPPR_2015(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """2015-method SPPR: a matrix-inversion (Leontief-style) formulation.

        Detritus columns are dissolved by reassigning the PP-derived fraction of each detritus
        flow back onto the PP groups, leaving only living compartments. The production-normalized
        transaction matrix A then yields the production-requirement matrix L = (I - A)^-1, whose
        PP columns give the per-group SPPR; a balancing detritus SPPR is added back at the end.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: (SPPR, A, L), the per-group SPPR
            (groups x PP+Import sources), the production-normalized transaction matrix A, and
            the production-requirement matrix L = (I - A)^-1.
        """
        groups_data = self.get_groups_df()
        production = self.p.copy()
        PP_seq = list(self.get_Import_seq()) + list(self.get_PP_seq())
        DET_seqs = self.get_DET_seq()  # list, may have >1 elements
        Z = self.get_Z(DET_as_PP=False)

        # combine part of each DET that is PP into PP row:
        Z_without_DET = Z.copy()
        for det_j in DET_seqs:
            det_row_total = Z_without_DET.loc[det_j, :].sum()
            percent_of_det_that_is_PP = (Z_without_DET.loc[det_j, PP_seq] / det_row_total).fillna(0)
            for i in PP_seq: # this is what is implemented in the article
                Z_without_DET.loc[:, i] += percent_of_det_that_is_PP[i] * Z_without_DET.loc[:, det_j]
        Z_without_DET = Z_without_DET.drop(index=DET_seqs, columns=DET_seqs)

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
        L[seq_to_drop] = 0
        new_index = L.index.union(seq_to_drop)
        L = L.reindex(new_index).fillna(0)
        for d in seq_to_drop:  # self-basal diagonal only; block assignment leaked SPPR=1 into dropped PP_seq columns
            L.loc[d, d] = 1
        SPPR = L.loc[PP_seq, :].T

        # add back sppr_det that makes model balanced:
        SPPR = SPPR.reindex(self.catch.index, fill_value=0)
        det_seqs = self.get_DET_seq()
        if det_seqs:
            for i in PP_seq:
                SPPR.loc[det_seqs, i] = self.M0.loc[PP_seq][i] / (self.M0 + self.egestion).sum()

        return SPPR, A, L

    def _build_det_BC(self, sppr_basis: pd.DataFrame, non_DET_sppr: pd.Series, DET_seq: list, DC: pd.DataFrame, TE_option: str) -> tuple[np.ndarray, np.ndarray]:
        """Build the detritus recycling system (I - B) x = c for the GE / With Egestion modes.

        Each detritus group l has an unknown scaling factor x_l; recycled detritus production
        depends on the SPPR already attributed to other compartments, giving the coupled system:

            x_l     = scaling factor for detritus group l
            c[l]    = m_eff_l . non_DET_sppr       (production drawn from non-detritus sources)
            B[l, j] = m_eff_l . sppr_basis[det_j]  (recursive dependence on detritus basis j)

        where m_eff_l = M0*fracs/q_l                                 (GE)
                      = M0*fracs/q_l + DC.T @ (egestion*fracs/q_l)   (With Egestion)
        and fracs = det_fate[:, det_l] is the fraction of each group's flow_to_det reaching
        det_l. In a multi-DET model, a present det_fate matrix lacking a column for det_l means
        nothing feeds det_l -> fracs = 0 (a warning is emitted), NOT whole-flow. For single-DET
        (or no det_fate at all) fracs defaults to 1 (legacy behavior).

        Args:
            sppr_basis (pd.DataFrame): the current SPPR basis (per group x basal source) whose
                detritus columns the recycling depends on.
            non_DET_sppr (pd.Series): per-group SPPR contribution from all non-detritus sources.
            DET_seq (list): seq IDs of the detritus groups being solved.
            DC (pd.DataFrame): diet-composition matrix (used for the egestion term).
            TE_option (str): 'GE' or 'With Egestion'; selects whether the egestion term is added.

        Returns:
            tuple[np.ndarray, np.ndarray]: (B, c), the k x k recycling matrix and the length-k
            source vector of the system (I - B) x = c, where k = len(DET_seq).
        """
        k = len(DET_seq)
        det_fate = getattr(self, '_det_fate', None)
        B = np.zeros((k, k))
        c = np.zeros(k)
        for li, det_l in enumerate(DET_seq):
            q_l = self.q[det_l] if self.q[det_l] > 0 else 1.0
            if det_fate is not None and det_l in det_fate.columns:
                fracs = det_fate[det_l].reindex(self.M0.index).fillna(0)
            elif det_fate is not None and k > 1:
                warnings.warn(
                    f"det_fate has no column for DET group {det_l} "
                    f"({self.seq2name.get(det_l)}); treating its detritus inflow as 0.",
                    RuntimeWarning,
                )
                fracs = pd.Series(0.0, index=self.M0.index)
            else:
                fracs = pd.Series(1.0, index=self.M0.index)
            m_l = (self.M0 * fracs / q_l).fillna(0)
            if TE_option == 'With Egestion':
                e_l = (self.egestion * fracs / q_l).fillna(0)
                m_eff_l = m_l + DC.T @ e_l
            else:  # GE
                m_eff_l = m_l
            # Align operands to M0.index by label (0-filling groups absent from the basis).
            # No-op for SPPR_new / as_PP (basis already spans all groups); for the as_DC
            # symbolic basis it drops the Import rows cleanly instead of raising on misalignment.
            m_eff_l = m_eff_l.reindex(self.M0.index).fillna(0)
            nds = non_DET_sppr.reindex(self.M0.index).fillna(0)
            c[li] = float(m_eff_l @ nds)
            for ji, det_j in enumerate(DET_seq):
                basis_j = sppr_basis[det_j].reindex(self.M0.index).fillna(0)
                B[li, ji] = float(m_eff_l @ basis_j)
        return B, c

    @staticmethod
    def _spectral_radius(M: np.ndarray) -> float:
        """Return the spectral radius (largest absolute eigenvalue) of M.

        Used to test whether the detritus recycling matrix B is subcritical: when rho(B) < 1 a
        finite, nonnegative solution to the recycling system (I - B) x = c exists, and the
        Neumann series sum(B^k) converges. rho >= 1 signals divergent recycling.

        Args:
            M (np.ndarray): a square matrix (the detritus recycling matrix B).

        Returns:
            float: the largest absolute eigenvalue of M, or 0.0 if M is empty.
        """
        M = np.asarray(M, dtype=float)
        if M.size == 0:
            return 0.0
        return float(np.max(np.abs(np.linalg.eigvals(M))))

    def _resolve_det_param(self, param: float | dict | None, DET_seq: list, default: float) -> np.ndarray:
        """Resolve a per-DET parameter into a float array aligned with DET_seq.

        Turns the user-facing det_theta / det_external_sppr knobs into per-DET vectors aligned
        with the detritus columns being scaled.

        Args:
            param (float | dict | None): a scalar broadcast to all DET groups, or a dict keyed
                by DET group seq (int) or DET group name (str). None falls back to `default`.
                Dict entries missing for a given DET group also fall back to `default`.
            DET_seq (list): seq IDs of the detritus groups, defining the output order.
            default (float): value used for None param and for missing dict keys.

        Returns:
            np.ndarray: a length-len(DET_seq) float array of per-DET values.
        """
        if param is None:
            param = default
        if np.isscalar(param):
            return np.full(len(DET_seq), float(param), dtype=float)
        out = np.full(len(DET_seq), float(default), dtype=float)
        for i, d in enumerate(DET_seq):
            if d in param:
                out[i] = float(param[d])
            elif self.seq2name.get(d) in param:
                out[i] = float(param[self.seq2name[d]])
        return out

    @staticmethod
    def _collapse_det_scaling(SPPR: pd.DataFrame, DET_seq: list, non_DET_sppr: pd.Series, M0: pd.Series, egestion: pd.Series, q: pd.Series,
                              flow_to_det: pd.Series, DC: pd.DataFrame, TE_option: str, det_fate: Optional[pd.DataFrame] = None,
                              theta: Optional[np.ndarray] = None, ext: Optional[np.ndarray] = None, det_open_mode: str = 'none') -> tuple[pd.DataFrame, float, float, float]:
        """Fallback DET scaling: pool all detritus into one compartment and solve a 1-D problem.

        Used when the coupled detritus system is unstable. Sums the individual detritus basis
        vectors into one combined pool, solves the scalar self-consistency x = a + b*x, and
        multiplies every DET column of SPPR by that single factor. The combined denominator
        q_combined >> any individual q_j, which usually brings b well below 1.

        M0 / egestion are weighted by fate_to_modeled = sum over the modeled DET columns of
        det_fate, so the numerator counts only material actually entering the modeled detritus
        system (consistent with the fate-weighted q_combined). With det_fate rows summing to 1
        this weight is 1 (a no-op on the real models). Openness is applied to the pooled scalar
        using the mean theta / ext across DET groups.

        Args:
            SPPR (pd.DataFrame): the SPPR basis; its DET columns are scaled in place.
            DET_seq (list): seq IDs of the detritus groups.
            non_DET_sppr (pd.Series): per-group SPPR from non-detritus sources.
            M0 (pd.Series): natural-mortality flow per group.
            egestion (pd.Series): egestion flow per group.
            q (pd.Series): consumption per group (used to form q_combined).
            flow_to_det (pd.Series): M0 + egestion per group (fallback denominator).
            DC (pd.DataFrame): diet-composition matrix (egestion term, With Egestion only).
            TE_option (str): 'GE' or 'With Egestion'.
            det_fate (Optional[pd.DataFrame]): per-group fate fractions to each DET column.
                Defaults to None.
            theta (Optional[np.ndarray]): per-DET availability/retention used for openness;
                its mean is applied. Defaults to None.
            ext (Optional[np.ndarray]): per-DET external SPPR used for source dilution; its mean
                is applied. Defaults to None.
            det_open_mode (str): 'none', 'recycling_loss' (damp b), or 'source_dilution' (damp b
                and dilute a toward ext). Defaults to 'none'.

        Returns:
            tuple[pd.DataFrame, float, float, float]: (SPPR, sppr_det, a, b), the scaled SPPR,
            the pooled scaling scalar sppr_det = a/(1-b), and the solved coefficients a and b.

        Raises:
            ValueError: if b >= 1 (pooled recycling still diverges).
        """
        q_combined = float(sum(float(q[d]) for d in DET_seq if float(q[d]) > 0))
        if q_combined <= 0:
            q_combined = float(flow_to_det.sum())

        if det_fate is not None:
            fate_to_modeled = (det_fate.reindex(index=M0.index, columns=list(DET_seq))
                               .fillna(0).sum(axis=1).clip(lower=0.0, upper=1.0))
        else:
            fate_to_modeled = pd.Series(1.0, index=M0.index)

        m_combined = (M0 * fate_to_modeled / q_combined).fillna(0)

        if TE_option == 'With Egestion':
            e_combined = (egestion * fate_to_modeled / q_combined).fillna(0)
            m_eff = m_combined + DC.T @ e_combined
        else:  # GE
            m_eff = m_combined

        SPPR_combined_raw = SPPR[list(DET_seq)].sum(axis=1)

        a = float(m_eff @ non_DET_sppr)
        b = float(m_eff @ SPPR_combined_raw)

        if det_open_mode != 'none' and theta is not None:
            th = float(np.mean(theta))
            if det_open_mode == 'recycling_loss':
                b = th * b
            elif det_open_mode == 'source_dilution':
                ex = float(np.mean(ext)) if ext is not None else 0.0
                a = th * a + ex * (1.0 - th)
                b = th * b

        if b >= 1.0:
            raise ValueError(
                f"pooled-detritus fallback also diverges (b={b:.4f} >= 1). "
                "Detrital cycling is too strong for SPPR_new."
            )

        sppr_det = a / (1.0 - b)
        for det_j in DET_seq:
            SPPR[det_j] *= sppr_det
        return SPPR, sppr_det, a, b

    def _solve_det_scaling(self, B: np.ndarray, c_vec: np.ndarray, DET_seq: list, SPPR: pd.DataFrame, non_DET_sppr: pd.Series, DC: pd.DataFrame, TE_option: str,
                           det_collapse_mode: str = 'never', det_open_mode: str = 'none',
                           det_theta: float | dict = 1.0, det_external_sppr: float | dict = 0.0,
                           tol: float = 1e-10, cond_threshold: float = 1e10) -> pd.DataFrame:
        """Resolve the detritus recycling system: apply openness, choose solve-vs-pool, scale SPPR.

        Applies the openness transform to (B, c), tests stability via spectral radius and matrix
        conditioning, decides whether to solve the coupled system directly or fall back to the
        pooled scaling, scales the DET columns of SPPR in place, and records a diagnostic
        dict on self.detritus_resolution_info.

        Args:
            B (np.ndarray): k x k recycling matrix from _build_det_BC.
            c_vec (np.ndarray): length-k source vector from _build_det_BC.
            DET_seq (list): seq IDs of the detritus groups (length k).
            SPPR (pd.DataFrame): SPPR basis whose DET columns are scaled in place.
            non_DET_sppr (pd.Series): per-group SPPR from non-detritus sources (for the pool).
            DC (pd.DataFrame): diet-composition matrix (egestion term in the pool path).
            TE_option (str): 'GE' or 'With Egestion'.
            det_collapse_mode (str): how the solve-vs-pool decision is made:
                'never' -> always solve the coupled system directly (may return negative SPPR
                but never raises; the Monte-Carlo samplers rely on this to reject unstable draws);
                'auto' -> pool only if the system is unstable (spectral radius rho >= 1-tol or
                condition number > cond_threshold);
                'always' -> always use the pooled scaling. Defaults to 'never'.
            det_open_mode (str): openness model: 'none' (closed recycling); 'recycling_loss'
                (a fraction of recycled detritus is lost, damping B by diag(theta)); or
                'source_dilution' (damp B by diag(theta) AND dilute the source toward an external
                SPPR). Defaults to 'none'.
            det_theta (float | dict): detritus availability/retention fraction; a single float
                for all DET groups or a dict keyed by DET seq (int) or name (str). 1.0 reproduces
                the closed system. Defaults to 1.0.
            det_external_sppr (float | dict): external SPPR assigned to diluted material under
                'source_dilution'; same scalar-or-dict form as det_theta. Defaults to 0.0.
            tol (float): tolerance used in the rho >= 1 - tol instability test. Defaults to 1e-10.
            cond_threshold (float): condition-number threshold above which (I - B) is deemed
                ill-conditioned. Defaults to 1e10.

        Returns:
            pd.DataFrame: SPPR with its DET columns scaled.

        Raises:
            ValueError: if det_open_mode or det_collapse_mode is not a recognized value, or if
                the pooled fallback also diverges.

        Notes:
            Openness transform on the recycling system (theta, ext aligned to DET_seq):
                none            : B_open = B,             c_open = c
                recycling_loss  : B_open = diag(theta) B, c_open = c
                source_dilution : B_open = diag(theta) B, c_open = theta*c + ext*(1-theta)
            With theta=1 / ext=0 all three reduce to the original closed system.
        """
        k = len(DET_seq)
        theta = self._resolve_det_param(det_theta, DET_seq, 1.0)
        ext = self._resolve_det_param(det_external_sppr, DET_seq, 0.0)

        if det_open_mode == 'none':
            B_open, c_open = B.copy(), c_vec.copy()
        elif det_open_mode == 'recycling_loss':
            B_open, c_open = np.diag(theta) @ B, c_vec.copy()
        elif det_open_mode == 'source_dilution':
            B_open = np.diag(theta) @ B
            c_open = theta * c_vec + ext * (1.0 - theta)
        else:
            raise ValueError("det_open_mode must be 'none', 'recycling_loss', or 'source_dilution'")

        IminusB = np.eye(k) - B_open
        rho = self._spectral_radius(B_open)
        try:
            cond = float(np.linalg.cond(IminusB))
        except np.linalg.LinAlgError:
            cond = np.inf

        if det_collapse_mode == 'always':
            use_collapse, reason = True, 'mode_always'
        elif det_collapse_mode == 'auto':
            if rho >= 1.0 - tol:
                use_collapse, reason = True, 'spectral_radius_ge_1'
            elif cond > cond_threshold:
                use_collapse, reason = True, 'ill_conditioned'
            else:
                use_collapse, reason = False, None
        elif det_collapse_mode == 'never':
            use_collapse, reason = False, None
        else:
            raise ValueError("det_collapse_mode must be 'never', 'auto', or 'always'")

        det_fate = getattr(self, '_det_fate', None)
        flow_to_det = (self.M0 + self.egestion).fillna(0)

        if use_collapse:
            SPPR, scalar, a_p, b_p = self._collapse_det_scaling(
                SPPR, DET_seq, non_DET_sppr, self.M0, self.egestion, self.q,
                flow_to_det, DC, TE_option, det_fate=det_fate,
                theta=theta, ext=ext, det_open_mode=det_open_mode)
            self.detritus_resolution_info = {
                'method': 'pooled_detritus_scaling', 'reason': reason,
                'rho_B': rho, 'cond_IminusB': cond,
                'det_seq': list(DET_seq), 'det_names': [self.seq2name[d] for d in DET_seq],
                'open_mode': det_open_mode, 'theta': theta.tolist(),
                'collapse_scalar': scalar, 'collapse_a': a_p, 'collapse_b': b_p,
            }
        else:
            # Directly solve the coupled recycling system (I - B) x = c for the per-DET scaling
            # factors, then multiply each detritus column by its own factor.
            x_vec = np.linalg.solve(IminusB, c_open)
            for i, det_j in enumerate(DET_seq):
                SPPR[det_j] *= x_vec[i]
            self.detritus_resolution_info = {
                'method': 'multi_detritus' if k > 1 else 'single_detritus', 'reason': reason,
                'rho_B': rho, 'cond_IminusB': cond,
                'det_seq': list(DET_seq), 'det_names': [self.seq2name[d] for d in DET_seq],
                'open_mode': det_open_mode, 'theta': theta.tolist(),
                'x_vec': x_vec.tolist(),
            }
        return SPPR

    def SPPR_new(self, TE: Optional[pd.DataFrame] = None, TE_option: str = 'GE', DET_TE_vals: float = 1,
                 det_collapse_mode: str = 'never', det_open_mode: str = 'none',
                 det_theta: float | dict = 1.0, det_external_sppr: float | dict = 0.0,
                 fix_EE_0_cases: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Primary numeric SPPR solver via the nullspace of L = A - I.

        Builds the per-edge weight matrix A = DC/TE, replaces basal (producer) rows with
        identity rows, and finds the steady-state SPPR as the nullspace of L = A - I, RREF-
        normalized so each column is anchored to one basal source. Detritus columns are then
        resolved by TE_option: 'TE' scales each detritus column by its direct PP+Import inflow
        share (times theta); 'GE' / 'With Egestion' build the coupled recycling system
        (I - B) x = c via _build_det_BC and resolve it (with optional openness / collapse)
        through _solve_det_scaling.

        Args:
            TE (Optional[pd.DataFrame]): an explicit TE matrix (e.g. a Monte-Carlo sample);
                if None it is built from TE_option. Defaults to None.
            TE_option (str): transfer-efficiency mode when TE is None; one of 'GE', 'TE',
                'With Egestion', 'global'. Defaults to 'GE'.
            DET_TE_vals (float): TE assigned to detritus rows when building the TE matrix.
                Defaults to 1.
            det_collapse_mode (str): detritus solve-vs-pool strategy: 'never' (always solve the
                coupled system; may return negative SPPR but never raises), 'auto' (pool only if
                unstable -- spectral radius >= 1 or ill-conditioned), or 'always' (always pool).
                Defaults to 'never'.
            det_open_mode (str): detritus recycling openness: 'none' (closed recycling),
                'recycling_loss' (damp recycling B by diag(theta)), or 'source_dilution' (damp B
                by diag(theta) and dilute the source toward an external SPPR). Defaults to 'none'.
            det_theta (float | dict): detritus availability/retention fraction; a float applied
                to all detritus groups or a dict keyed by DET seq (int) or name (str). 1.0
                reproduces the closed system. Defaults to 1.0.
            det_external_sppr (float | dict): external SPPR assigned to diluted material under
                'source_dilution'; same scalar-or-dict form as det_theta. Defaults to 0.0.
            fix_EE_0_cases (bool): if True (default), under TE_option='TE' re-credit the PP
                consumed by EE=0 dead-end groups (transfer efficiency te=(p-M0)/q == 0 because
                M0==p, so they are zeroed out of the nullspace) back to the detritus pool, since
                their production is 100% other-mortality and physically flows to detritus. This
                closes the global PP balance (inflow == outflow) that those severed groups would
                otherwise break. Only active for single-DET models under 'TE'; a no-op when no
                EE=0 groups exist. A RuntimeWarning is emitted whenever EE=0 groups are present.
                Does NOT address near-singular (0 < EE << 1) groups, which are an inherent
                singularity of the TE method. Defaults to True.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: (SPPR, A, L), the per-group SPPR
            (groups x basal sources), the per-edge weight matrix A, and L = A - I.

        Raises:
            ValueError: if L has an empty nullspace (no steady-state solution).
            Exception: if TE_option is not one of the supported strings.
        """
        # get DC:
        DC = self.get_DC(DET_as_PP=True, normalize=False)
        
        # get TE matrix:
        if TE is not None:
            GE = TE.copy()
        else:
            GE = self.get_TE(TE_option=TE_option, DET_values=DET_TE_vals, as_matrix=True)

        # calculate A matrix and turn to symbolic matrix:
        A = (DC / GE).fillna(0).values.copy()
        A[GE.values == 0] = 0
        for i in range(len(DC)):# Replace producer rows with identity rows
            if DC.values[i, :].sum() == 0:
                A[i, :] = sm.zeros(1, len(DC))
                A[i, i] = 1
        A = pd.DataFrame(A, index=DC.index, columns=DC.columns)
        A, new_index, new_columns = move_scattered_identity(A)
        A = mat_from_np(A)

        # calculate L matrix
        L = A - sm.eye(A.rows)
        
        # solve for SPPR
        ns = L.nullspace()

        if len(ns) == 0:
            raise ValueError("No steady-state solution found. Check matrix connectivity.")
        
        # Stack the nullspace vectors as rows and RREF them so each basis row is pivoted on a
        # single basal source (a leading 1), yielding one SPPR column per basal source.
        ns = [np.array(s).T.flatten() for s in ns]
        M = sm.Matrix(ns)
        rref_matrix, _ = M.rref()
        # Drop all-zero rows, cast the symbolic RREF rows to float, and transpose so columns are
        # the basal sources and rows are the groups.
        basis = [np.array(rref_matrix.row(i).evalf()).astype(float).flatten()
                    for i in range(rref_matrix.rows)
                    if not rref_matrix.row(i).is_zero]
        basis = np.array(basis, dtype=float).T

        # turn back to DataFrames:
        SPPR = pd.DataFrame(basis, index=new_index)
        # Name each column by the group holding the pivot 1 (its anchoring basal source).
        SPPR = SPPR.rename(columns=lambda c: SPPR.index.values[SPPR[c] == 1][0])
        A = pd.DataFrame(np.array(A.tolist(), dtype=float), index=new_index, columns=new_columns)
        L = pd.DataFrame(np.array(L.tolist(), dtype=float), index=new_index, columns=new_columns)

        # find sppr_det if it is in the output:
        DET_seq = self.get_DET_seq()
        det_fate = getattr(self, '_det_fate', None)
        flow_to_det = (self.M0 + self.egestion).fillna(0)
        non_DET_sppr = SPPR.drop(columns=DET_seq).sum(axis=1)

        # Unified detritus resolution for single- and multi-DET models.
        # GE / With Egestion build the coupled recycling system (I - B) x = c via the shared
        # _build_det_BC helper and resolve it through _solve_det_scaling (which applies the
        # openness transform, decides solve-vs-pool, and scales the DET columns in place).
        # TE has no recycling matrix, so each DET column is scaled by its direct PP+Import
        # inflow share, with availability theta applied as a plain multiplier.
        if TE_option == 'TE':
            PP_Import_seq = list(self.get_PP_seq() + self.get_Import_seq())
            theta = self._resolve_det_param(det_theta, DET_seq, 1.0)
            # EE=0 dead-ends: consumer groups whose entire TE row is 0 (te=(p-M0)/q == 0 because
            # M0==p) get zeroed out of the nullspace by A[GE==0]=0, severing their mortality->
            # detritus recycling. The PP they consume is then neither transferred up (nothing
            # preys on them) nor returned to detritus, so the global PP budget leaks. Detect them
            # from the actual TE matrix in use (so a Monte-Carlo TE sample is handled too).
            deadend_mask = (GE == 0).all(axis=1)
            deadend_seq = [g for g in deadend_mask.index[deadend_mask]
                           if g not in PP_Import_seq and g not in list(DET_seq)]
            if deadend_seq:
                names = ", ".join(f"{g} ({self.seq2name.get(g, '?')})" for g in deadend_seq)
                if fix_EE_0_cases and len(DET_seq) == 1:
                    warnings.warn(
                        f"{len(deadend_seq)} group(s) have TE=0 because EE=0 (M0=p): {names}. "
                        f"Their consumed PP is re-credited to detritus (fix_EE_0_cases=True) so "
                        f"the TE SPPR stays PP-balanced.",
                        RuntimeWarning,
                    )
                else:
                    reason = (" (multi-DET: re-credit not applied)" if len(DET_seq) > 1
                              else "; pass fix_EE_0_cases=True to correct")
                    warnings.warn(
                        f"{len(deadend_seq)} group(s) have TE=0 because EE=0 (M0=p): {names}. "
                        f"Their consumed PP is dropped, so the TE SPPR will not balance{reason}.",
                        RuntimeWarning,
                    )
            # Near-singular groups (0 < EE << 1, i.e. te just above 0) are NOT zeroed, but their
            # SPPR ~ 1/te blows up (an inherent singularity of the TE method). They are not fixed
            # here; worse, when an EE=0 dead-end feeds on one, the re-credit's detritus-coupling
            # term (b) inflates and the resulting detritus SPPR can exceed 1 implausibly. Warn so
            # an inflated sppr_det is explained rather than silent.
            te_row = GE.iloc[:, 0] if GE.shape[1] else pd.Series(dtype=float)
            near_singular = [g for g in te_row.index
                             if g not in PP_Import_seq and g not in list(DET_seq)
                             and 0 < abs(float(te_row[g])) < 1e-3]
            if near_singular:
                names = ", ".join(f"{g} ({self.seq2name.get(g, '?')})" for g in near_singular)
                warnings.warn(
                    f"{len(near_singular)} group(s) have near-zero TE (0 < EE << 1): {names}. "
                    f"Their SPPR is near-singular (inherent to the TE method); the detritus SPPR "
                    f"may be inflated and is not corrected by fix_EE_0_cases.",
                    RuntimeWarning,
                )
            # Re-credit only makes sense for a single detritus pool (multi-DET pools couple
            # through det_fate and are out of scope here -- see MULTIDET_TE_DISCREPANCY.md).
            apply_fix = fix_EE_0_cases and bool(deadend_seq) and len(DET_seq) == 1
            for i, det_l in enumerate(DET_seq):
                q_l = self.q[det_l] if self.q[det_l] > 0 else flow_to_det.sum()
                if det_fate is not None and det_l in det_fate.columns:
                    fr = det_fate[det_l].reindex(flow_to_det.index).fillna(0)
                    inflow = (flow_to_det[PP_Import_seq] * fr[PP_Import_seq]).sum()
                elif det_fate is not None and len(DET_seq) > 1:
                    # Point 6: a present det_fate missing this DET column means nothing feeds it.
                    inflow = 0.0
                else:
                    inflow = flow_to_det[PP_Import_seq].sum()
                if apply_fix:
                    # Each dead-end g sends its full production (=M0, since M0=p) to detritus,
                    # carrying the PP it consumed: q_g * (DC[g] . sppr_tot), where
                    #   sppr_tot = non_DET_sppr + m * SPPR[:, det_l]   (m = the DET multiplier).
                    # Splitting the diet dot-product into its PP/Import part (a) and its
                    # detritus part (b) makes the detritus inflow linear in m, with the
                    # availability damping theta applied to the final pool SPPR:
                    #   m = theta * (inflow + a) / (q_l - theta * b).
                    a = sum(self.q[g] * (DC.loc[g] * non_DET_sppr).sum() for g in deadend_seq)
                    b = sum(self.q[g] * (DC.loc[g] * SPPR[det_l]).sum() for g in deadend_seq)
                    m = theta[i] * (inflow + a) / (q_l - theta[i] * b)
                else:
                    m = (inflow / q_l) * theta[i]
                SPPR[det_l] *= m
        elif TE_option in ('GE', 'With Egestion'):
            B, c_vec = self._build_det_BC(SPPR, non_DET_sppr, DET_seq, DC, TE_option)
            SPPR = self._solve_det_scaling(
                B, c_vec, DET_seq, SPPR, non_DET_sppr, DC, TE_option,
                det_collapse_mode=det_collapse_mode, det_open_mode=det_open_mode,
                det_theta=det_theta, det_external_sppr=det_external_sppr)
        else:
            raise Exception("TE_option should be in ['GE', 'TE', 'With Egestion', 'global']")

        return SPPR, A, L
    
    def _SPPR_symbolic_helper_diet_import_as_PP(self, TE: Optional[pd.DataFrame], TE_option: str, DET_TE_vals: float, sppr_det_value: Optional[float],
                                                det_collapse_mode: str = 'never', det_open_mode: str = 'none',
                                                det_theta: float | dict = 1.0, det_external_sppr: float | dict = 0.0,
                                                fix_EE_0_cases: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, list, list]:
        """Symbolic SPPR helper, "diet import as PP" variant.

        Imported diet is treated like an extra primary-production source (its own free SPPR
        symbol fixed to 1). Builds the per-group SPPR equations A x - x = 0 symbolically, solves
        the non-detritus block, then the detritus block, and returns both the symbolic solution
        and a numeric basis matrix (sppr_mat). Detritus columns are scaled either by the exact
        symbolic solution (default path: det_open_mode='none' and det_collapse_mode='never') or
        through the shared numeric recycling solver (_build_det_BC + _solve_det_scaling) when
        openness/collapse is requested. The detritus diet rows are built to mirror SPPR_new's
        open-system detritus handling (det_fate-weighted inflow normalized by the routed q[det]),
        so this helper reproduces SPPR_new for matching TE_options.

        Args:
            TE (Optional[pd.DataFrame]): explicit TE matrix; if None it is built from TE_option.
            TE_option (str): transfer-efficiency / detritus-DC mode; one of 'GE', 'TE',
                'With Egestion', 'global'.
            DET_TE_vals (float): TE assigned to detritus rows when building the TE matrix.
            sppr_det_value (Optional[float]): if not None, every detritus column is scaled by
                this fixed value instead of being solved.
            det_collapse_mode (str): 'never', 'auto', or 'always' (see SPPR_new). Defaults to 'never'.
            det_open_mode (str): 'none', 'recycling_loss', or 'source_dilution' (see SPPR_new).
                Defaults to 'none'.
            det_theta (float | dict): detritus availability/retention; float or dict keyed by
                DET seq/name. Defaults to 1.0.
            det_external_sppr (float | dict): external SPPR for 'source_dilution'. Defaults to 0.0.
            fix_EE_0_cases (bool): mirrors SPPR_new's flag. When True (default) and TE_option='TE',
                single-DET EE=0 dead-end groups have their consumed PP re-credited to the detritus
                pool (the same (inflow+a)/(q_l-theta*b) factor SPPR_new uses); when False the plain
                inflow/q_l scalar is used. No effect for GE / With Egestion or for multi-DET.
                Defaults to True.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, list, list]: (sppr_symbolic, sppr_mat, equations,
            variables), where sppr_symbolic is the symbolic per-group solution, sppr_mat is the
            numeric basis matrix (groups x basal sources), equations is the full symbolic system,
            and variables is the ordered list of SPPR symbols.

        Raises:
            Exception: if TE_option is not one of the supported strings.
        """
        DET_seq = self.get_DET_seq()
        non_DET_seq = [i for i in self.GE.index if i not in DET_seq]
        Regular_seq = self.get_Regular_seq()
        Import_seq = self.get_Import_seq()
        PP_seq = self.get_PP_seq()

        # get DC and TE matrix:
        DC = self.get_DC(DET_as_PP=False, normalize=False)
        GE = self.get_TE(TE_option=TE_option, DET_values=DET_TE_vals, as_matrix=True)
        flow2det = (self.M0 + self.egestion).sum()
        if TE_option == 'TE':
            # Factor-0 TE semantics (see MULTIDET_TE_DISCREPANCY.md): detritus pools are basal
            # sources, exactly as SPPR_new treats them (DET_as_PP=True). Zero the WHOLE detritus
            # diet row so every DET pool becomes a free basal symbol; its PP content is assigned
            # afterwards by SPPR_new's direct PP+Import inflow scalar (the TE branch of the
            # DET-column scaling below). Zeroing only the Regular columns instead would leave a
            # consumer-fed (secondary) pool as a free source scaled by 1 -- i.e. treating
            # secondary detritus as 100% primary production -- which is the multi-DET TE
            # discrepancy this fixes.
            DC.loc[DET_seq, :] = 0
        elif TE_option == 'GE':
            det_fate_mat = getattr(self, '_det_fate', None)
            for det_j in DET_seq:
                q_j = self.q[det_j] if self.q[det_j] > 0 else flow2det / max(len(DET_seq), 1)
                if det_fate_mat is not None and det_j in det_fate_mat.columns and q_j > 0:
                    fracs = det_fate_mat[det_j].reindex(self.M0.index).fillna(0)
                    DC.loc[det_j, :] = ((self.M0 * fracs) / q_j).fillna(0)
                else:
                    DC.loc[det_j, :] = ((self.M0) / flow2det).fillna(0)
        elif TE_option == 'With Egestion':
            m = (self.M0 / flow2det).fillna(0)
            e = (self.egestion / flow2det).fillna(0)
            for det_j in DET_seq:
                DC.loc[det_j, :] = 0
            for det_j in DET_seq:
                DC.loc[det_j, :] = (m + e @ DC)
        else:
            raise Exception("diet_import_option should be in one of ['GE', 'TE', 'With Egestion', 'global']")

        if TE is not None:
            GE = TE.copy()

        # calculate A matrix and turn to symbolic matrix:
        A = (DC / GE).fillna(0).values.copy()
        A[GE.values == 0] = 0
        for i in range(len(DC)):  # Replace producer rows with identity rows
            if DC.values[i, :].sum() == 0:
                A[i, :] = 0
                A[i, i] = 1
        A = pd.DataFrame(A, index=DC.index, columns=DC.columns)

        # define and order symbols:
        # solve linear equation:
        index = Regular_seq + DET_seq + Import_seq + PP_seq
        symbols_by_trophic_info = {
            'Regular': sm.symbols([f'SPPR_{self.seq2name[s]}'.replace(' ', '_') for s in Regular_seq]),
            'DET': sm.symbols([f'SPPR_{self.seq2name[s]}'.replace(' ', '_') for s in DET_seq]),
            'Import': sm.symbols([f'SPPR_{self.seq2name[s]}'.replace(' ', '_') for s in Import_seq]),
            'PP': sm.symbols([f'SPPR_{self.seq2name[s]}'.replace(' ', '_') for s in PP_seq]),
        }
        ordered_symbols = symbols_by_trophic_info['Regular'] + symbols_by_trophic_info['DET'] \
                        + symbols_by_trophic_info['Import'] + symbols_by_trophic_info['PP']
        sppr_vec = pd.DataFrame(ordered_symbols, index=index)
        # Steady-state SPPR equations: A x - x = 0. Split into a non-detritus block (solved
        # first) and a detritus block (solved after substituting the non-DET solution back in).
        equations = A @ sppr_vec - sppr_vec
        equation_DET = equations.loc[DET_seq]
        equations_non_DET = equations.loc[non_DET_seq]
        sol = sm.linsolve(equations_non_DET.squeeze().tolist(), ordered_symbols)
        sol_tuple1 = list(sol)[0]
        sol_dict1 = dict(zip(ordered_symbols, sol_tuple1))
        sppr = sppr_vec.replace(sol_dict1)

        # solve DET equation:
        equation_DET = equation_DET.apply(lambda col: col.map(lambda x: x.subs(sol_dict1) if hasattr(x, 'subs') else x))
        sol = sm.linsolve(list(equation_DET.values.ravel()), symbols_by_trophic_info['DET'])
        sol_tuple2 = list(sol)[0]
        sol_dict2 = dict(zip(symbols_by_trophic_info['DET'], sol_tuple2))
        sppr_det_symbol = sppr_vec.loc[DET_seq].values[0, 0]
        sppr_det = sol_dict2[sppr_det_symbol]

        # sub:
        target_free_seq = DET_seq + Import_seq + PP_seq
        target_free_symbols = symbols_by_trophic_info['DET'] + symbols_by_trophic_info['Import'] + symbols_by_trophic_info['PP']
        sppr = sppr_vec.replace(sol_dict1)
        subs_dict = {v: 1 for v in sppr_vec.squeeze().loc[target_free_seq]}
        sppr_symbolic = sppr.apply(lambda col: col.map(lambda x: x.evalf(subs=subs_dict) if hasattr(x, 'evalf') else x))

        sppr_mat, _ = sm.linear_eq_to_matrix(sol_tuple1, target_free_symbols)
        sppr_mat = sm.lambdify([], sppr_mat, 'numpy')() # Converts SymPy matrix to NumPy
        sppr_mat = pd.DataFrame(sppr_mat, index=index, columns=target_free_seq)
        if sppr_det_value is not None:
            sppr_mat[DET_seq] *= float(sppr_det_value)
        elif TE_option == 'TE':
            # Mirror SPPR_new's TE detritus branch EXACTLY (keep in sync with SPPR_new): the
            # detritus rows were zeroed above, so the sppr_mat DET columns are the raw basal
            # bases. Scale each by its direct PP+Import inflow share (inflow/q_l)*theta -- a
            # secondary pool with no direct PP/Import inflow gets factor 0. Single-DET EE=0
            # dead-ends get the same re-credit SPPR_new applies, gated by fix_EE_0_cases.
            DC_pp = self.get_DC(DET_as_PP=True, normalize=False)
            flow_to_det = (self.M0 + self.egestion).fillna(0)
            det_fate = getattr(self, '_det_fate', None)
            PP_Import_seq = list(PP_seq + Import_seq)
            theta = self._resolve_det_param(det_theta, DET_seq, 1.0)
            non_DET_sppr = sppr_mat.drop(columns=list(DET_seq), errors='ignore').sum(axis=1)
            deadend_mask = (GE == 0).all(axis=1)
            deadend_seq = [g for g in deadend_mask.index[deadend_mask]
                           if g not in PP_Import_seq and g not in list(DET_seq)]
            apply_fix = fix_EE_0_cases and bool(deadend_seq) and len(DET_seq) == 1
            for i, det_l in enumerate(DET_seq):
                q_l = self.q[det_l] if self.q[det_l] > 0 else flow_to_det.sum()
                if det_fate is not None and det_l in det_fate.columns:
                    fr = det_fate[det_l].reindex(flow_to_det.index).fillna(0)
                    inflow = (flow_to_det[PP_Import_seq] * fr[PP_Import_seq]).sum()
                elif det_fate is not None and len(DET_seq) > 1:
                    inflow = 0.0
                else:
                    inflow = flow_to_det[PP_Import_seq].sum()
                if apply_fix:
                    a = sum(self.q[g] * (DC_pp.loc[g] * non_DET_sppr).sum() for g in deadend_seq)
                    b = sum(self.q[g] * (DC_pp.loc[g] * sppr_mat[det_l]).sum() for g in deadend_seq)
                    m = theta[i] * (inflow + a) / (q_l - theta[i] * b)
                else:
                    m = (inflow / q_l) * theta[i]
                sppr_mat[det_l] *= m
        elif det_open_mode == 'none' and det_collapse_mode == 'never':
            # DEFAULT: keep the proven exact symbolic per-DET scaling (byte-identical to old).
            for det_j in DET_seq:
                det_sym = sppr_vec.loc[det_j].values.ravel()[0]
                sppr_mat[det_j] *= float(sol_dict2[det_sym].evalf(subs=subs_dict))
        else:
            # Non-default: build the numeric (I-B)x=c from the symbolic basis sppr_mat and
            # route through the shared openness/stability/collapse solver.
            non_DET_sppr = sppr_mat.drop(columns=list(DET_seq), errors='ignore').sum(axis=1)
            B, c_vec = self._build_det_BC(sppr_mat, non_DET_sppr, DET_seq, DC, TE_option)
            sppr_mat = self._solve_det_scaling(
                B, c_vec, DET_seq, sppr_mat, non_DET_sppr, DC, TE_option,
                det_collapse_mode=det_collapse_mode, det_open_mode=det_open_mode,
                det_theta=det_theta, det_external_sppr=det_external_sppr)

        equations = equations.squeeze().tolist()
        variables = ordered_symbols

        return sppr_symbolic, sppr_mat, equations, variables
    
    def _SPPR_symbolic_helper_diet_import_as_DC(self, TE: Optional[pd.DataFrame], TE_option: str, DET_TE_vals: float, sppr_det_value: Optional[float],
                                                det_collapse_mode: str = 'never', det_open_mode: str = 'none',
                                                det_theta: float | dict = 1.0, det_external_sppr: float | dict = 0.0,
                                                fix_EE_0_cases: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, list, list]:
        """Symbolic SPPR helper, "diet import as DC" variant.

        Imported diet is kept as a separate production source whose own SPPR (DIET_SPPR_*) is
        solved from a second linear system, so each imported group carries the
        diet-composition-weighted production it requires. Solves the non-detritus SPPR block,
        the detritus block, then the diet-import block, and folds the DIET_* contributions back
        into the single Import column of the numeric basis matrix. Detritus columns use the exact
        symbolic scaling by default, or the shared numeric recycling solver when openness/collapse
        is requested.

        Args:
            TE (Optional[pd.DataFrame]): explicit TE matrix; if None it is built from TE_option.
            TE_option (str): transfer-efficiency / detritus-DC mode; one of 'GE', 'TE',
                'With Egestion', 'global'.
            DET_TE_vals (float): TE assigned to detritus rows when building the TE matrix.
            sppr_det_value (Optional[float]): if not None, every detritus column is scaled by
                this fixed value instead of being solved.
            det_collapse_mode (str): 'never', 'auto', or 'always' (see SPPR_new). Defaults to 'never'.
            det_open_mode (str): 'none', 'recycling_loss', or 'source_dilution' (see SPPR_new).
                Defaults to 'none'.
            det_theta (float | dict): detritus availability/retention; float or dict keyed by
                DET seq/name. Defaults to 1.0.
            det_external_sppr (float | dict): external SPPR for 'source_dilution'. Defaults to 0.0.
            fix_EE_0_cases (bool): accepted for signature parity with the as_PP helper (SPPR_symbolic
                forwards a shared kwargs dict to both). The EE=0 re-credit is the as_PP TE factor-0
                semantics; this variant does not use it. Defaults to True.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, list, list]: (sppr_symbolic, sppr_mat, equations,
            variables), where sppr_symbolic is the symbolic per-group solution, sppr_mat is the
            numeric basis matrix (groups x basal sources, with diet-import folded into the Import
            column), equations is the combined SPPR + diet-import symbolic system, and variables
            is the ordered list of DIET_SPPR and SPPR symbols.

        Raises:
            Exception: if TE_option is not one of the supported strings.
        """
        DET_seq = self.get_DET_seq()
        Regular_seq = self.get_Regular_seq()
        Import_seq = self.get_Import_seq()
        PP_seq = self.get_PP_seq()

        # get DC and TE matrix:
        DC = self.get_DC(DET_as_PP=False, normalize=False)
        GE = self.get_TE(TE_option=TE_option, DET_values=DET_TE_vals, as_matrix=True)
        flow2det = (self.M0 + self.egestion).sum()
        if TE_option == 'TE':
            DC.loc[DET_seq, Regular_seq] = 0
        elif TE_option == 'GE':
            det_fate_mat = getattr(self, '_det_fate', None)
            for det_j in DET_seq:
                q_j = self.q[det_j] if self.q[det_j] > 0 else flow2det / max(len(DET_seq), 1)
                if det_fate_mat is not None and det_j in det_fate_mat.columns and q_j > 0:
                    fracs = det_fate_mat[det_j].reindex(self.M0.index).fillna(0)
                    DC.loc[det_j, :] = ((self.M0 * fracs) / q_j).fillna(0)
                else:
                    DC.loc[det_j, :] = ((self.M0) / flow2det).fillna(0)
        elif TE_option == 'With Egestion':
            m = (self.M0 / flow2det).fillna(0)
            e = (self.egestion / flow2det).fillna(0)
            for det_j in DET_seq:
                DC.loc[det_j, :] = 0
            for det_j in DET_seq:
                DC.loc[det_j, :] = (m + e @ DC)
        else:
            raise Exception("diet_import_option should be in one of ['GE', 'TE', 'With Egestion', 'global']")

        if TE is not None:
            GE = TE.copy()

        # calculate A matrix and turn to symbolic matrix:
        A = (DC / GE).fillna(0).values.copy()
        A[GE.values == 0] = 0
        for i in range(len(DC)):  # Replace producer rows with identity rows
            if DC.values[i, :].sum() == 0:
                A[i, :] = 0
                A[i, i] = 1
        A = pd.DataFrame(A, index=DC.index, columns=DC.columns)

        # define and order symbols:
        index = Regular_seq + DET_seq + PP_seq
        diet_sppr_symbols = sm.symbols([f'DIET_SPPR_{self.seq2name[s]}'.replace(' ', '_') for s in index])
        diet_sppr_vec = pd.DataFrame(diet_sppr_symbols, index=index, columns=Import_seq)
        A = A.loc[index, :]

        # solve linear equation:
        symbols_by_trophic_info = {
            'Regular': sm.symbols([f'SPPR_{self.seq2name[s]}'.replace(' ', '_') for s in Regular_seq]),
            'DET': sm.symbols([f'SPPR_{self.seq2name[s]}'.replace(' ', '_') for s in DET_seq]),
            'PP': sm.symbols([f'SPPR_{self.seq2name[s]}'.replace(' ', '_') for s in PP_seq]),
        }
        ordered_symbols = symbols_by_trophic_info['Regular'] + symbols_by_trophic_info['DET'] \
                        + symbols_by_trophic_info['PP']
        sppr_vec = pd.DataFrame(ordered_symbols, index=index)
        equations = ((A.loc[index, index] @ sppr_vec).squeeze() + (A.loc[index, Import_seq] * diet_sppr_vec).squeeze()) - sppr_vec.squeeze()
        equation_DET = equations.loc[DET_seq]
        equations_non_DET = equations.loc[[i for i in index if i not in DET_seq]]
        sol = sm.linsolve(equations_non_DET.squeeze().tolist(), ordered_symbols)
        sol_tuple1 = list(sol)[0]
        sol_dict1 = dict(zip(ordered_symbols, sol_tuple1))
        sppr = sppr_vec.replace(sol_dict1)

        # solve DET equation:
        equation_DET = equation_DET.apply(lambda x: x.subs(sol_dict1) if hasattr(x, 'subs') else x)
        sol = sm.linsolve(list(np.ravel(equation_DET.values)), symbols_by_trophic_info['DET'])
        sol_tuple2 = list(sol)[0]
        sol_dict2 = dict(zip(symbols_by_trophic_info['DET'], sol_tuple2))
        sppr_det_symbol = sppr_vec.loc[DET_seq].values[0, 0]
        sppr_det = sol_dict2[sppr_det_symbol]

        # solve diet import equation:
        # Each imported group's own SPPR (DIET_SPPR_*) is the diet-composition-weighted SPPR of
        # what it eats, so it satisfies its own linear system once the within-system SPPR (with
        # PP pinned to 1 and detritus substituted) is known.
        sppr = sppr_vec.replace(sol_dict1)
        subs_dict_PP = {v: 1 for v in sppr_vec.squeeze().loc[PP_seq]}
        subs_dict = subs_dict_PP | {ds: sol_dict2[ds].subs(subs_dict_PP) for ds in symbols_by_trophic_info['DET']}
        sppr = sppr.apply(lambda col: col.map(lambda x: x.evalf(subs=subs_dict) if hasattr(x, 'evalf') else x))
        equations_di = (DC.loc[index, index] @ sppr).squeeze() + (DC.loc[index, Import_seq] * diet_sppr_vec).squeeze() - diet_sppr_vec.squeeze()

        equations_di = equations_di.squeeze().tolist()
        diet_symbols = diet_sppr_vec.squeeze().tolist()

        sol = sm.linsolve(equations_di, diet_symbols)
        sol_tuple3 = list(sol)[0]
        sol_dict3 = dict(zip(diet_sppr_symbols, sol_tuple3))

        # sub:
        target_free_seq = PP_seq
        sppr = sppr_vec.replace(sol_dict1)
        subs_dict = {v: 1 for v in sppr_vec.squeeze().loc[PP_seq]} | sol_dict3
        sppr_symbolic = sppr.apply(lambda col: col.map(lambda x: x.evalf(subs=subs_dict) if hasattr(x, 'evalf') else x))

        target_free_symbols = diet_sppr_symbols + symbols_by_trophic_info['DET'] + symbols_by_trophic_info['PP']
        target_free_seq = DET_seq + PP_seq
        sppr_mat, _ = sm.linear_eq_to_matrix(sol_tuple1, target_free_symbols)
        sppr_mat = sm.lambdify([], sppr_mat, 'numpy')() # Converts SymPy matrix to NumPy
        cols = [f'DIET_{s}'.replace(' ', '_') for s in index] + target_free_seq
        sppr_mat = pd.DataFrame(sppr_mat, index=index, columns=cols)
        if sppr_det_value is not None:
            sppr_mat[DET_seq] *= float(sppr_det_value)
        elif det_open_mode == 'none' and det_collapse_mode == 'never':
            # DEFAULT: keep the proven exact symbolic per-DET scaling (byte-identical to old).
            for det_j in DET_seq:
                det_sym = sppr_vec.loc[det_j].values.ravel()[0]
                sppr_mat[det_j] *= float(sol_dict2[det_sym].evalf(subs=subs_dict))
        else:
            # Non-default: build the numeric (I-B)x=c from the symbolic basis sppr_mat and
            # route through the shared openness/stability/collapse solver. DET columns only;
            # the DIET_* import columns are handled separately below.
            non_DET_sppr = sppr_mat.drop(
                columns=list(DET_seq) + [f'DIET_{s}'.replace(' ', '_') for s in index],
                errors='ignore').sum(axis=1)
            B, c_vec = self._build_det_BC(sppr_mat, non_DET_sppr, DET_seq, DC, TE_option)
            sppr_mat = self._solve_det_scaling(
                B, c_vec, DET_seq, sppr_mat, non_DET_sppr, DC, TE_option,
                det_collapse_mode=det_collapse_mode, det_open_mode=det_open_mode,
                det_theta=det_theta, det_external_sppr=det_external_sppr)
        for s in index:
            s_symbol = diet_sppr_vec.loc[s].values[0]
            sppr_mat[f'DIET_{s}'.replace(' ', '_')] *= float(sol_dict3[s_symbol])
        sppr_mat[Import_seq[0]] = sppr_mat.loc[:, [f'DIET_{s}'.replace(' ', '_') for s in index]].sum(axis=1)
        sppr_mat.drop(columns=[f'DIET_{s}'.replace(' ', '_') for s in index], inplace=True)

        equations = equations.squeeze().tolist() + equations_di
        variabls = diet_sppr_symbols + ordered_symbols

        return sppr_symbolic, sppr_mat, equations, variabls

    def SPPR_symbolic(self, TE: Optional[pd.DataFrame] = None, TE_option: str = 'GE', diet_import_option: str = 'as_DC', DET_TE_vals: float = 1,
                      sppr_det_value: Optional[float] = None,
                      det_collapse_mode: str = 'never', det_open_mode: str = 'none',
                      det_theta: float | dict = 1.0, det_external_sppr: float | dict = 0.0,
                      fix_EE_0_cases: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, list, list]:
        """Symbolic SPPR solver: dispatch to the selected diet-import helper.

        Applies the same detritus knobs as SPPR_new and forwards them to whichever diet-import
        helper is chosen. The default path (det_open_mode='none', det_collapse_mode='never')
        keeps the exact sympy per-DET scaling.

        Args:
            TE (Optional[pd.DataFrame]): explicit TE matrix; if None it is built from TE_option.
                Defaults to None.
            TE_option (str): transfer-efficiency mode; one of 'GE', 'TE', 'With Egestion',
                'global'. Defaults to 'GE'.
            diet_import_option (str): how imported diet is handled: 'as_DC' (imported diet kept
                as a separate production source with its own DIET_SPPR) or 'as_PP' (imported
                diet treated as an extra primary-production source). Defaults to 'as_DC'.
            DET_TE_vals (float): TE assigned to detritus rows. Defaults to 1.
            sppr_det_value (Optional[float]): if set, detritus columns are scaled by this fixed
                value instead of being solved. Defaults to None.
            det_collapse_mode (str): 'never', 'auto', or 'always' (see SPPR_new). Defaults to 'never'.
            det_open_mode (str): 'none', 'recycling_loss', or 'source_dilution' (see SPPR_new).
                Defaults to 'none'.
            det_theta (float | dict): detritus availability/retention; float or dict keyed by
                DET seq/name. Defaults to 1.0.
            det_external_sppr (float | dict): external SPPR for 'source_dilution'. Defaults to 0.0.
            fix_EE_0_cases (bool): mirrors SPPR_new's flag of the same name. When True (default),
                the TE / diet_import_option='as_PP' path re-credits the PP consumed by single-DET
                EE=0 dead-end groups back to the detritus pool, matching SPPR_new exactly. Only
                consumed by the as_PP TE factor-0 path (see _SPPR_symbolic_helper_diet_import_as_PP);
                a no-op for GE / With Egestion, for multi-DET, and for the as_DC variant.
                Defaults to True.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, list, list]: (sppr_symbolic, sppr_mat, equations,
            variables) from the selected helper (see _SPPR_symbolic_helper_diet_import_as_DC /
            _as_PP).
        """
        kwargs = dict(TE=TE, TE_option=TE_option, DET_TE_vals=DET_TE_vals,
                      sppr_det_value=sppr_det_value, det_collapse_mode=det_collapse_mode,
                      det_open_mode=det_open_mode, det_theta=det_theta,
                      det_external_sppr=det_external_sppr, fix_EE_0_cases=fix_EE_0_cases)
        if diet_import_option == 'as_PP':
            return self._SPPR_symbolic_helper_diet_import_as_PP(**kwargs)
        elif diet_import_option == 'as_DC':
            return self._SPPR_symbolic_helper_diet_import_as_DC(**kwargs)

    def _sample_SPPR_new_forced_balance(self, TE_option: str = 'TE', sppr_det: Optional[float] = None) -> tuple[pd.DataFrame, float]:
        """Run SPPR_new and force exact global balance by solving the single detritus SPPR.

        Treats the detritus scaling as one symbolic unknown x = sppr_det, then solves the scalar
        equation that forces the PP inflow to exactly equal the export outflow
        (catch + growth + net_migration weighted by SPPR), and scales the detritus columns by the
        solved value.

        Args:
            TE_option (str): transfer-efficiency mode; one of 'GE', 'TE', 'With Egestion',
                'global'. Defaults to 'GE'.
            sppr_det (Optional[float]): unused placeholder; the value is solved and returned.
                Defaults to None.

        Returns:
            tuple[pd.DataFrame, float]: (sppr, sppr_det), the balanced SPPR DataFrame and the
            single detritus scaling value that achieves balance.
        """
        sppr, _, _ = self.SPPR_new(
            TE_option=TE_option
        )

        x = sm.symbols('x')  # x = sppr_det
        sppr_det_vec = (sppr.loc[:, self.get_DET_seq()] * x).sum(axis=1)
        sppr_pp = sppr.loc[:, self.get_PP_seq()].sum(axis=1)
        sppr_symbolic = sppr_det_vec + sppr_pp
        
        inflow = self.p[self.get_PP_seq()].sum()
        outflow = ((self.catch + self.growth + self.net_migration)* sppr_symbolic).sum()

        sppr_det = float(sm.solve(inflow - outflow, x)[0])
        
        # multiply sppr[det] column by sppr_det
        sppr.loc[:, self.get_DET_seq()] *= sppr_det

        return sppr, sppr_det

    def monte_carlo_SPPR(self, n_samples: int = 1000, TE_error_percent: float = 10, TE_error_cut_percent: float = 20,
                            TE_option: str = 'GE', DET_TE_vals: float = 1, kind: str = 'new', diet_import_option: str = 'as_DC', silent: bool = True,
                            det_collapse_mode: str = 'never', det_open_mode: str = 'none',
                            det_theta: float | dict = 1.0, det_external_sppr: float | dict = 0.0) -> tuple:
        """Monte-Carlo uncertainty propagation over transfer efficiency.

        Repeatedly resamples the TE matrix from gamma distributions centred on the model TEs
        (relative error TE_error_percent, clipped at +/- TE_error_cut_percent), recomputes SPPR
        via SPPR_new or SPPR_symbolic, discards any sample that produces a negative SPPR (an
        unstable / non-physical draw), and averages the accepted samples.

        Args:
            n_samples (int, optional): number of SPPR samples to draw. Defaults to 1000.
            TE_error_percent (float, optional): TE standard deviation as a percentage of its mean
                (the gamma CV). Defaults to 10.
            TE_error_cut_percent (float, optional): clip band for each TE sample as a percentage
                of its mean. Defaults to 20.
            TE_option (str, optional): transfer-efficiency mode; one of 'GE', 'TE',
                'With Egestion', 'global'. Defaults to 'GE'.
            DET_TE_vals (float, optional): TE assigned to detritus rows. Defaults to 1.
            kind (str, optional): which solver to resample: 'new' (SPPR_new) or 'symbolic'
                (SPPR_symbolic). Defaults to 'new'.
            diet_import_option (str, optional): 'as_DC' or 'as_PP', passed to SPPR_symbolic when
                kind='symbolic'. Defaults to 'as_DC'.
            silent (bool, optional): suppress progress bars / prints. Defaults to True.
            det_collapse_mode (str, optional): 'never', 'auto', or 'always' (see SPPR_new);
                forwarded to every SPPR call. Defaults to 'never'.
            det_open_mode (str, optional): 'none', 'recycling_loss', or 'source_dilution' (see
                SPPR_new). Defaults to 'none'.
            det_theta (float | dict, optional): detritus availability/retention; float or dict
                keyed by DET seq/name. Defaults to 1.0.
            det_external_sppr (float | dict, optional): external SPPR for 'source_dilution'.
                Defaults to 0.0.

        Returns:
            tuple: (mean_sppr, accepted_samples_array, rejection_fraction, equations, variables),
            where mean_sppr is a pd.DataFrame averaged over accepted samples,
            accepted_samples_array is the np.ndarray of accepted SPPR samples,
            rejection_fraction is the float share of discarded samples, and equations/variables
            are the symbolic system from SPPR_symbolic (both None when kind='new').

        Raises:
            Exception: if kind is not 'new' or 'symbolic'.
        """

        # define basis sequence:
        PP_seq = self.get_PP_seq()
        Import_seq = self.get_Import_seq()
        DET_seq = self.get_DET_seq()
        if TE_option == 'GE':
            basis_seq = list(PP_seq) + list(DET_seq) + list(Import_seq)
        else:
            basis_seq = list(PP_seq) + list(Import_seq)
        basis_seq = sorted(basis_seq, reverse=True)

        # choose TE matrix:
        TE_means = self.get_TE(TE_option=TE_option, DET_values=DET_TE_vals, as_matrix=False)

        def sample_TE(TE_error_percent, TE_error_cut_percent):  # TE samplers as gamma distributions:
            # Draw one TE matrix: gamma with mean = model TE and CV = TE_error, clipped to the
            # +/- TE_error_cut_percent band, broadcast to a full matrix, with basal rows pinned to 1.
            # mean = shape * scale = TE
            # variance = shape * scale^2
            # std = sqrt(shape) * scale
            # std / mean = 1/sqrt(shape) = TE_error (given)
            # shape = 1/(TE_error^2)
            # scale = TE / shape = TE * TE_error^2
            TE_error = TE_error_percent / 100
            shape = 1/(TE_error**2)

            # Gamma keeps the sampled TE strictly positive (unlike a normal) with the requested
            # mean and CV; shape/scale chosen above so mean=TE_means and std/mean=TE_error.
            sampler = lambda: gamma.rvs(a=shape, scale=TE_means/shape)
            TE_sample = sampler()

            # Clip each TE to the +/- TE_error_cut_percent band to drop extreme tail draws.
            TE_high = (TE_means * (1 + TE_error_cut_percent/100)).values
            TE_low = (TE_means * (1 - TE_error_cut_percent/100)).values
            TE_sample[TE_sample >= TE_high] = TE_high[TE_sample >= TE_high]
            TE_sample[TE_sample <= TE_low] = TE_low[TE_sample <= TE_low]

            # Broadcast the per-group TE vector into a full n x n matrix, then pin basal rows to 1.
            TE_sample = pd.DataFrame([TE_sample]*self.n_groups, index=TE_means.index, columns=TE_means.index).T

            TE_sample.loc[basis_seq, :] = 1

            return TE_sample

        # Detritus openness/collapse knobs forwarded unchanged to every SPPR call below.
        det_kwargs = dict(det_collapse_mode=det_collapse_mode, det_open_mode=det_open_mode,
                          det_theta=det_theta, det_external_sppr=det_external_sppr)

        # initialize collectors:
        if kind == 'new':
            sppr, _, _ = self.SPPR_new(TE=None, TE_option=TE_option, DET_TE_vals=DET_TE_vals, **det_kwargs)
        elif kind == 'symbolic':
            _, sppr, e, v = self.SPPR_symbolic(TE=None, TE_option=TE_option, DET_TE_vals=DET_TE_vals, diet_import_option=diet_import_option, **det_kwargs)
        else:
            raise Exception(f'kind = {kind}')

        index = sppr.index
        columns = sppr.columns
        sppr_array = np.zeros((n_samples, len(index), len(columns)))
        not_counted_counter = 0
        counted_rows_array = np.ones(n_samples).astype(bool)

        # perform monte-carlo:
        for i in tqdm(range(n_samples), disable=silent, desc="monte-carlo on TE"):
            TE_sample = sample_TE(TE_error_percent, TE_error_cut_percent)
            if  kind == 'new':
                sppr, _, _ = self.SPPR_new(TE=TE_sample, TE_option=TE_option, DET_TE_vals=DET_TE_vals, **det_kwargs)
                # sppr = sppr.sort_index(ascending=False)
            elif kind == 'symbolic':
                _, sppr, _, _ = self.SPPR_symbolic(TE=TE_sample, TE_option=TE_option, DET_TE_vals=DET_TE_vals, diet_import_option=diet_import_option, **det_kwargs)
            # turn to numpy and collect:
            sppr = sppr.values
            # Rejection step: a negative SPPR means the resampled TE drove the detritus
            # recycling system unstable / non-physical, so drop this draw from the average.
            if np.any(sppr < -1e-10):
                not_counted_counter += 1
                counted_rows_array[i] = False
                continue
            sppr_array[i, :, :] = sppr

        # take average SPPR over the accepted (non-rejected) samples only:
        sppr = np.mean(sppr_array[counted_rows_array], axis=0)

        if not silent:
            print(f'    proportion of un-counted calculations: {not_counted_counter}/{n_samples}')

        # back to dataframe:
        sppr = pd.DataFrame(sppr, index=index, columns=columns)

        if kind == 'new':
            return sppr, sppr_array[counted_rows_array], not_counted_counter/n_samples, None, None
        else:
            return sppr, sppr_array[counted_rows_array], not_counted_counter/n_samples, e, v

    def monte_carlo_SPPR_2(self, n_samples: int = 1000, TE_error_percent: float = 10, TE_error_cut_percent: float = 20,
                         TE_option: str = 'GE', DET_TE_vals: float = 1, kind: str = 'new', silent: bool = True,
                         det_collapse_mode: str = 'never', det_open_mode: str = 'none',
                         det_theta: float | dict = 1.0, det_external_sppr: float | dict = 0.0) -> tuple[pd.DataFrame, np.ndarray, float]:
        """Variant of monte_carlo_SPPR supporting only kind='new'.

        Pre-allocates the sample array from the model's (n_groups x n_PP) shape rather than from
        the first SPPR call, but is otherwise the same gamma-resampling / negative-rejection /
        averaging loop as monte_carlo_SPPR.

        Args:
            n_samples (int, optional): number of SPPR samples to draw. Defaults to 1000.
            TE_error_percent (float, optional): TE standard deviation as a percentage of its mean
                (the gamma CV). Defaults to 10.
            TE_error_cut_percent (float, optional): clip band for each TE sample as a percentage
                of its mean. Defaults to 20.
            TE_option (str, optional): transfer-efficiency mode; one of 'GE', 'TE',
                'With Egestion', 'global'. Defaults to 'GE'.
            DET_TE_vals (float, optional): TE assigned to detritus rows. Defaults to 1.
            kind (str, optional): must be 'new'. Defaults to 'new'.
            silent (bool, optional): suppress progress bars / prints. Defaults to True.
            det_collapse_mode (str, optional): 'never', 'auto', or 'always' (see SPPR_new).
                Defaults to 'never'.
            det_open_mode (str, optional): 'none', 'recycling_loss', or 'source_dilution' (see
                SPPR_new). Defaults to 'none'.
            det_theta (float | dict, optional): detritus availability/retention; float or dict
                keyed by DET seq/name. Defaults to 1.0.
            det_external_sppr (float | dict, optional): external SPPR for 'source_dilution'.
                Defaults to 0.0.

        Returns:
            tuple[pd.DataFrame, np.ndarray, float]: (mean_sppr, accepted_samples_array,
            rejection_fraction), the per-group SPPR averaged over accepted samples, the array of
            accepted SPPR samples, and the share of discarded samples.

        Raises:
            Exception: if kind is not 'new'.
        """

        # define basis sequence:
        PP_seq = self.get_PP_seq()
        Import_seq = self.get_Import_seq()
        DET_seq = self.get_DET_seq()
        if TE_option == 'GE':
            basis_seq = list(PP_seq) + list(DET_seq) + list(Import_seq)
        else:
            basis_seq = list(PP_seq) + list(Import_seq)
        basis_seq = sorted(basis_seq, reverse=True)
            
        # choose TE matrix:
        TE_means = self.get_TE(TE_option=TE_option, DET_values=DET_TE_vals, as_matrix=False)

        def sample_TE(TE_error_percent, TE_error_cut_percent):  # TE samplers as gamma distributions:
            # mean = shape * scale = TE
            # variance = shape * scale^2
            # std = sqrt(shape) * scale
            # std / mean = 1/sqrt(shape) = TE_error (given)
            # shape = 1/(TE_error^2)
            # scale = TE / shape = TE * TE_error^2
            TE_error = TE_error_percent / 100
            shape = 1/(TE_error**2)

            # Gamma keeps the sampled TE strictly positive with the requested mean and CV.
            sampler = lambda: gamma.rvs(a=shape, scale=TE_means/shape)
            TE_sample = sampler()

            # Clip each TE to the +/- TE_error_cut_percent band to drop extreme tail draws.
            TE_high = (TE_means * (1 + TE_error_cut_percent/100)).values
            TE_low = (TE_means * (1 - TE_error_cut_percent/100)).values
            TE_sample[TE_sample >= TE_high] = TE_high[TE_sample >= TE_high]
            TE_sample[TE_sample <= TE_low] = TE_low[TE_sample <= TE_low]

            # Broadcast the per-group TE vector into a full matrix, then pin basal rows to 1.
            TE_sample = pd.DataFrame([TE_sample]*n_groups, index=TE_means.index, columns=TE_means.index).T

            TE_sample.loc[basis_seq, :] = 1

            return TE_sample

        # Detritus openness/collapse knobs forwarded unchanged to every SPPR call below.
        det_kwargs = dict(det_collapse_mode=det_collapse_mode, det_open_mode=det_open_mode,
                          det_theta=det_theta, det_external_sppr=det_external_sppr)

        # initialize collectors:
        sppr, _, _ = self.SPPR_new(TE=None, TE_option=TE_option, DET_TE_vals=DET_TE_vals, **det_kwargs)
        index = sppr.index
        columns = sppr.columns
        n_PP = len(columns)
        n_groups = self.n_groups
        sppr_array = np.zeros((n_samples, n_groups, n_PP))
        not_counted_counter = 0
        counted_rows_array = np.ones(n_samples).astype(bool)

        # perform monte-carlo:
        for i in tqdm(range(n_samples), disable=silent):
            TE_sample = sample_TE(TE_error_percent, TE_error_cut_percent)
            if  kind == 'new':
                sppr, _, _ = self.SPPR_new(TE=TE_sample, TE_option=TE_option, DET_TE_vals=DET_TE_vals, **det_kwargs)
                # sppr = sppr.sort_index(ascending=False)
            else:
                raise Exception(f'kind = {kind}')
            # turn to numpy and collect:
            sppr = sppr.values
            # Reject negative (unstable / non-physical) draws so they don't bias the mean.
            if np.any(sppr < -1e-10):
                not_counted_counter += 1
                counted_rows_array[i] = False
                continue
            sppr_array[i, :, :] = sppr

        # take average SPPR over the accepted (non-rejected) samples only:
        sppr = np.mean(sppr_array[counted_rows_array], axis=0)

        if not silent:
            print(f'    proportion of un-counted calculations: {not_counted_counter}/{n_samples}')

        # back to dataframe:
        sppr = pd.DataFrame(sppr, index=index, columns=columns)
        return sppr, sppr_array[counted_rows_array], not_counted_counter/n_samples
    
    # class methods:
    @classmethod
    def rename_results(cls, results: list | pd.DataFrame | pd.Series, renaming_dict: dict) -> list | pd.DataFrame | pd.Series:
        """Relabel the index (and columns) of one or more SPPR-style results.

        Typically used to map between group seq IDs and group names. Each result is first sorted
        by descending index (and columns for DataFrames), then renamed.

        Args:
            results (list | pd.DataFrame | pd.Series): a single result or a list of results to
                relabel.
            renaming_dict (dict): mapping applied to index labels (and column labels for
                DataFrames), e.g. seq2name or name2seq.

        Returns:
            list | pd.DataFrame | pd.Series: the relabeled result(s), returned as the same type
            (single object or list) that was passed in.
        """
        is_list = isinstance(results, list)
        results = results if isinstance(results, list) else [results]
        for i in range(len(results)):
            r = results[i]
            if isinstance(r, pd.DataFrame):
                r = r.sort_index(ascending=False).sort_index(axis=1, ascending=False).rename(index=renaming_dict, columns=renaming_dict)
            elif isinstance(r, pd.Series):
                r = r.sort_index(ascending=False).rename(index=renaming_dict)
            results[i] = r
        if not is_list:
            return results[0]
        return results
