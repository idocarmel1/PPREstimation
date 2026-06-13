import numpy as np
import pandas as pd
import sympy as sm
from tqdm.notebook import tqdm
from scipy.stats import gamma
import igraph as ig
from scipy.optimize import minimize

from ModelData import ModelData
from utils import mat_from_np, remove_cycles, move_scattered_identity
from copy import deepcopy

class PPRCalculator:

    # constructors:
    def __init__(self, model_number, underdetermined=False, zero_catch=True, zero_biomass_accum=True, default_gs=True, weight_flow=1.0, weight_guess=1.0):
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
    def from_dict(cls, data_dict, underdetermined=False, zero_catch=True, zero_biomass_accum=True, default_gs=True, weight_flow=1.0, weight_guess=1.0):
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
    def from_modeldata(cls, modeldata: ModelData, underdetermined=False, zero_catch=True, zero_biomass_accum=True, default_gs=True, weight_flow=1.0, weight_guess=1.0):
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

    def _sort(self):
        for name, value in vars(self).items():
            if isinstance(value, pd.Series):
                setattr(self, name, value.sort_index(ascending=False))
            elif isinstance(value, pd.DataFrame) and name != '_groups_df':
                setattr(self, name, value.sort_index(ascending=False).sort_index(ascending=False, axis=1))
        return self

    def _fill_properties(self, groups_df):

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
        self.EE = groups_df['ee'].fillna(1).copy()
        self.GE = groups_df['ge'].fillna(1).copy()
        self.GS = groups_df['gs'].fillna(0).copy()
        self.TL = groups_df['tl'].fillna(1).copy()

        # balance:
        self.n_balance_runs = 0
        self.is_balanced, _, _ = self.is_model_balanced()
        self.balanced_model = self.balance_model(change_production=False)
    
    @classmethod
    def apply_ecopath_defaults(cls, df, DC, det_fate=None, zero_catch=False, zero_biomass_accum=False, default_gs=False):
        """
        Applies Ecopath defaults and ensures flows are synced with ratios.
        Assumes p, q, ee, and catch are given. assumes p, q not given for detritus group.
            calculates M0 = p*(1-ee)
        Assumes gs is 0.2 for regular groups and 0 for pp, det, diet_import.
            calculates egestion = q*gs.
        Assumes biomas_accum, net_migration are 0 where currently nan.
        Calculates once cell missing from p, q, egestion, respiration.
        Calculates predation = Z.sum(axis=1).
            Calculates once cell missing from M0, predation, catch, biomass_accum, net_migration.
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
        if det_fate is not None and len(det_idx) > 1:
            for det_j in det_idx:
                if det_j in det_fate.columns:
                    df.loc[det_j, 'q'] = (df['flow_to_det'] * det_fate[det_j].reindex(df.index).fillna(0)).sum()
                else:
                    df.loc[det_j, 'q'] = df['flow_to_det'].sum()
        else:
            df.loc[is_det, 'q'] = df['flow_to_det'].sum()
        df.loc[is_det, 'p'] = df.loc[is_det, 'q']
        df.loc[is_det, 'biomass_accum'] = df.loc[is_det, 'p'] - (df.loc[is_det, 'predation'] + df.loc[is_det, 'net_migration'])

        return df

    def apply_lim(self, df, weight_flow=1.0, weight_guess=0.0):
        """
        Uses Linear Inverse Modeling (SLSQP) to fill in missing mass-balance 
        variables by dynamically registering active constraints and coordinating guesses.
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

            # rebalance DET rows:
            det_seqs = self.get_DET_seq()
            det_fate = getattr(self, '_det_fate', None)
            if det_fate is not None and len(det_seqs) > 1:
                for det_j in det_seqs:
                    if det_j in det_fate.columns:
                        df_final.loc[det_j, 'q'] = (df_final['flow_to_det'] * det_fate[det_j].reindex(df_final.index).fillna(0)).sum()
                    else:
                        df_final.loc[det_j, 'q'] = df_final['flow_to_det'].sum()
            else:
                df_final.loc[det_seqs, 'q'] = df_final['flow_to_det'].sum()
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
    def is_model_balanced(self):

        # production = catch + predation + growth + net_migration + M0
        production = self.catch + self.predation + self.growth + self.net_migration + self.M0
        p_is_balanced = all(np.isclose(production, self.p))

        # consumpotion = production + egestion + respiration
        consumption = production + self.egestion + self.respiration
        q_is_balanced =  all(np.isclose(consumption, self.q, atol=1e-6))

        model_is_balanced = p_is_balanced and q_is_balanced

        return model_is_balanced, production, consumption
    
    def balance_model(self, change_production=False):
        """rebalance by changing growth, and net_migration (if change_production=False) or production (otherwise)"""
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
    
    def is_sppr_balanced(self, sppr, diet_import_equations=None):
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
    def get_model(self):
        if hasattr(self, "_model"):
            return self._model
        else:
            return None
    
    def get_groups_df(self):
        return self._groups_df.copy().sort_index(ascending=False)
    
    def get_DC(self, DET_as_PP=True, normalize=False):
        """
        Args:
            DET_as_PP (bool, optional): if True, DET row is set to 1. otherwise, it is (M0 + egestion)/(flow2det).
                Defaults to True.
            normalize (bool, optional): if True, DC rows sum to 1. else, they keep original sum.
                Defaults to normalize.
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
    
    def get_Z(self, DET_as_PP=False):
        """get Z matrix. if DET_as_PP is False (default), DET rows are flow_to_det split by det_fate. otherwise set to 0"""
        Z = self._DC.mul(self._groups_df['q'].fillna(0), axis='index')
        DET_seq = self.get_DET_seq()
        if not DET_as_PP:
            flow_to_det = (self.M0 + self.egestion).fillna(0)
            det_fate = getattr(self, '_det_fate', None)
            if det_fate is not None and len(DET_seq) > 1:
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
    
    def get_DET_seq(self):
        return sorted(self._groups_df.index[self._groups_df['trophic_info'] == 'DET'].values)
    
    def get_PP_seq(self):
        return sorted(self._groups_df.index[self._groups_df['trophic_info'] == 'PP'].values)
    
    def get_Regular_seq(self):
        return sorted(self._groups_df.index[self._groups_df['trophic_info'] == 'Regular'].values)
    
    def get_Import_seq(self):
        return sorted(self._groups_df.index[self._groups_df['trophic_info'] == 'Import'].values)
    
    def get_TE(self, TE_option: str, DET_values=1, as_matrix=True, global_TE='mean'):
        """
        Args:
            TE_option (str): should be one of ['GE', 'TE', 'With Egestion', 'global']. if 'global', global_value must be set.
            DET_values (int, optional): TE of Detritus (from detritus to others). Defaults to 1.
            as_matrix (bool, optional): whether to return pd.Dataframe (nxn) or pd.Series (nx1). Defaults to True.
            global_TE (str or int, optional): global value for TE in case of TE_option == 'global'. 
                Defaults to 'mean', and then it is the mean TE of catch if sum(catch)!=0, else to mean of biomass.
        Returns:
            if as_matrix: TE DataFrame, else TE Series.
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
    
    def get_TL(self, break_cycles: bool, DET_as_PP: bool, TE_option='With Egestion'):

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

    def get_PPR(self, sppr, only_inner=False):
        sppr = sppr.copy().fillna(0) # sppr is an output of an SPPR calculating method from this class.
        sppr = PPRCalculator.rename_results(sppr, self.name2seq)
        sppr = sppr.reindex(self.catch.index, fill_value=0)
        sppr = sppr.replace(np.inf, 0)

        if only_inner and (set(self.get_Import_seq()).issubset(set(sppr.columns))):
            sppr = sppr.drop(columns=self.get_Import_seq())

        if isinstance(sppr, pd.DataFrame):
            return self.catch.dot(sppr).to_frame().T
        else:
            return self.catch.dot(sppr)
    
    def get_NPP(self, only_inner=True):
        if only_inner:
            return self.p[self.get_PP_seq()].sum()
        else:
            raise Exception('not implemented yet')
    
    def get_PPR2NPP_ratio(self, sppr):
        return self.get_PPR(sppr, only_inner=True).sum(axis=1).sum() / self.get_NPP(only_inner=True)
        
    ##########################################################################################
    ############################## SPPR calculating methods ##################################
    ##########################################################################################
    def SPPR_1986(self):
        if all(self.catch == 0):
            return pd.DataFrame(0, index=self.GE.index, columns=['sppr'])
        TE = 0.1
        TL = np.average(self.get_TL(break_cycles=True, DET_as_PP=True), weights=self.catch)  # break cylces in DC, force DET to be of TL=1 before TL calculation.
        SPPR = TE ** (1 - TL)
        SPPR = pd.DataFrame(SPPR, index=self.GE.index, columns=['sppr'])
        return SPPR
    
    def SPPR_1995(self, global_TE=0.1):
        """
        Args:
            global_TE (float, optional): 'mean' or float. Defaults to 0.1.
        """
        TE = self.get_TE(TE_option='global', global_TE=global_TE, as_matrix=False)
        TL = self.get_TL(break_cycles=True, DET_as_PP=True)  # break cylces in DC and force DET to be of TL=1.
        SPPR = TE ** (1 - TL)
        return SPPR.to_frame(name='sppr')
    
    def SPPR_1995_TL_fix(self, global_TE=0.1):
        TE = self.get_TE(TE_option='global', global_TE=global_TE, as_matrix=False)
        TL = self.get_TL(break_cycles=True, DET_as_PP=True)  # break cylces in DC and force DET to be of TL=1.
        TL_fraction = TL % 1
        TL_int = TL.astype(int)
        sppr = (1-TL_fraction) * (1/TE)**(TL_int-1) + TL_fraction * (1/TE)**(TL_int)
        sppr = sppr.fillna(1)
        return sppr.to_frame(name='sppr')
        
    def SPPR_EwE(self, TE_option: str, use_EE=True, return_paths=True, silent=True):
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

                    prod = 1.0
                    for i in range(len(path) - 1):
                        prod *= A_list[path[i]][path[i + 1]]

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

            # Calculate products using reduceat (C-level segmented product)
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

            # Vectorized 'Scatter-Add'
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

    def SPPR_EwE_Ido(self, TE_option: str, global_TE='mean', use_EE=True):
        """
        Args:
            TE_option (str): should be one of ['GE', 'TE', 'With Egestion', 'global']. if 'global', global_value must be set.
            global_TE (str or int, optional): global value for TE in case of TE_option == 'global'. 
                Defaults to 'mean', and then it is the mean TE of catch if sum(catch)!=0, else to mean of biomass.
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
        
        # normalize by the first primary producer:
        ns = [np.array(s).T.flatten() for s in ns]
        M = sm.Matrix(ns)
        rref_matrix, _ = M.rref()
        basis = [np.array(rref_matrix.row(i).evalf()).astype(float).flatten()
                    for i in range(rref_matrix.rows)
                    if not rref_matrix.row(i).is_zero]
        basis = np.array(basis, dtype=float).T

        # turn back to DataFrames:
        SPPR = pd.DataFrame(basis, index=new_index)
        SPPR = SPPR.rename(columns=lambda c: SPPR.index.values[SPPR[c] == 1][0])
        if use_EE:
            SPPR = SPPR.mul(self.EE, axis='index')
        A = pd.DataFrame(np.array(A.tolist(), dtype=float), index=new_index, columns=new_columns)
        L = pd.DataFrame(np.array(L.tolist(), dtype=float), index=new_index, columns=new_columns)

        return SPPR, A, L

    def SPPR_2015(self, only_pp_det=True):
        only_pp_det=True

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
            if only_pp_det:  # this is what is implemented in the article
                for i in PP_seq:
                    Z_without_DET.loc[:, i] += percent_of_det_that_is_PP[i] * Z_without_DET.loc[:, det_j]
            else:
                Z_without_DET.loc[:, PP_seq] += Z_without_DET.loc[:, det_j]
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
        L.loc[seq_to_drop, seq_to_drop] = 1
        SPPR = L.loc[PP_seq, :].T

        # add back sppr_det that makes model balanced:
        SPPR = SPPR.reindex(self.catch.index, fill_value=0)
        det_seqs = self.get_DET_seq()
        if det_seqs:
            for i in PP_seq:
                SPPR.loc[det_seqs, i] = self.M0.loc[PP_seq][i] / (self.M0 + self.egestion).sum()

        return SPPR, A, L
        
    def SPPR_new(self, TE=None, TE_option='GE', DET_TE_vals=1):
        # get DC:
        DC = self.get_DC(DET_as_PP=True, normalize=False)
        
        # get TE matrix:
        if TE is not None:
            GE = TE.copy()
        else:
            GE = self.get_TE(TE_option=TE_option, DET_values=DET_TE_vals, as_matrix=True)

        # calculate A matrix and turn to symbolic matrix:
        A = (DC / GE).fillna(0).values
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
        
        # normalize by the first primary producer:
        ns = [np.array(s).T.flatten() for s in ns]
        M = sm.Matrix(ns)
        rref_matrix, _ = M.rref()
        basis = [np.array(rref_matrix.row(i).evalf()).astype(float).flatten()
                    for i in range(rref_matrix.rows)
                    if not rref_matrix.row(i).is_zero]
        basis = np.array(basis, dtype=float).T

        # turn back to DataFrames:
        SPPR = pd.DataFrame(basis, index=new_index)
        SPPR = SPPR.rename(columns=lambda c: SPPR.index.values[SPPR[c] == 1][0])
        A = pd.DataFrame(np.array(A.tolist(), dtype=float), index=new_index, columns=new_columns)
        L = pd.DataFrame(np.array(L.tolist(), dtype=float), index=new_index, columns=new_columns)

        # find sppr_det if it is in the output:
        DET_seq = self.get_DET_seq()
        det_fate = getattr(self, '_det_fate', None)
        flow_to_det = (self.M0 + self.egestion).fillna(0)
        non_DET_sppr = SPPR.drop(columns=DET_seq).sum(axis=1)

        if len(DET_seq) == 1:
            det_j = DET_seq[0]
            flow2det = self.q[det_j] if self.q[det_j] > 0 else flow_to_det.sum()
            if TE_option == 'TE':
                PP_Import_seq = list(self.get_PP_seq() + self.get_Import_seq())
                sppr_det = flow_to_det[PP_Import_seq].sum() / flow2det
                SPPR[DET_seq] *= sppr_det
            elif TE_option == 'GE':
                m = (self.M0 / flow2det).fillna(0)
                x = sm.symbols('x')
                sppr = non_DET_sppr + x * SPPR[det_j]
                sppr_det = float(sm.solve(x - (m @ sppr), x)[0])
                SPPR[DET_seq] *= sppr_det
            elif TE_option == 'With Egestion':
                m = (self.M0 / flow2det).fillna(0)
                e = (self.egestion / flow2det).fillna(0)
                x = sm.symbols('x')
                sppr = non_DET_sppr + x * SPPR[det_j]
                sppr_det = float(sm.solve(x - (m @ sppr + (DC @ sppr) @ e), x)[0])
                SPPR[DET_seq] *= sppr_det
            else:
                raise Exception("TE_option should be in ['GE', 'TE', 'With Egestion', 'global']")
        else:
            # Multiple DET groups: solve coupled linear system (I - B) x = c
            # x_l = sppr_det for DET group l
            # B[l,j] = m_l @ SPPR[:,DET_j],  c[l] = m_l @ non_DET_sppr
            k = len(DET_seq)
            PP_Import_seq = list(self.get_PP_seq() + self.get_Import_seq())

            def _get_m_l(det_l):
                q_l = self.q[det_l] if self.q[det_l] > 0 else 1.0
                if det_fate is not None and det_l in det_fate.columns:
                    fracs = det_fate[det_l].reindex(self.M0.index).fillna(0)
                    return (self.M0 * fracs / q_l).fillna(0)
                return (self.M0 / q_l).fillna(0)

            if TE_option == 'TE':
                for det_l in DET_seq:
                    q_l = self.q[det_l] if self.q[det_l] > 0 else 1.0
                    if det_fate is not None and det_l in det_fate.columns:
                        fracs = det_fate[det_l].reindex(flow_to_det.index).fillna(0)
                        inflow_from_PP_Import = (flow_to_det[PP_Import_seq] * fracs[PP_Import_seq]).sum()
                    else:
                        inflow_from_PP_Import = flow_to_det[PP_Import_seq].sum()
                    SPPR[det_l] *= inflow_from_PP_Import / q_l
            elif TE_option == 'GE':
                B = np.zeros((k, k))
                c_vec = np.zeros(k)
                for li, det_l in enumerate(DET_seq):
                    m_l = _get_m_l(det_l)
                    c_vec[li] = float(m_l @ non_DET_sppr)
                    for ji, det_j in enumerate(DET_seq):
                        B[li, ji] = float(m_l @ SPPR[det_j])
                x_vec = np.linalg.solve(np.eye(k) - B, c_vec)
                for i, det_j in enumerate(DET_seq):
                    SPPR[det_j] *= x_vec[i]
            elif TE_option == 'With Egestion':
                B = np.zeros((k, k))
                c_vec = np.zeros(k)
                for li, det_l in enumerate(DET_seq):
                    q_l = self.q[det_l] if self.q[det_l] > 0 else 1.0
                    if det_fate is not None and det_l in det_fate.columns:
                        fracs = det_fate[det_l].reindex(self.M0.index).fillna(0)
                        m_l = (self.M0 * fracs / q_l).fillna(0)
                        e_l = (self.egestion * fracs / q_l).fillna(0)
                    else:
                        m_l = (self.M0 / q_l).fillna(0)
                        e_l = (self.egestion / q_l).fillna(0)
                    m_eff_l = m_l + DC.T @ e_l
                    c_vec[li] = float(m_eff_l @ non_DET_sppr)
                    for ji, det_j in enumerate(DET_seq):
                        B[li, ji] = float(m_eff_l @ SPPR[det_j])
                x_vec = np.linalg.solve(np.eye(k) - B, c_vec)
                for i, det_j in enumerate(DET_seq):
                    SPPR[det_j] *= x_vec[i]
            else:
                raise Exception("TE_option should be in ['GE', 'TE', 'With Egestion', 'global']")

        return SPPR, A, L
    
    def _SPPR_symbolic_helper_diet_import_as_PP(self, TE, TE_option, DET_TE_vals, sppr_det_value):
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
        if sppr_det_value is None:
            for det_j in DET_seq:
                det_sym = sppr_vec.loc[det_j].values.ravel()[0]
                sppr_mat[det_j] *= float(sol_dict2[det_sym].evalf(subs=subs_dict))
        else:
            sppr_mat[DET_seq] *= float(sppr_det_value)

        equations = equations.squeeze().tolist()
        variables = ordered_symbols

        return sppr_symbolic, sppr_mat, equations, variables
    
    def _SPPR_symbolic_helper_diet_import_as_DC(self, TE, TE_option, DET_TE_vals, sppr_det_value):
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
        if sppr_det_value is None:
            for det_j in DET_seq:
                det_sym = sppr_vec.loc[det_j].values.ravel()[0]
                sppr_mat[det_j] *= float(sol_dict2[det_sym].evalf(subs=subs_dict))
        else:
            sppr_mat[DET_seq] *= float(sppr_det_value)
        for s in index:
            s_symbol = diet_sppr_vec.loc[s].values[0]
            sppr_mat[f'DIET_{s}'.replace(' ', '_')] *= float(sol_dict3[s_symbol])
        sppr_mat[Import_seq[0]] = sppr_mat.loc[:, [f'DIET_{s}'.replace(' ', '_') for s in index]].sum(axis=1)
        sppr_mat.drop(columns=[f'DIET_{s}'.replace(' ', '_') for s in index], inplace=True)

        equations = equations.squeeze().tolist() + equations_di
        variabls = diet_sppr_symbols + ordered_symbols

        return sppr_symbolic, sppr_mat, equations, variabls

    def SPPR_symbolic(self, TE=None, TE_option='GE', diet_import_option='as_DC', DET_TE_vals=1, sppr_det_value=None):
        if diet_import_option == 'as_PP':
            return self._SPPR_symbolic_helper_diet_import_as_PP(TE=TE, TE_option=TE_option, DET_TE_vals=DET_TE_vals, sppr_det_value=sppr_det_value)
        elif diet_import_option == 'as_DC':
            return self._SPPR_symbolic_helper_diet_import_as_DC(TE=TE, TE_option=TE_option, DET_TE_vals=DET_TE_vals, sppr_det_value=sppr_det_value)

    def _sample_SPPR_new_forced_balance(self, TE=None, sppr_det=None):
        sppr, _, _ = self.SPPR_new(
            DET_modeling='as_PP', DET_TE_vals=1, TE=TE
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

    def monte_carlo_SPPR(self, n_samples=1000, TE_error_percent=10, TE_error_cut_percent=20,
                            TE_option='GE', DET_TE_vals=1, kind='new', diet_import_option='as_DC', silent=True):
        """

        Args:
            n_samples (int, optional): number of sppr samples. Defaults to 1000.
            TE_error_percent (float, optional): percentage of TE std relative to it's mean. Defaults to 0.1.
            TE_error_cut_percent (float, optional): cut value to TE in percentage relative to it's mean. Defaults to 0.2.
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

            sampler = lambda: gamma.rvs(a=shape, scale=TE_means/shape)
            TE_sample = sampler()

            TE_high = (TE_means * (1 + TE_error_cut_percent/100)).values
            TE_low = (TE_means * (1 - TE_error_cut_percent/100)).values
            TE_sample[TE_sample >= TE_high] = TE_high[TE_sample >= TE_high]
            TE_sample[TE_sample <= TE_low] = TE_low[TE_sample <= TE_low]

            TE_sample = pd.DataFrame([TE_sample]*self.n_groups, index=TE_means.index, columns=TE_means.index).T

            TE_sample.loc[basis_seq, :] = 1

            return TE_sample

        # initialize collectors:
        if kind == 'new':
            sppr, _, _ = self.SPPR_new(TE=None, TE_option=TE_option, DET_TE_vals=DET_TE_vals)
        elif kind == 'symbolic':
            _, sppr, e, v = self.SPPR_symbolic(TE=None, TE_option=TE_option, DET_TE_vals=DET_TE_vals, diet_import_option=diet_import_option)
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
                sppr, _, _ = self.SPPR_new(TE=TE_sample, TE_option=TE_option, DET_TE_vals=DET_TE_vals)
                # sppr = sppr.sort_index(ascending=False)
            elif kind == 'symbolic':
                _, sppr, _, _ = self.SPPR_symbolic(TE=TE_sample, TE_option=TE_option, DET_TE_vals=DET_TE_vals, diet_import_option=diet_import_option)
            # turn to numpy and collect:
            sppr = sppr.values
            if np.any(sppr < -1e-10):
                not_counted_counter += 1
                counted_rows_array[i] = False
                continue
            sppr_array[i, :, :] = sppr
                
        # take average SPPR:
        sppr = np.mean(sppr_array[counted_rows_array], axis=0)

        if not silent:
            print(f'    proportion of un-counted calculations: {not_counted_counter}/{n_samples}')

        # back to dataframe:
        sppr = pd.DataFrame(sppr, index=index, columns=columns)

        if kind == 'new':
            return sppr, sppr_array[counted_rows_array], not_counted_counter/n_samples, None, None
        else:
            return sppr, sppr_array[counted_rows_array], not_counted_counter/n_samples, e, v

    def monte_carlo_SPPR_2(self, n_samples=1000, TE_error_percent=10, TE_error_cut_percent=20,
                         TE_option='GE', DET_TE_vals=1, kind='new', silent=True):
        """

        Args:
            n_samples (int, optional): number of sppr samples. Defaults to 1000.
            TE_error_percent (float, optional): percentage of TE std relative to it's mean. Defaults to 0.1.
            TE_error_cut_percent (float, optional): cut value to TE in percentage relative to it's mean. Defaults to 0.2.
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

            sampler = lambda: gamma.rvs(a=shape, scale=TE_means/shape)
            TE_sample = sampler()

            TE_high = (TE_means * (1 + TE_error_cut_percent/100)).values
            TE_low = (TE_means * (1 - TE_error_cut_percent/100)).values
            TE_sample[TE_sample >= TE_high] = TE_high[TE_sample >= TE_high]
            TE_sample[TE_sample <= TE_low] = TE_low[TE_sample <= TE_low]

            TE_sample = pd.DataFrame([TE_sample]*n_groups, index=TE_means.index, columns=TE_means.index).T

            TE_sample.loc[basis_seq, :] = 1

            return TE_sample

        # initialize collectors:
        sppr, _, _ = self.SPPR_new(TE=None, TE_option=TE_option, DET_TE_vals=DET_TE_vals)
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
                sppr, _, _ = self.SPPR_new(TE=TE_sample, TE_option=TE_option, DET_TE_vals=DET_TE_vals)
                # sppr = sppr.sort_index(ascending=False)
            else:
                raise Exception(f'kind = {kind}')
            # turn to numpy and collect:
            sppr = sppr.values
            if np.any(sppr < -1e-10):
                not_counted_counter += 1
                counted_rows_array[i] = False
                continue
            sppr_array[i, :, :] = sppr
            
        # take average SPPR:
        sppr = np.mean(sppr_array[counted_rows_array], axis=0)

        if not silent:
            print(f'    proportion of un-counted calculations: {not_counted_counter}/{n_samples}')

        # back to dataframe:
        sppr = pd.DataFrame(sppr, index=index, columns=columns)
        return sppr, sppr_array[counted_rows_array], not_counted_counter/n_samples
    
    # class methods:
    @classmethod
    def rename_results(cls, results: list, renaming_dict):
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
