import os
import sys

import pandas as pd
import numpy as np
from tqdm import tqdm

from PPRCalculator import PPRCalculator
from utils import species_groups

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


df = pd.DataFrame(species_groups)
model_numbers = sorted(list(df['model_number'].unique()))
collect_mc_sppr = True
n_samples = 100


def main():
    path = 'real_models/'
    saved_models = [int(f.split('[')[1].split(']')[0]) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

    # saved_models = []
    # model_numbers = [227, 462]

    for model_number in tqdm(model_numbers):
        print(model_number)
        # print(model_number)
        if model_number in saved_models:
            continue
        # create model calculator:
        try:
            model = PPRCalculator(model_number)

            # create name for the model
            year = model.get_model().model_year
            name = model.get_model().model_name
            model_name = f'[{model_number}] {name} {year}'
            # print(model_name)
            # don't treat models with more than one DET row yet:
            if len(model.get_DET_seq()) > 1:
                # print(f'    model {model_name} failed - has more than 1 DET row')
                # print('-'*20)
                continue
            # don't treat models where catch is 0:
            if model.catch.sum() == 0:
                # print(f'    model {model_name} has 0 catch :(')
                # print('-'*20)
                # continue
                pass

            # collect simple spprs:
            # print(f'running {model_name}...')
            results = dict()
            spacial_balance = dict()
            results['SPPR_1986'] = model.SPPR_1986()
            results['SPPR_1995_mTL_global_TE_01'] = model.SPPR_1995(global_TE=0.1)  # This is the 1995 model
            results['SPPR_1995_TL2_global_TE_01'] = model.SPPR_1995_TL_fix(global_TE=0.1) # should be overestimation for omnivorous species
            results['SPPR_1995_mTL_global_mTE'] = model.SPPR_1995(global_TE='mean')  
            results['SPPR_1995_TL2_global_mTE'] = model.SPPR_1995_TL_fix(global_TE='mean') # should be overestimation for omnivorous species
            results['SPPR_1995_TL_Jensened_global_mTE'], _, _ = model.SPPR_EwE_Ulanowicz(TE_option='global', global_TE='mean', use_EE=False)  # by Jensen: should be between SPPR_1995_mTL_global_TE and SPPR_1995_TL2_global_TE
            results['SPPR_1995_TL_TE_Jensened'], _, _ = model.SPPR_EwE_Ulanowicz(TE_option='TE', global_TE=None, use_EE=False)  # by Jensen: should be between SPPR_1995_mTL_global_TE and SPPR_1995_TL2_global_TE
            results['SPPR_EwE'], _, _ = model.SPPR_EwE(TE_option='TE', use_EE=True, return_paths=True, silent=False)  # This is the EwE model
            results['SPPR_2015'], _, _ = model.SPPR_2015(only_pp_det=True)  # this is 2015's model
            results['SPPR_new_2015'], _, _ = model.SPPR_new(TE_option='TE', DET_TE_vals=1)  # should be like SPPR_2015
            results['SPPR_new_full'], _, _ = model.SPPR_new(TE_option='With Egestion', DET_TE_vals=1)  # TE includes egestion
            results['SPPR_new_GE'], _, _  = model.SPPR_new(TE_option='GE', DET_TE_vals=1)  # egestion counts as respiration
            _, results['SPPR_symbolic_GE_import_as_PP'], _, _  = model.SPPR_symbolic(TE_option='GE', diet_import_option='as_PP')
            _, results['SPPR_symbolic_GE_import_as_DC'], e, v  = model.SPPR_symbolic(TE_option='GE', diet_import_option='as_DC')
            spacial_balance['SPPR_symbolic_GE_import_as_DC'] = model.is_sppr_balanced(results['SPPR_symbolic_GE_import_as_DC'], diet_import_equations=(e, v))
            _, results['SPPR_symbolic_TE_import_as_PP'], _, _  = model.SPPR_symbolic(TE_option='TE', diet_import_option='as_PP')
            _, results['SPPR_symbolic_TE_import_as_DC'], e, v  = model.SPPR_symbolic(TE_option='TE', diet_import_option='as_DC')
            spacial_balance['SPPR_symbolic_TE_import_as_DC'] = model.is_sppr_balanced(results['SPPR_symbolic_TE_import_as_DC'], diet_import_equations=(e, v))
            _, results['SPPR_symbolic_full_import_as_PP'], _, _  = model.SPPR_symbolic(TE_option='With Egestion', diet_import_option='as_PP')
            _, results['SPPR_symbolic_full_import_as_DC'], e, v  = model.SPPR_symbolic(TE_option='With Egestion', diet_import_option='as_DC')
            spacial_balance['SPPR_symbolic_full_import_as_DC'] = model.is_sppr_balanced(results['SPPR_symbolic_full_import_as_DC'], diet_import_equations=(e, v))
            
            if collect_mc_sppr:
                # results[f'SPPR_mc_GE_{5}_{10}'], _ = model.monte_carlo_SPPR(
                #     n_samples=n_samples, TE_error_percent=5, TE_error_cut_percent=10, kind='symbolic', TE_option='GE', silent=True
                # )
                results[f'SPPR_mc_GE_{10}_{20}'], _, _, _, _ = model.monte_carlo_SPPR(
                    n_samples=n_samples, TE_error_percent=10, TE_error_cut_percent=20, kind='new', TE_option='GE', silent=False
                )
                # results[f'SPPR_mc_full_{10}_{20}'], _, _, _, _ = model.monte_carlo_SPPR(
                #     n_samples=n_samples, TE_error_percent=10, TE_error_cut_percent=20, kind='new', TE_option='With Egestion', silent=True
                # )
                # results[f'SPPR_mc_GE_{20}_{50}'], _ = model.monte_carlo_SPPR(
                #     n_samples=n_samples, TE_error_percent=20, TE_error_cut_percent=50, kind='symbolic', TE_option='GE', silent=True
                # )

            sanity_checks = dict()
            for method_name, sppr in results.items():
                if method_name not in spacial_balance.keys():
                    is_sppr_balanced, inflow, outflow = model.is_sppr_balanced(sppr)
                else:
                    is_sppr_balanced, inflow, outflow = spacial_balance[method_name]
                is_model_balanced, p, q = model.is_model_balanced()
                sanity_checks[method_name] = {
                    'model is balanced': is_model_balanced,
                    'SPPR did not explode': bool(np.all(sppr >= -1e-10)),
                    'SPPR is balanced': is_sppr_balanced,
                    'inflow': inflow,
                    'outflow': outflow,
                    'DC rows sum to 1': bool(np.all(np.isclose(model._DC.sum(axis=1)[model.get_Regular_seq()], 1, atol=1e-5, rtol=1e-5))),
                    'catch is not negative': all(model.catch >= 0)
                }

            # save results into groups_df:
            groups_df = model.get_groups_df()
            index = groups_df.index
            groups_df['ge'] = model.GE
            groups_df['catch'] = model.catch
            groups_df['p'] = model.p
            groups_df['q'] = model.q
            groups_df['predation'] = model.predation
            groups_df['M0'] = model.M0
            groups_df['egestion'] = model.egestion
            groups_df['biomass_accum'] = model.growth
            for k, v in results.items():
                v = v.reindex(index, fill_value=1)
                if isinstance(v, pd.Series):
                    v = v.to_frame()
                if len(v.columns) > 1:
                    v_orig = v.copy()
                    v['sum all'] = v_orig.sum(axis=1)
                    v['sum PP'] = v_orig[model.get_PP_seq()].sum(axis=1)
                    v['sum inner'] = v_orig.drop(columns=model.get_Import_seq()).sum(axis=1)
                    v = v.rename(columns=model.seq2name)
                    for col in v.columns:
                        groups_df[k + ' (' + str(col) + ')'] = v[col]
                else:
                    groups_df[k + ' (sum all)'] = v.iloc[:, 0]
                    groups_df[k + ' (sum inner)'] = v.iloc[:, 0]
            
            # save:
            with pd.ExcelWriter(f'real_models/{model_name}.xlsx') as writer:
                groups_df.to_excel(writer, sheet_name='sppr')
                pd.DataFrame(sanity_checks).to_excel(writer, sheet_name='sanity_checks')
            # print('-'*20)
            saved_models.append(model_number)

        except KeyboardInterrupt:
            break

        except:
            # print(f'    model {model_number} failed for unknown reason')
            # print('-'*20)
            continue

    failed_models = [i for i in model_numbers if i not in saved_models]
    print(f'failed on models ({len(failed_models)}/{len(model_numbers)}): {failed_models}')
    print(f'success on models ({len(saved_models)}/{len(model_numbers)}): {saved_models}')


if __name__ == "__main__":
    main()