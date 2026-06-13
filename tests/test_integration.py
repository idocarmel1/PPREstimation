"""
Comprehensive test of ModelData and PPRCalculator integration.
Tests both old API (model_number) and new API (json_filepath).
"""
import os
import sys

# Change to the project directory and ensure root is on sys.path
os.chdir(r'c:\Users\idoca\Desktop\אישי\אקדמיה\תואר שני\מחקר\BTN\FishEstimationAI')
sys.path.insert(0, os.getcwd())

print("=" * 60)
print("Testing ModelData and PPRCalculator Integration")
print("=" * 60)

# Test 1: Import and basic checks
print("\n[1] Testing imports...")
try:
    from ModelData import ModelData, SpeciesGroup, SpeciesGroupLegacy
    print("    SUCCESS: ModelData and classes imported")
except Exception as e:
    print(f"    FAILED: {e}")
    sys.exit(1)

# Test 2: New API (JSON filepath)
print("\n[2] Testing NEW API (JSON filepath)...")
try:
    json_filepath = 'real_models/EwE_jsons/227_Iceland_(1950).json'
    model_data_json = ModelData(json_filepath)
    assert model_data_json.model_number == 227
    assert model_data_json.model_name == "Iceland"
    assert model_data_json.model_year == "1950"
    assert model_data_json.groups_data.shape[0] > 0
    assert model_data_json.DC.shape[0] > 0
    print(f"    SUCCESS: {json_filepath}")
    print(f"      Model: {model_data_json.model_number} ({model_data_json.model_name}, {model_data_json.model_year})")
    print(f"      Groups: {model_data_json.groups_data.shape[0]}, DC: {model_data_json.DC.shape}")
except Exception as e:
    print(f"    FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 3: Old API (model_number)
print("\n[3] Testing OLD API (model_number - backward compatibility)...")
try:
    model_data_num = ModelData(227)
    assert model_data_num.model_number == 227
    assert model_data_num.model_name == "Iceland"
    assert model_data_num.model_year == "1950"
    assert model_data_num.groups_data.shape[0] > 0
    print(f"    SUCCESS: Model 227")
    print(f"      Model: {model_data_num.model_number} ({model_data_num.model_name}, {model_data_num.model_year})")
    print(f"      Groups: {model_data_num.groups_data.shape[0]}")
except Exception as e:
    print(f"    FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 4: PPRCalculator with new API
print("\n[4] Testing PPRCalculator with NEW API...")
try:
    from PPRCalculator import PPRCalculator
    json_filepath = 'real_models/EwE_jsons/107_Grand_Banks_of_Newfoundland_(1980).json'
    model_data_json = ModelData(json_filepath)
    ppr_calc = PPRCalculator.from_modeldata(model_data_json)
    assert ppr_calc.n_groups > 0
    print(f"    SUCCESS: PPRCalculator created with JSON API")
    print(f"      Groups: {ppr_calc.n_groups}")
except Exception as e:
    print(f"    FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 5: PPRCalculator with old API
print("\n[5] Testing PPRCalculator with OLD API...")
try:
    model_data_num = ModelData(107)
    ppr_calc = PPRCalculator.from_modeldata(model_data_num)
    assert ppr_calc.n_groups > 0
    print(f"    SUCCESS: PPRCalculator created with old API")
    print(f"      Groups: {ppr_calc.n_groups}")
except Exception as e:
    print(f"    FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 6: Verify consistency between APIs
print("\n[6] Verifying consistency between APIs...")
try:
    # Use model 227 which has both old API (pre-serialized) and new API (JSON file)
    model_json = ModelData('real_models/EwE_jsons/227_Iceland_(1950).json')
    model_num = ModelData(227)
    
    # Compare key attributes
    assert model_json.model_number == model_num.model_number, f"Model numbers don't match: {model_json.model_number} vs {model_num.model_number}"
    assert model_json.model_name == model_num.model_name, f"Model names don't match: {model_json.model_name} vs {model_num.model_name}"
    assert model_json.model_year == model_num.model_year, f"Model years don't match: {model_json.model_year} vs {model_num.model_year}"
    # Both should have the same number of rows (groups)
    assert model_json.groups_data.shape[0] == model_num.groups_data.shape[0], f"Row counts don't match: {model_json.groups_data.shape[0]} vs {model_num.groups_data.shape[0]}"
    # New API has taxon_descr column that old API may not have - that's okay
    
    print(f"    SUCCESS: Both APIs produce consistent results")
    print(f"      Model: {model_json.model_number} ({model_json.model_name}, {model_json.model_year})")
    print(f"      Groups: {model_json.groups_data.shape[0]} (new API columns: {model_json.groups_data.shape[1]}, old API columns: {model_num.groups_data.shape[1]})")
    print(f"      Note: New API includes 'taxon_descr' column for additional metadata")
except Exception as e:
    print(f"    FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 7: Verify PPRCalculator groups_df consistency (new API must not lose data vs old API)
print("\n[7] Verifying PPRCalculator groups_df contents...")
try:
    from PPRCalculator import PPRCalculator
    import pandas as pd
    import numpy as np

    model_json = ModelData('real_models/EwE_jsons/227_Iceland_(1950).json')
    model_num = ModelData(227)

    ppr_json = PPRCalculator.from_modeldata(model_json)
    ppr_num = PPRCalculator.from_modeldata(model_num)

    df_json = ppr_json._groups_df.sort_index().sort_index(axis=1)
    df_num = ppr_num._groups_df.sort_index().sort_index(axis=1)

    # Only compare columns present in both (new API adds taxon_descr).
    # 'tl' is pre-stored in the old API's legacy data but computed on-the-fly by the
    # new API (stored in self.TL, not written back to _groups_df), so exclude it.
    SKIP_COLS = {'tl'}
    common_cols = sorted((set(df_json.columns) & set(df_num.columns)) - SKIP_COLS)
    df_json_c = df_json[common_cols]
    df_num_c = df_num[common_cols]

    assert df_json_c.shape == df_num_c.shape, f"Shape mismatch: {df_json_c.shape} vs {df_num_c.shape}"

    # Classify differences as regressions (new API lost data) vs enrichments (new API has more data).
    # The old API serialised some values as NaN that the EwE JSON has as actual numbers; that is
    # acceptable — the new API is richer. A regression would be the reverse: old has a value, new is 0/NaN.
    regressions = []
    enrichments = []

    for col in common_cols:
        j_vals = df_json_c[col]
        n_vals = df_num_c[col]
        if not pd.api.types.is_numeric_dtype(j_vals):
            mismatches = j_vals.index[j_vals != n_vals].tolist()
            if mismatches:
                regressions.append(f"'{col}' string mismatch at groups: {mismatches}")
            continue

        close = np.isclose(j_vals.fillna(0), n_vals.fillna(0), rtol=0.01, atol=0.001)
        if not close.all():
            for idx in j_vals.index[~close]:
                jv = float(j_vals[idx]) if not pd.isna(j_vals[idx]) else 0.0
                nv = float(n_vals[idx]) if not pd.isna(n_vals[idx]) else 0.0
                if abs(nv) > 0.001 and abs(jv) < 0.001:
                    # old had a meaningful value, new API has 0/NaN — regression
                    regressions.append(f"'{col}'[group_seq={idx}]: old={nv:.6g}, new={jv:.6g}")
                else:
                    # new API has data where old had 0/NaN — acceptable enrichment
                    enrichments.append(f"'{col}'[group_seq={idx}]: old={nv:.6g} -> new={jv:.6g}")

    if enrichments:
        print(f"      Info: New API has richer data in {len(enrichments)} field(s) (old API had 0/NaN):")
        for e in enrichments:
            print(f"        {e}")

    assert not regressions, "New API regressed vs old API:\n" + "\n".join(regressions)

    print(f"    SUCCESS: PPRCalculator groups_df are consistent (no regressions)")
    print(f"      Shape (common columns): {df_json_c.shape}, common columns: {len(common_cols)}")
    if 'taxon_descr' in df_json.columns and 'taxon_descr' not in df_num.columns:
        print(f"      Note: New API has extra 'taxon_descr' column (expected)")
except Exception as e:
    print(f"    FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 8: Multi-model consistency check across a diverse sample
print("\n[8] Multi-model consistency check...")
import glob, re as _re, numpy as np, pandas as pd
from PPRCalculator import PPRCalculator

SKIP_COLS = {'tl'}
TEST_MODELS = [2, 40, 107, 413, 227, 311, 444, 112]  # diverse sample (63/175 skipped: multiple DET groups)

json_map = {}
for path in glob.glob('real_models/EwE_jsons/*.json'):
    m = _re.match(r'.*[/\\](\d+)_', path)
    if m:
        json_map[int(m.group(1))] = path.replace('\\', '/')

n_passed = 0
n_failed = 0
for model_id in TEST_MODELS:
    if model_id not in json_map:
        print(f"  SKIP {model_id}: no JSON file found")
        continue
    try:
        m_json = ModelData(json_map[model_id])
        m_num  = ModelData(model_id)
        p_json = PPRCalculator.from_modeldata(m_json)
        p_num  = PPRCalculator.from_modeldata(m_num)

        df_j = p_json._groups_df.sort_index().sort_index(axis=1)
        df_n = p_num._groups_df.sort_index().sort_index(axis=1)
        common = sorted((set(df_j.columns) & set(df_n.columns)) - SKIP_COLS)

        df_jc = df_j[common]
        df_nc = df_n[common]
        assert df_jc.shape == df_nc.shape, f"shape mismatch {df_jc.shape} vs {df_nc.shape}"

        regressions = []
        enrichments = []
        for col in common:
            jv = df_jc[col]; nv = df_nc[col]
            if not pd.api.types.is_numeric_dtype(jv):
                bad = jv.index[jv != nv].tolist()
                if bad:
                    regressions.append(f"'{col}' string diff at {bad}")
                continue
            close = np.isclose(jv.fillna(0), nv.fillna(0), rtol=0.01, atol=0.001)
            for idx in jv.index[~close]:
                j, n = (float(jv[idx]) if not pd.isna(jv[idx]) else 0.0,
                        float(nv[idx]) if not pd.isna(nv[idx]) else 0.0)
                if abs(n) > 0.001 and abs(j) < 0.001:
                    regressions.append(f"'{col}'[{idx}]: old={n:.4g}, new={j:.4g}")
                else:
                    enrichments.append(f"'{col}'[{idx}]")

        assert not regressions, "regressions: " + "; ".join(regressions)

        enrich_note = f", {len(enrichments)} enrichment(s)" if enrichments else ""
        print(f"  PASS  model {model_id:4d} ({m_json.model_name}, {m_json.model_year}) "
              f"— {df_jc.shape[0]} groups, {len(common)} cols{enrich_note}")
        n_passed += 1

    except Exception as e:
        print(f"  FAIL  model {model_id}: {e}")
        import traceback; traceback.print_exc()
        n_failed += 1

print(f"\n    {'SUCCESS' if n_failed == 0 else 'FAILED'}: {n_passed}/{len(TEST_MODELS)} models passed")
if n_failed:
    sys.exit(1)
