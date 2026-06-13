"""
Debug script to find where the processing differs between APIs
"""
import os, sys
os.chdir(r'c:\Users\idoca\Desktop\אישי\אקדמיה\תואר שני\מחקר\BTN\FishEstimationAI')
sys.path.insert(0, os.getcwd())

from ModelData import ModelData
from PPRCalculator import PPRCalculator
import pandas as pd

# Create both models
model_json = ModelData('real_models/EwE_jsons/227_Iceland_(1950).json')
model_num = ModelData(227)

print("=" * 60)
print("STEP 1: Compare initial groups_data")
print("=" * 60)

# Manually check after PPRCalculator initialization but before full processing
groups_df_json = model_json.groups_data.sort_index(ascending=False).copy()
groups_df_num = model_num.groups_data.sort_index(ascending=False).copy()

# Find common columns
common_cols = sorted(set(groups_df_json.columns) & set(groups_df_num.columns))

print(f"Group 21 catch before processing:")
print(f"  JSON: {groups_df_json.loc[21, 'catch']}")
print(f"  Old:  {groups_df_num.loc[21, 'catch']}")

print(f"\nGroup 22 ({model_json.seq2name.get(22)}) before processing:")
print(f"  JSON: {groups_df_json.loc[22, 'catch']}")
print(f"  Old:  {groups_df_num.loc[22, 'catch']}")
print(f"  JSON tl: {groups_df_json.loc[22, 'tl']}")
print(f"  Old tl: {groups_df_num.loc[22, 'tl']}")

# Check if they're equal on common columns
df_json_common = groups_df_json[common_cols]
df_num_common = groups_df_num[common_cols]
print(f"\nCommon column DataFrames equal: {df_json_common.equals(df_num_common)}")

print("\n" + "=" * 60)
print("STEP 2: Compare after full PPRCalculator creation")
print("=" * 60)

# Create full PPRCalculator objects
ppr_json = PPRCalculator.from_modeldata(model_json)
ppr_num = PPRCalculator.from_modeldata(model_num)

df_json_after = ppr_json._groups_df
df_num_after = ppr_num._groups_df

print(f"Group 21 catch after processing:")
print(f"  JSON: {df_json_after.loc[21, 'catch']}")
print(f"  Old:  {df_num_after.loc[21, 'catch']}")

# Check differences
common_cols_after = sorted(set(df_json_after.columns) & set(df_num_after.columns))
df_json_common_after = df_json_after[common_cols_after]
df_num_common_after = df_num_after[common_cols_after]

print(f"\nCommon column DataFrames equal: {df_json_common_after.equals(df_num_common_after)}")

# Find which columns differ
print("\nColumn-by-column comparison for group 21:")
for col in common_cols_after:
    json_val = df_json_common_after.loc[21, col]
    old_val = df_num_common_after.loc[21, col]
    if pd.isna(json_val) and pd.isna(old_val):
        status = "OK (both NaN)"
    elif json_val == old_val:
        status = "OK"
    else:
        status = f"DIFFER: {json_val} vs {old_val}"
        if not pd.isna(json_val) and not pd.isna(old_val) and isinstance(json_val, (int, float)) and isinstance(old_val, (int, float)):
            pct_diff = abs(json_val - old_val) / max(abs(old_val), abs(json_val), 1e-10) * 100
            status += f" ({pct_diff:.2f}%)"
    print(f"  {col:20s}: {status}")
