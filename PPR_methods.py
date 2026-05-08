import sympy as sm
import numpy as np
import pandas as pd
from tqdm import tqdm

from ratios import utils_sympy
import utils


# all methods calculate local PPR (not including import) for the purpose of not counting import twice when I run on the global scale.


def method_1995_species(
        Z: np.array,
        import_vec: np.array,
        DET_idx: int,
        catch_vec: np.array  # catch_vec[DET_idx] = 0
):
    # Initialize vectors to zero if not provided
    import_vec = import_vec if import_vec is not None else np.zeros(Z.rows)

    # add import row and column to Z
    # remove cycles from Z
    # turn base and det rows to 0 (import row is also primary producer)
    # calc TL
    Z = utils_sympy.mat_from_np(utils.remove_cycles(Z))
    TL = sm.matrix2numpy(utils_sympy.TL_from_Z(Z, DET_idx=DET_idx, import_vec=import_vec).evalf()).flatten()

    # define TL = [0.1]_n where n=len(Z)
    TE = 0.1

    # calculate import_j = fraction of import in j's diet
    row_sums = np.array(Z.evalf()).astype(float).sum(axis=1)
    import_fraction = np.where(row_sums != 0, import_vec / (row_sums + import_vec), 0)

    # calc SPPR per species by SPPR = TE**(1-TL) * (1-import_j)
    sppr = (TE ** (1-TL)) * (1-import_fraction)

    # calculate PPR = sum(catch * SPPR)
    return sum(catch_vec * sppr)


def method_1995_groups(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,        species2groups: dict | pd.DataFrame,
        TL_averaging_method: str,
        catch_vec: np.array | pd.DataFrame
):
    """1995's method, aggregating species to groups.

    Args:
        Z (pd.DataFrame): used to calculate TLs.
        species2groups (dict | pd.DataFrame): used to translate species to group.
        TL_averaging_method (str): average, weighted_average, special_average.
        catch_vec (pd.DataFrame): map between species to catch weight.
    """
    # Initialize vectors to zero if not provided
    import_vec = import_vec if import_vec is not None else np.zeros(Z.rows)

    # add import row and column to Z
    # remove cycles from Z
    # turn base and det rows to 0 (import row is also primary producer)
    # calc TL
    Z = utils_sympy.mat_from_np(utils.remove_cycles(Z))
    TL = sm.matrix2numpy(utils_sympy.TL_from_Z(Z, DET_idx=DET_idx, import_vec=import_vec).evalf()).flatten()

    # define TE = [0.1]_n where n=len(Z)
    TE = np.ones_like(Z) * 0.1

    # aggregate TLs to groups

    # calculate import_j = fraction of import in j's diet
    row_sums = np.array(Z.evalf()).astype(float).sum(axis=1)
    import_fraction = np.where(row_sums != 0, import_vec / (row_sums + import_vec), 0)
    # calc SPPR per species by SPPR = TE**(1-TL) * (1-import_j)
    # calculate PPR = sum(catch * SPPR)
    pass


def method_1995_TL_fix(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    def SPPR(meanTL):
        TE = 0.1  # global TE, also possible to use 0.119 from 2020's article.
        TL_fraction = meanTL % 1
        TL_int = int(meanTL)
        sppr = (1-TL_fraction) * (1/TE)**(TL_int-1) + TL_fraction * (1/TE)**(TL_int)
        return sppr
    # Initialize vectors to zero if not provided
    import_vec = import_vec if import_vec is not None else np.zeros(Z.rows)

    # add import row and column to Z
    # remove cycles from Z
    # turn base and det rows to 0 (import row is also primary producer)
    # calc TL
    Z = utils_sympy.mat_from_np(utils.remove_cycles(Z))
    TL = sm.matrix2numpy(utils_sympy.TL_from_Z(Z, DET_idx=DET_idx, import_vec=import_vec).evalf()).flatten()

    # calc SPPR by seperation of TL, assuming each fish eats only 2 TLs under itself.
    sppr = np.array([SPPR(x) for x in TL])

    # calculate import_j = fraction of import in j's diet
    row_sums = np.array(Z.evalf()).astype(float).sum(axis=1)
    import_fraction = np.where(row_sums != 0, import_vec / (row_sums + import_vec), 0)

    # calculate PPR = sum(catch * SPPR * (1-import_j))
    return sum(catch_vec * sppr * (1-import_fraction))


def method_full_structure_no_cycles(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    # remove cycles from Z
    Z = utils_sympy.mat_from_np(utils.remove_cycles(Z))
    # add import row and column to Z
    if import_vec is not None:
        import_vec = utils_sympy.mat_from_np(import_vec)
        Z = Z.row_insert(1, sm.zeros(1, Z.rows))
        import_vec = import_vec.row_insert(0, sm.zeros(1, import_vec.cols))
        Z = Z.col_insert(1, utils_sympy.mat_from_np(import_vec))
    
    # turn base and det rows to 0 (import row is also primary producer)
    if DET_idx is not None:
        Z[DET_idx, :] = sm.zeros(1, Z.rows)

    # define TE = [0.1]_n where n=len(Z)
    TE = np.ones_like(Z) * 0.1

    # calculate SPPR by matrix inversion
    sppr = utils_sympy.SPPR(Z, TE=TE)
    # calculate PPR = sum(catch * SPPR * (1-import_j))
    pass


def method_full_structure_no_cycles_differential_TE(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    # add import row and column to Z
    # remove cycles from Z
    # turn base and det rows to 0 (import row is also primary producer)
    # calculate TE from Z
    # calculate SPPR by matrix inversion
    # calculate PPR = sum(catch * SPPR * (1-import_j))
    pass


def method_full_structure(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    # calculate DC from Z and import
    # define TE = [0.1]_n where n=len(Z)
    # calculate SPPR by matrix inversion
    # calculate PPR = sum(catch * SPPR)
    pass


def method_full_structure_differential_TE(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    # calculate DC from Z and import
    # calculate TE from Z
    # calculate SPPR by matrix inversion
    # calculate PPR = sum(catch * SPPR)
    pass


def method_2015_full_structure_differential_TE_DET_from_pp(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    # calculate DC from Z and import
    # calculate TE from Z
    # calculate SPPR by matrix inversion
    # calculate PPR = sum(catch * SPPR)
    pass


def method_full_structure_differential_TE_DET_in_web(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    # calculate DC from Z and import
    # calculate TE from Z
    # calculate SPPR by matrix inversion
    # calculate PPR = sum(catch * SPPR)
    pass


def method_2015_full_structure_TE_mc_DET_from_pp(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    # calculate DC from Z and import
    # calculate TE from Z
    # calculate SPPR by matrix inversion
    # calculate PPR = sum(catch * SPPR)
    pass


def method_new_full_structure_TE_mc_DET_in_web(
        Z: np.array | pd.DataFrame,
        import_vec: np.array | pd.DataFrame,
        export_vec: np.array | pd.DataFrame,
        growth_vec: np.array | pd.DataFrame,
        DET_idx: int | str,
        catch_vec: np.array | pd.DataFrame
):
    # calculate DC from Z and import
    # calculate TE from Z
    # calculate SPPR by matrix inversion
    # calculate PPR = sum(catch * SPPR)
    pass