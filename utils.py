import numpy as np
import sympy as sm
import pandas as pd

def mat_from_np(x, resolution=None):
    """Convert numpy array to sympy Matrix with optional rational resolution."""
    x = np.array(x)

    if resolution is not None:
        if x.ndim == 1:
            rationals_list = [sm.Rational(str(element)).limit_denominator(resolution) for element in x]
        elif x.ndim == 2:
            rationals_list = [[sm.Rational(str(element)).limit_denominator(resolution) for element in row] for row in x]
        else:
            raise Exception(f'array has too many dimensions: {x.ndim}')
    else:
        if x.ndim == 1:
            rationals_list = [sm.Rational(str(element)) for element in x]
        elif x.ndim == 2:
            rationals_list = [[sm.Rational(str(element)) for element in row] for row in x]
        else:
            raise Exception(f'array has too many dimensions: {x.ndim}')
    return sm.Matrix(rationals_list)


def move_scattered_identity(df):
    """
    Identifies rows and columns that form an identity matrix,
    even if they are scattered, and moves them to the top-left.
    """
    identity_pairs = [] # List of (row_label, col_label)
    
    # Iterate through rows to find unit vectors
    for row_label, row in df.iterrows():
        # Check if row has exactly one '1' and the rest are '0'
        if np.sum(row == 1) == 1 and np.sum(row == 0) == (len(row) - 1):
            col_idx = np.where(row == 1)[0][0]
            col_label = df.columns[col_idx]
            identity_pairs.append((row_label, col_label))

    if not identity_pairs:
        return df, df.index, df.columns

    # Separate the pairs into row and column lists
    id_rows = [p[0] for p in identity_pairs]
    id_cols = [p[1] for p in identity_pairs]

    # Create the new order (Identity components first)
    new_row_order = id_rows + [r for r in df.index if r not in id_rows]
    new_col_order = id_cols + [c for c in df.columns if c not in id_cols]

    # Apply to both DataFrames
    return df.reindex(index=new_row_order, columns=new_col_order), new_row_order, new_col_order

def _find_all_cycles(matrix):
    n = len(matrix)
    cycles = []

    def dfs(start_node, current_path):
        current_node = current_path[-1]
        for neighbor in range(n):
            if matrix[current_node][neighbor] > 0:
                if neighbor == start_node:
                    cycles.append(list(current_path))
                elif neighbor not in current_path and neighbor > start_node:
                    dfs(start_node, current_path + [neighbor])

    for i in range(n):
        dfs(i, [i])
    return cycles


def _get_circuit_probability(cycle, matrix):
    prob = 1.0
    for i in range(len(cycle)):
        u = cycle[i]
        v = cycle[(i + 1) % len(cycle)]
        total_output = np.sum(matrix[u, :])
        if total_output > 0:
            prob *= matrix[u][v] / total_output
        else:
            return 0.0
    return prob


def _remove_cycles_nexus(original_matrix):
    """Ulanowicz nexus-based cycle removal (Phases 1+2)."""
    residual = np.copy(original_matrix).astype(float).T
    cycled_flow_matrix = np.zeros_like(residual)

    while True:
        all_cycles = _find_all_cycles(residual)
        if not all_cycles:
            break

        cycle_data = []
        for cycle in all_cycles:
            flows = [residual[cycle[i]][cycle[(i+1)%len(cycle)]] for i in range(len(cycle))]
            min_flow = min(flows)
            critical_arc_idx = flows.index(min_flow)
            u, v = cycle[critical_arc_idx], cycle[(critical_arc_idx+1)%len(cycle)]
            cycle_data.append({
                'nodes': cycle,
                'min_flow': min_flow,
                'crit_arc': (u, v),
                'prob': _get_circuit_probability(cycle, residual)
            })

        global_min_flow = min(c['min_flow'] for c in cycle_data)
        nexus_cycles = [c for c in cycle_data if c['min_flow'] == global_min_flow]
        total_prob = sum(c['prob'] for c in nexus_cycles)

        for cycle_info in nexus_cycles:
            fraction = (cycle_info['prob'] / total_prob) if total_prob > 0 else (1 / len(nexus_cycles))
            assigned_flow = global_min_flow * fraction
            nodes = cycle_info['nodes']
            for i in range(len(nodes)):
                u, v = nodes[i], nodes[(i+1)%len(nodes)]
                residual[u][v] -= assigned_flow
                cycled_flow_matrix[u][v] += assigned_flow

    return residual.T, cycled_flow_matrix.T


def remove_cycles(Z_input, new=False):
    if new:
        return remove_cycles_new(Z_input)
    else:
        return remove_cycles_old(Z_input)

def remove_cycles_new(Z_input):
    # Work on a copy to avoid mutating the original matrix
    Z = Z_input.copy()
    
    # 1. Identify input type and convert to numpy for the logic
    is_pandas = isinstance(Z_input, pd.DataFrame)
    is_sm_matrix = isinstance(Z, sm.Matrix)
    if is_pandas:
        Z = Z_input.to_numpy(dtype=float, copy=True)
    elif is_sm_matrix:
        Z = sm.matrix2numpy(Z)
    else:
        Z = np.array(Z_input, dtype=float, copy=True)

    Z, _ = _remove_cycles_nexus(Z)

    # 3. Re-wrap in DataFrame if necessary
    if is_pandas:
        return pd.DataFrame(Z, index=Z_input.index, columns=Z_input.columns)
    if is_sm_matrix:
        return mat_from_np(Z)
    
    return Z

def remove_cycles_old(Z_input):
    """
    Removes cycles from a flow matrix Z using the Ulanowicz method.
    Z[i, j] represents flow from prey i to consumer j.
    """
    # Work on a copy to avoid mutating the original matrix
    Z = Z_input.copy()
    
    # 1. Identify input type and convert to numpy for the logic
    is_pandas = isinstance(Z_input, pd.DataFrame)
    is_sm_matrix = isinstance(Z, sm.Matrix)
    if is_pandas:
        Z = Z_input.to_numpy(dtype=float, copy=True)
    elif is_sm_matrix:
        Z = sm.matrix2numpy(Z)
    else:
        Z = np.array(Z_input, dtype=float, copy=True)
    
    def find_cycle(matrix):
        """Finds a single cycle using DFS."""
        num_nodes = matrix.shape[0]
        visited = [False] * num_nodes
        stack = []
        parent = [-1] * num_nodes

        def dfs(u):
            visited[u] = True
            stack.append(u)
            for v in range(num_nodes):
                if matrix[u, v] > 0:
                    if v in stack: # Cycle detected
                        cycle_path = stack[stack.index(v):]
                        return cycle_path
                    if not visited[v]:
                        res = dfs(v)
                        if res: return res
            stack.pop()
            return None

        for i in range(num_nodes):
            if not visited[i]:
                cycle = dfs(i)
                if cycle: return cycle
        return None

    while True:
        cycle = find_cycle(Z)
        if not cycle:
            break  # No more cycles found
        
        # 1. Find the "Weakest Link" (Minimum flow in the cycle)
        # Create pairs of (from_node, to_node) for the cycle
        links = []
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            links.append((u, v))
            
        min_flow = min(Z[u, v] for u, v in links)
        
        # 2. Subtract the bottleneck flow from all links in the cycle
        for u, v in links:
            Z[u, v] -= min_flow
            # Clean up tiny floating point errors
            if Z[u, v] < 1e-20:
                Z[u, v] = 0

    # 3. Re-wrap in DataFrame if necessary
    if is_pandas:
        return pd.DataFrame(Z, index=Z_input.index, columns=Z_input.columns)
    if is_sm_matrix:
        return mat_from_np(Z)

    return Z


