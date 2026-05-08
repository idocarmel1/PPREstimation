import numpy as np

def find_all_cycles(matrix):
    """
    Phase 1: Enumerates all simple cycles using a depth-first 
    backtracking search[cite: 58, 68].
    """
    n = len(matrix)
    cycles = []

    def dfs(start_node, current_path):
        current_node = current_path[-1]
        
        # Check all possible connections from the current node [cite: 64]
        for neighbor in range(n):
            if matrix[current_node][neighbor] > 0:
                # If we return to the start, a simple cycle is found [cite: 69]
                if neighbor == start_node:
                    cycles.append(list(current_path))
                # Otherwise, continue deeper if neighbor isn't in path [cite: 67]
                elif neighbor not in current_path and neighbor > start_node:
                    dfs(start_node, current_path + [neighbor])

    # Iterate through each node as a potential starting point [cite: 59, 73]
    for i in range(n):
        dfs(i, [i])
    return cycles

def get_circuit_probability(cycle, matrix):
    """
    Calculates the probability of a quantum completing the cycle[cite: 97].
    The probability is the product of link fractions[cite: 103].
    """
    prob = 1.0
    for i in range(len(cycle)):
        u = cycle[i]
        v = cycle[(i + 1) % len(cycle)]
        total_output = np.sum(matrix[u, :])
        if total_output > 0:
            prob *= (matrix[u][v] / total_output) # [cite: 99, 101]
        else:
            return 0.0
    return prob

def remove_cycles(original_matrix):
    """
    Phase 2: Matrix decomposition by grouping cycles into nexuses 
    based on shared critical arcs[cite: 2, 93, 111].
    """
    residual = np.copy(original_matrix).astype(float).T
    cycled_flow_matrix = np.zeros_like(residual)
    
    while True:
        all_cycles = find_all_cycles(residual)
        if not all_cycles:
            break
            
        # 1. Identify critical arc (smallest flow) for each cycle 
        cycle_data = []
        for cycle in all_cycles:
            flows = [residual[cycle[i]][cycle[(i+1)%len(cycle)]] for i in range(len(cycle))]
            min_flow = min(flows)
            # Store the cycle, its critical arc value, and its arc indices
            critical_arc_idx = flows.index(min_flow)
            u, v = cycle[critical_arc_idx], cycle[(critical_arc_idx+1)%len(cycle)]
            cycle_data.append({
                'nodes': cycle,
                'min_flow': min_flow,
                'crit_arc': (u, v),
                'prob': get_circuit_probability(cycle, residual)
            })

        # 2. Identify the smallest critical arc across the entire network [cite: 106]
        global_min_flow = min(c['min_flow'] for c in cycle_data)
        
        # 3. Find cycles sharing this smallest critical arc (the nexus) [cite: 94, 109]
        nexus_cycles = [c for c in cycle_data if c['min_flow'] == global_min_flow]
        total_prob = sum(c['prob'] for c in nexus_cycles)

        # 4. Distribute global_min_flow flow among nexus cycles by probability [cite: 96, 110]
        for cycle_info in nexus_cycles:
            fraction = (cycle_info['prob'] / total_prob) if total_prob > 0 else (1/len(nexus_cycles))
            assigned_flow = global_min_flow * fraction
            
            nodes = cycle_info['nodes']
            for i in range(len(nodes)):
                u, v = nodes[i], nodes[(i+1)%len(nodes)]
                residual[u][v] -= assigned_flow
                cycled_flow_matrix[u][v] += assigned_flow

    return residual.T, cycled_flow_matrix.T

# import utils
# # Example usage with a 3x3 matrix
# # M[i][j] is flow from node i to node j
# M = np.array([
#     [0, 1, 0],
#     [0.01, 0.49, 0.5],
#     [0.65, 0.25, 0.1]
# ])

# acyclic = utils.remove_cycles(M)
# print("Acyclic Matrix:\n", acyclic)
# print()

# acyclic, cycled = remove_cycles(M)
# print("Acyclic Matrix:\n", acyclic)
# print("Cycled Flow Matrix:\n", cycled)