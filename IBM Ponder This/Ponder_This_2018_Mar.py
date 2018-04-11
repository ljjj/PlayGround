# initialization
debug = False
vowels = 5 # number of vowels
switches = 6 # number of switches
states = [-1 for _ in range(2**switches)] # assign each state of the switches to a certain vowel
masks = [1 << m for m in range(switches)] # used to toggle switches

# find all neighbors to constrain search
all_neighbors = [sorted([i ^ m for m in masks]) for i in range(2**switches)]

# approach: DFS, backtracking
def dfs(states, k): # assuming states[k] == -1
    if k == 2**switches: # check result validity
        for i in range(2**switches):
            covered = set([states[n] for n in all_neighbors[i]] + [states[i]]) # assuming all states filled
            if len(covered) != vowels: return False
        return True
    
    allowed = set(range(vowels))
    
    # constrain based on previous states i
    for i in [x for x in all_neighbors[k] if x < k]:
        # find which states are covered by neighbors of i
        covered = [states[n] for n in all_neighbors[i]] + [states[i]]
        missing = covered.count(-1)
        covered = set(covered)
        covered.remove(-1) # missing >= 1 because states[k] == -1
        
        # cannot cover all vowels even if all missing neighbors fill correctly
        if len(covered) + missing < vowels:
            if debug: print(k, states, allowed, covered, 'backtrack')
            return False
        
        # k (and all missing neighbors of i) must cover the missing vowels
        if len(covered) + missing == vowels:
            allowed &= set(range(vowels)) - covered
        
        # no states[k] value can satisfy all constraints
        if not allowed:
            if debug: print(k, states, allowed, covered, 'backtrack')
            return False
    
    if debug: print(k, states, allowed, 'next')
    # go through all constrained possible values
    for a in allowed:
        states[k] = a
        if dfs(states, k+1): return True
    
    # no solution exists in this branch
    states[k] = -1
    return False

dfs(states, 0)
result = "".join(["AEIOU"[s] for s in states])