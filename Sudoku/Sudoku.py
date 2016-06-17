##This program solves a Sudoku puzzle by representing each configuration
##of the 9 places of a certain digit with a single number, and then
##eliminates possibilities everytime a possible configuration is updated.

import copy
import numpy as np

# translate between the configuration state 0~6**6 and the apparent puzzle state [0~8]*9 for a single digit
# the puzzle state: the index shows which block, and the value shows where in the block the digit is
def config2puzzle(c):
    p = []
    p.append(c % 9) # block 0 has 9 possibilities
    c //= 9
    
    temp = c % 6 # block 1 has 6 possilibities
    if temp // 3 >= p[0] // 3:
        temp += 3 # shift down a row
    p.append(temp)
    c //= 6
    
    temp = c % 3  # block 2 has 3 possibilities
    temp += 3*(3-p[0]//3-p[1]//3) # find the row
    p.append(temp)
    c //= 3
    
    temp = c % 6 # block 3 has 6 possilibities
    temp += temp // 2 # the direction of constraint is now vertical, temp is in [0,1,3,4,6,7]
    if temp % 3 >= p[0] % 3:
        temp += 1 # shift right a column
    p.append(temp)
    c //= 6
    
    temp = c % 4 # block 4 has 4 possibilities
    temp += temp // 2 # vertical constraint, temp is in [0,1,3,4]
    if temp // 3 >= p[3] // 3:
        temp += 3 # shift down a row
    if temp % 3 >= p[1] % 3:
        temp += 1 # shift right a column
    p.append(temp)
    c //= 4
    
    temp = c % 2 # block 5 has 2 possibilities
    temp += 3*(3-p[3]//3-p[4]//3) # find the row
    if temp % 3 >= p[2] % 3:
        temp += 1 # shift right a column
    p.append(temp)
    c //= 2
    
    temp = c % 3 # block 6 has 3 possibilities
    temp *= 3 # vertical constraints, temp is in [0,3,6]
    temp += 3- p[0]%3 - p[3]%3 # find the column
    p.append(temp)
    c //= 3
    
    temp = c % 2 # block 7 has 2 possibilities
    temp *= 3 # vertical constraints, temp is in [0,3]
    if temp // 3 >= p[6] // 3:
        temp += 3 # shift down a row
    temp += 3- p[1]%3 - p[4]%3 # find the column
    p.append(temp)
    c //= 2
    assert(c == 0)
    
    temp = 0 # final possibility is fixed
    temp += 3*(3-p[6]//3-p[7]//3) # find the row
    temp += 3- p[2]%3 - p[5]%3 # find the column
    p.append(temp)
    return p


# returns true if the puzzle state is valid
def puzzleValid(p):
    for e in p:
        if e not in range(9):
            return False
    for i in range(3):
        if p[3*i]//3 + p[3*i+1]//3 + p[3*i+2]//3 != 3:
            return False
    for i in range(3):
        if p[i]%3 + p[i+3]%3 + p[i+6]%3 != 3:
            return False
    return True


# reverse steps of config2puzzle
def puzzle2config(p):
    c = 0
    temp = p[7] - (3- p[1]%3 - p[4]%3)
    if temp // 3 >= p[6] // 3:
        temp -= 3
    assert(temp % 3 == 0)
    temp //= 3
    assert(temp in range(2))
    c += temp
    
    c *= 3
    temp = p[6] - (3- p[0]%3 - p[3]%3)
    assert(temp % 3 == 0)
    temp //= 3
    assert(temp in range(3))
    c += temp
    
    c *= 2
    temp = p[5]
    if temp % 3 >= p[2] % 3:
        temp -= 1
    temp -= 3*(3-p[3]//3-p[4]//3)
    assert(temp in range(2))
    c += temp
    
    c *= 4
    temp = p[4]
    if temp % 3 >= p[1] % 3:
        temp -= 1
    if temp // 3 >= p[3] // 3:
        temp -= 3
    temp -= temp // 3
    assert(temp in range(4))
    c += temp
    
    c *= 6
    temp = p[3]
    if temp % 3 >= p[0] % 3:
        temp -= 1
    temp -= temp // 3
    assert(temp in range(6))
    c += temp
    
    c *= 3
    temp = p[2] - 3*(3-p[0]//3-p[1]//3)
    assert(temp in range(3))
    c += temp
    
    c *= 6
    temp = p[1]
    if temp // 3 >= p[0] // 3:
        temp -= 3
    assert(temp in range(6))
    c += temp
    
    c *= 9
    c += p[0]
    return c


# remove possible configurations of digit d by existing puzzle states of all digits
# confs is changed while puzzs is not
def removeConfigs(d, confs, puzzs):
    anyRemoved = False
    for c in confs[d]: # loop over all configurations of d
        p = config2puzzle(c)
        removed = False
        for i in range(9): # check each block i
            for e in range(9): # check each digit
                if puzzs[e][i] >= 0: # the digit exists in the block
                    if (e == d and puzzs[e][i] != p[i]) or (e != d and puzzs[e][i] == p[i]):
                        # d exists not here, or the place is taken by another other digit
                        confs[d].remove(c)
                        removed = True
                        anyRemoved = True
                        break
            if removed:
                break                    
    return anyRemoved


# eliminate possbilities in configurations based on puzzle states
# confs is changed while puzzs is not
def loopRemoveConfigs(confs, puzzs):
    i = 0
    keepLooping = [True]*9
    while any(keepLooping):
        keepLooping[i] = removeConfigs(i, confs, puzzs)
        i = (i+1) % 9


# use DFS to find the solution and record it in the global variable
def fillingSearch(confsIn, puzzsIn):
    global Results
    if any(map(lambda x: len(x)==0, confsIn)): # no solution
        return
    if all(map(lambda x: len(x)==1, confsIn)): # a potential solution is found
        result = presentSolution(np.squeeze(confsIn))
        if all(result.flatten()): # no conflicts in configuration states
            Results.append(result)
        return
    # find the (first) digit with fewest possibilities (>1)
    minL = max(map(len, confsIn)) + 1
    sd = -1
    for i in range(9):
        Li = len(confsIn[i])
        if Li > 1 and Li < minL:
            minL = Li
            sd = i
    for sc in confsIn[sd]: # search into each configuration of this digit
        puzzs = copy.deepcopy(puzzsIn) # puzzle will be updated locally
        confs = copy.deepcopy(confsIn) # configs will be updated locally
        sp = config2puzzle(sc)
        puzzs[sd] = sp # update the puzzle
        loopRemoveConfigs(confs, puzzs) # remove configs based on new puzzle
        fillingSearch(confs, puzzs) # next-level search


# show solution on the Sudoku grid, assuming sol is a list of 9 numbers of configuration states
def presentSolution(sol):
    result = np.zeros((9,9), dtype = int)
    for i in range(9):
        puzz = config2puzzle(sol[i]) # convert solution (as configuration) to puzzle states
        for j in range(9): # convert puzzle states to apparent views
            result[(j//3)*3 + puzz[j]//3][(j%3)*3 + puzz[j]%3] = i + 1
    return result


# initialization
Configs = [range(6**6) for _ in range(9)]
Puzzles = [[-1 for _ in range(9)] for _ in range(9)]
Results = []

# take input of the puzzle in the form of 9 lines of numbers
# the nubmers are separated by white spaces where empty entries are denoted with dots
for testNo in range(4):
    f = open('test'+str(testNo)+'.txt','r')
    for row in range(9):
        line = f.readline()
        # line = raw_input()
        line = line.split()
        if len(line) != 9:
            raise Exception("each line should have 9 entries")
        for col in range(9):
            entry = line[col]
            if entry not in map(str,range(1,10))+['.']:
                raise Exception("each entry should either be a digit or a dot")
            if entry != '.': # translate input to (incomplete) puzzle states
                entry = int(entry) - 1
                Puzzles[entry][(row//3)*3 + col//3] = (row%3)*3 + col%3
    f.close()

    # main loop: kill most configurations based on initial input
    loopRemoveConfigs(Configs, Puzzles)

    # search by trying filling out digits
    fillingSearch(Configs, Puzzles)

    # save solutions
    import pickle
    f = open('solution'+str(testNo)+'.pckl','wb')
    pickle.dump(Results,f)
    f.close()

    f = open('solution'+str(testNo)+'.txt','w')
    for res in Results:
        f.write(str(res))
        f.write('\n')
    f.close()