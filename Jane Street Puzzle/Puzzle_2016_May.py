# DFS approach, remember to deep copy the result
import numpy, copy, itertools

def placeHook(n0,x0,y0): # fill in the n0-by-n0 hook
    global validGrids
    if n0 == 0: # done!
        # check sums again
        for i in range(N):
            rs = 0
            for j in range(N):
                rs += grid[i][j]
            if rs != rowSums[i]:
                return
        for j in range(N):
            cs = 0
            for i in range(N):
                cs += grid[i][j]
            if cs != colSums[j]:
                return
        validGrids.append(copy.deepcopy(grid))
    elif n0 == 1: # simply fill in the single digit
        placeDigit([x0],[y0],x0,y0,n0,1,0)
    else:
        for hookOrient in range(4): # loop over hook orientation
            if hookOrient == 0:
                xInd = list(range(x0,x0+n0)) + [x0 for x in range(1,n0)]
                yInd = [y0 for y in range(n0)] + list(range(y0+1,y0+n0))
                xNext = x0+1
                yNext = y0+1
            elif hookOrient == 1:
                xInd = list(range(x0,x0+n0)) + [x0 for x in range(1,n0)]
                yInd = [y0+n0-1 for y in range(n0)] + list(range(y0,y0+n0-1))
                xNext = x0+1
                yNext = y0
            elif hookOrient == 2:
                xInd = list(range(x0,x0+n0)) + [x0+n0-1 for x in range(1,n0)]
                yInd = [y0 for y in range(n0)] + list(range(y0+1,y0+n0))
                xNext = x0
                yNext = y0+1
            elif hookOrient == 3:
                xInd = list(range(x0,x0+n0)) + [x0+n0-1 for x in range(1,n0)]
                yInd = [y0+n0-1 for y in range(n0)] + list(range(y0,y0+n0-1))
                xNext = x0
                yNext = y0
            placeDigit(xInd,yInd,xNext,yNext,n0,1,0) # start searching the first digit at 0 index

def placeDigit(xs,ys,xN,yN,n,k,ind): # search with the k'th digit of n at ind
    global grid
    if ind+n-k >= 2*n-1: # exit search if it's not possible to fill in all digits
        return
    cx = xs[ind]
    cy = ys[ind]
    grid[cx][cy] = n # place the digit
    if checkSums(cx,cy): # only continue to the next digit if this place is valid
        if k == n: # all n digits are placed, go to next hook
            placeHook(n-1,xN,yN)
        else:
            placeDigit(xs,ys,xN,yN,n,k+1,ind+1) # put the next digit of n at the next ind
    grid[cx][cy] = 0 # remove the digit after search
    placeDigit(xs,ys,xN,yN,n,k,ind+1) # search the k'th digit at the next ind

def checkSums(cx,cy): # check if the row and column sums are satisfied
    rs = 0
    for j in range(N):
        rs += grid[cx][j]
    if rs > rowSums[cx]:
        return False
    cs = 0
    for i in range(N):
        cs += grid[i][cy]
    if cs > colSums[cy]:
        return False
    return True

# initialize
print('Initializing ...')
N = 9 # grid size
rowSums = [45,44,4,48,7,14,47,43,33]
colSums = [36,5,47,35,17,30,21,49,45]
##N = 4
##rowSums = [7,11,4,8]
##colSums = [10,3,5,12]
grid = [[0 for x in range(N)] for y in range(N)]
validGrids = []

# find the completed grid
placeHook(N,0,0)

# save the grid
import pickle
f = open('Puzzle_2016_May_result.pckl','wb')
pickle.dump(validGrids,f)
f.close()
f = open('Puzzle_2016_May_result.txt','w')
grid = numpy.zeros((N,N), dtype = int)
for v in range(len(validGrids)):
    for i in range(N):
        for j in range(N):
            grid[i][j] = validGrids[v][i][j]
    f.write(str(grid))
    f.write('\n')
f.close()

# find the largest product
largeP = 0
for v in range(len(validGrids)):
    grid = validGrids[v]
    for p in itertools.permutations(range(N)):
        ys = list(p)
        prod = 1
        for i in range(N):
            prod *= grid[i][ys[i]]
            if prod == 0:
                break
        if prod > largeP:
            largeP = prod
print(largeP)
