# approach the tetrahedron puzzle with depth first search
import numpy, sys, copy


# define the tetrahedron class and its operations
class Tetrahedron:
    def __init__(self, i, j, base, left, right, horizontal):
        self.i = i # current row on the grid
        self.j = j # current column on the grid
        self.orientup = (j % 2 == 0) # current orientation of the base triangle
        self.base = base # number on the base triangle
        self.left = left # number on the triangle touching the left edge of the base
        self.right = right # number on the triangle touching the right edge of the base
        self.horizontal = horizontal # number on the triangle touching the horizontal edge of the base

    def rollLeft(self):
        self.j -= 1
        self.orientup = not self.orientup
        temp = self.base
        self.base = self.left
        self.left = self.horizontal
        self.horizontal = self.right
        self.right = temp
        
    def rollRight(self):
        self.j += 1
        self.orientup = not self.orientup
        temp = self.base
        self.base = self.right
        self.right = self.horizontal
        self.horizontal = self.left
        self.left = temp
    
    def rollUp(self):
        if self.orientup is True:
            sys.stderr.write('Cannot roll up at (' + str(self.i) + ', ' + str(self.j) + ')')
            sys.exit(1)
        self.i -= 1
        self.j -= 1
        self.orientup = not self.orientup
        temp = self.base
        self.base = self.horizontal
        self.horizontal = temp
    
    def rollDown(self):
        if self.orientup is False:
            sys.stderr.write('Cannot roll down at (' + str(self.i) + ', ' + str(self.j) + ')')
            sys.exit(1)
        self.i += 1
        self.j += 1
        self.orientup = not self.orientup
        temp = self.base
        self.base = self.horizontal
        self.horizontal = temp


# initialize
M, N = 8, 16 # grid size
grid =[[ 71,  1,  1, 71,711,711, 17,711, 71,  1, 71, 17, 71,711,711,  1],
       [ 71,711, 17, 71,711, 17,  1, 71, 17, 17,711,711, 17,711,  1, 71],
       [ 17,711,  1,711, 71,  1,711, 17,711,  1,711, 17, 71,  1, 17,711],
       [  1, 71,  1, 71, 17,711, 71, 71, 17,711,  1, 17, 71,711, 71,  1],
       [  1, 71,711, 71, 17,  1,711, 17, 71,  1, 71, 17, 71,711, 71, 17],
       [ 71,711,  1, 71, 17,  1,711,  1,711, 71,  1, 71, 17,711, 71,  1],
       [ 71,  1,711, 17, 17,  1,711, 17, 71,  1,711, 17,  1,711,711, 17],
       [  1,711,  1, 17,711, 71,711,  1, 71, 17,711, 71,711, 17, 17,  1]]
visited = numpy.zeros((M, N), dtype=bool) # Boolean to show if this place has been visited
for j in xrange(N):
    if j % 2 == 1:
        visited[0,j] = True
currSum = 0 # current sum
minSum = numpy.inf # lowest sum found so far
currPath = [] # current path in the form: [[i,j,base],...]
optPath = [] # path corresponding to lowest sum


def searchCurrent(tet): # determine how to register the current position
    global grid
#     print tet.i, tet.j, tet.base, tet.left, tet.right, tet.horizontal
    # check if the current position is valid
    if tet.i >= M or tet.j >= N or tet.i < 0 or tet.j < 0:
        sys.stderr.write('Tetrahedron fall off grid')
        sys.exit(1)
    # check if the current tetrahedron base agrees with the grid value
    currGrid = grid[tet.i][tet.j]
#     print currGrid, tet.base
    if currGrid != tet.base:
        return
    # figure out the grid value and the tet top and proceed searching
    # afterwards remove the changes
    searchNext(tet)


def searchNext(tet): # register the current position and roll to next
    global visited, currSum, minSum, currPath, optPath
    # register current to path and sum, mark as visited
    currPath.append([tet.i, tet.j, tet.base])
    currSum += tet.base
    visited[tet.i][tet.j] = True
    # if finished, save to result
    if tet.i == M-1 and tet.j % 2 == 0:
        if currSum < minSum:
            minSum = currSum
            optPath = copy.deepcopy(currPath) # avoid copy by reference
            print optPath
            print(minSum)
    else:
    # roll to possible next
        if tet.j > 0 and not visited[tet.i][tet.j-1]:
            tet.rollLeft()
            searchCurrent(tet)
            tet.rollRight()
        if tet.j < N-1 and not visited[tet.i][tet.j+1]:
            tet.rollRight()
            searchCurrent(tet)
            tet.rollLeft()
        if tet.i > 0 and tet.j > 0 and tet.j % 2 == 1 and not visited[tet.i-1][tet.j-1]:
            tet.rollUp()
            searchCurrent(tet)
            tet.rollDown()
        if tet.i < M-1 and tet.j < N-1 and tet.j % 2 == 0 and not visited[tet.i+1][tet.j+1]:
            tet.rollDown()
            searchCurrent(tet)
            tet.rollUp()
    # remove changes
    currPath.pop()
    currSum -= tet.base
    visited[tet.i][tet.j] = False


for j in xrange(N):
    if j % 2 == 1:
        for base in [1,17,71,711]:
            for left in [1,17,71,711]:
                if left != base:
                    for right in [1,17,71,711]:
                        if right != base and right != left:
                            horizontal = 1 + 17 + 71 + 711 - base - left - right
#                             print 0, j, base, left, right, horizontal
                            tetrahedron = Tetrahedron(0, j, base, left, right, horizontal)
#                             print visited
                            searchCurrent(tetrahedron)


# save the result
import pickle
f = open('Puzzle_2016_Oct_result.pckl','w')
pickle.dump(minSum,f)
pickle.dump(optPath,f)
f.close()

f = open('Puzzle_2016_Oct_result.txt','w')
f.write('Sum: '+str(minSum)+'\n\n')
f.write('Path:\n')
# produce the travel path
travelGrid = numpy.zeros((M,N), dtype = int)
for i in range(len(optPath)):
    travelGrid[optPath[i][0]][optPath[i][1]] = optPath[i][2]
f.write(str(travelGrid))
f.close()
