# approach the travel puzzle with depth first search
import numpy, sys, copy

# define the dice class and its operations
class Dice:
    def __init__(self):
        self._u = 0 # current number up (top)
        self._d = 0 # current number down (bottom)
        self._l = 0 # current number left
        self._r = 0 # current number right
        self._f = 0 # current number front
        self._b = 0 # current number back
        self.i = 0 # current row on the grid
        self.j = 0 # current column on the grid

    def setTop(self, n): # set the current up face
        self._u = n

    def getTop(self): # get the current up face
        return(self._u)

    def dump(self): # output the faces of the dice
        print('Top   : '+str(self._u))
        print('Bottom: '+str(self._d))
        print('Left  : '+str(self._l))
        print('Right : '+str(self._r))
        print('Front : '+str(self._f))
        print('Back  : '+str(self._b))
        
    def tipLeft(self): # tip the dice left
        temp = self._u
        self._u = self._r
        self._r = self._d
        self._d = self._l
        self._l = temp
        self.j = self.j - 1

    def tipRight(self): # tip the dice right
        temp = self._u
        self._u = self._l
        self._l = self._d
        self._d = self._r
        self._r = temp
        self.j = self.j + 1

    def tipFront(self): # tip the dice front
        temp = self._u
        self._u = self._b
        self._b = self._d
        self._d = self._f
        self._f = temp
        self.i = self.i + 1

    def tipBack(self): # tip the dice back
        temp = self._u
        self._u = self._f
        self._f = self._d
        self._d = self._b
        self._b = temp
        self.i = self.i - 1

# initialize
N = 12 # grid size
grid =[[1,5,4,4,6,1,1,4,1,3,7,5],
       [3,0,0,0,0,0,0,0,0,0,0,1],
       [4,0,6,4,1,8,1,4,2,1,0,3],
       [7,0,1,0,0,0,0,0,0,1,0,2],
       [1,0,1,0,6,1,6,2,0,2,0,1],
       [8,0,4,0,1,0,0,8,0,3,0,5],
       [4,0,2,0,5,0,0,3,0,5,0,2],
       [8,0,5,0,1,1,2,3,0,4,0,6],
       [6,0,1,0,0,0,0,0,0,3,0,6],
       [3,0,6,3,6,5,4,3,4,5,0,1],
       [6,0,0,0,0,0,0,0,0,0,0,3],
       [2,1,6,6,4,5,2,1,1,1,7,1]]
visited = numpy.zeros((N, N), dtype=bool) # Boolean to show if this place has been visited
currScore = 1 # current score
maxScore = 1 # maximum score found so far
currPath = [] # current path in the form: [[i,j,top],...]
optPath = [] # path corresponding to maximum score
optDice = Dice() # dice that achieves maximum score
dice = Dice()

def searchCurrent(dice): # determine how to register current position
    global grid
    # check if current position is valid
    if dice.i >= N or dice.j >= N or dice.i < 0 or dice.j < 0:
        sys.stderr.write('Dice fall off grid')
        sys.exit(1)
    # check if current dice top and grid value is valid
    currGrid = grid[dice.i][dice.j]
    if currGrid != 0 and dice.getTop() != 0 and currGrid != dice.getTop():
        return
    # figure out the grid value and the dice top and proceed searching
    # afterwards remove the changes
    if currGrid != 0 and dice.getTop() != 0:
        searchNext(dice)
    elif currGrid != 0 and dice.getTop() == 0:
        dice.setTop(currGrid)
        searchNext(dice)
        dice.setTop(0)
    elif currGrid == 0 and dice.getTop() != 0:
        grid[dice.i][dice.j] = dice.getTop()
        searchNext(dice)
        grid[dice.i][dice.j] = 0
    else:
        for m in range(9):
            grid[dice.i][dice.j] = m+1
            dice.setTop(m+1)
            searchNext(dice)
            grid[dice.i][dice.j] = 0
            dice.setTop(0)

def searchNext(dice): # register current position and tip to next
    global visited, currScore, maxScore, currPath, optPath, optDice
    # register current to path and score, mark as visited
    currPath.append([dice.i, dice.j, dice.getTop()])
    currScore = currScore * dice.getTop()
    visited[dice.i][dice.j] = True
    # if finished, save to result
    if dice.i == N-1 and dice.j == N-1 and currScore > maxScore:
        maxScore = currScore
        optPath = copy.deepcopy(currPath) # avoid copy by reference
        optDice = copy.deepcopy(dice)
        print(maxScore)
        dice.dump()
    else:
    # tip to possible next
        if dice.j > 1 and not visited[dice.i][dice.j-1]:
            dice.tipLeft()
            searchCurrent(dice)
            dice.tipRight()
        if dice.j < N-1 and not visited[dice.i][dice.j+1]:
            dice.tipRight()
            searchCurrent(dice)
            dice.tipLeft()
        if dice.i > 1 and not visited[dice.i-1][dice.j]:
            dice.tipBack()
            searchCurrent(dice)
            dice.tipFront()
        if dice.i < N-1 and not visited[dice.i+1][dice.j]:
            dice.tipFront()
            searchCurrent(dice)
            dice.tipBack()
    # remove changes
    currPath.pop()
    currScore = currScore / dice.getTop()
    visited[dice.i][dice.j] = False

searchCurrent(dice)

# save the result
import pickle
f = open('Puzzle_2016_Feb_result.pckl','w')
pickle.dump(maxScore,f)
pickle.dump(optDice,f)
pickle.dump(optPath,f)
f.close()

f = open('Puzzle_2016_Feb_result.txt','w')
f.write('Score: '+str(maxScore)+'\n\n')
f.write('Dice:\n')
f.write('  '+str(optDice._b)+'\n')
f.write(str(optDice._l)+' '+str(optDice._u)+' '+str(optDice._r)+'\n')
f.write('  '+str(optDice._f)+'\n')
f.write('  '+str(optDice._d)+'\n\n')
f.write('Path:\n')
# produce the travel path
travelGrid = numpy.zeros((N,N), dtype = int)
for i in range(len(optPath)):
    travelGrid[optPath[i][0]][optPath[i][1]] = optPath[i][2]
f.write(str(travelGrid))
f.close()
