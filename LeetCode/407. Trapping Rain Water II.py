from heapq import *

class Solution:
    def setBoundary(self, boundary_states, boundary_heights, heightMap, i, j, h):
        boundary_states[i][j] = 1
        w = 0
        if heightMap[i][j] < h:
            w = h - heightMap[i][j]
            heightMap[i][j] = h # filled water forms effective boundary
        heappush(boundary_heights, (heightMap[i][j],i,j))
        return w
        
    def trapRainWater(self, heightMap):
        """
        :type heightMap: List[List[int]]
        :rtype: int
        """
        if not heightMap: return 0
        m = len(heightMap)
        n = len(heightMap[0])
        
        # set up initial boundary
        boundary_heights = [(heightMap[0][j],0,j) for j in range(1,n-1)] + [(heightMap[m-1][j],m-1,j) for j in range(1,n-1)] + [(heightMap[i][0],i,0) for i in range(1,m-1)] + [(heightMap[i][n-1],i,n-1) for i in range(1,m-1)]
        boundary_states = [[-1] + [1 for j in range(1,n-1)] + [-1]] + [[1] + [0 for j in range(1,n-1)] + [1] for i in range(1,m-1)] + [[-1] + [1 for j in range(1,n-1)] + [-1]] # 1 for boundary, 0 for inside, -1 for outside
        heapify(boundary_heights)
        
        # shrink boundary until to a point
        water = 0
        while boundary_heights:
            h, i, j = heappop(boundary_heights)
            boundary_states[i][j] = -1 # set state to be outside boundary
            # form new boundary
            if i > 0 and not boundary_states[i-1][j]: water += self.setBoundary(boundary_states, boundary_heights, heightMap, i-1, j, h)
            if i < m-1 and not boundary_states[i+1][j]: water += self.setBoundary(boundary_states, boundary_heights, heightMap, i+1, j, h)
            if j > 0 and not boundary_states[i][j-1]: water += self.setBoundary(boundary_states, boundary_heights, heightMap, i, j-1, h)
            if j < n-1 and not boundary_states[i][j+1]: water += self.setBoundary(boundary_states, boundary_heights, heightMap, i, j+1, h)
        
        return water