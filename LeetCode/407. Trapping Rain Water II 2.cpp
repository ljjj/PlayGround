class Solution {
private:
    typedef priority_queue<vector<int>, vector<vector<int>>, greater<vector<int>>> pq;
    int setBoundary(vector<vector<int>>& boundary_states, pq& boundary_heights, vector<vector<int>>& heightMap, int i, int j, int h) {
        boundary_states[i][j] = 1;
        int w = 0;
        if (heightMap[i][j] < h) {
            w = h - heightMap[i][j];
            heightMap[i][j] = h; // filled water forms effective boundary
        }
        boundary_heights.push(vector<int>{heightMap[i][j],i,j});
        return w;
    }
public:
    int trapRainWater(vector<vector<int>>& heightMap) {
        if (heightMap.empty()) return 0;
        int m = heightMap.size();
        int n = heightMap[0].size();
        
        // set up initial boundary
        pq boundary_heights;
        vector<vector<int>> boundary_states(m, vector<int>(n,0)); // 1 for boundary, 0 for inside, -1 for outside
        boundary_states[0][0] = -1;
        boundary_states[0][n-1] = -1;
        boundary_states[m-1][0] = -1;
        boundary_states[m-1][n-1] = -1;
        for (int j = 1; j < n-1; j++) {
            boundary_states[0][j] = 1;
            boundary_heights.push(vector<int>{heightMap[0][j],0,j});
            boundary_states[m-1][j] = 1;
            boundary_heights.push(vector<int>{heightMap[m-1][j],m-1,j});
        }
        for (int i = 1; i < m-1; i++) {
            boundary_states[i][0] = 1;
            boundary_heights.push(vector<int>{heightMap[i][0],i,0});
            boundary_states[i][n-1] = 1;
            boundary_heights.push(vector<int>{heightMap[i][n-1],i,n-1});
        }
        
        // shrink boundary until to a point
        int water = 0;
        while (!boundary_heights.empty()) {
            vector<int> t = boundary_heights.top();
            boundary_heights.pop();
            int h = t[0], i = t[1], j = t[2];
            boundary_states[i][j] = -1; // set state to be outside boundary
            // form new boundary
            if (i > 0 && !boundary_states[i-1][j]) water += setBoundary(boundary_states, boundary_heights, heightMap, i-1, j, h);
            if (i < m-1 && !boundary_states[i+1][j]) water += setBoundary(boundary_states, boundary_heights, heightMap, i+1, j, h);
            if (j > 0 && !boundary_states[i][j-1]) water += setBoundary(boundary_states, boundary_heights, heightMap, i, j-1, h);
            if (j < n-1 && !boundary_states[i][j+1]) water += setBoundary(boundary_states, boundary_heights, heightMap, i, j+1, h);
        }
        return water;
    }
};