class Solution {
private:
    int setBoundary(vector<vector<int>>& boundary_states, priority_queue<int, vector<int>, greater<int>>& boundary_heights, unordered_map<int, vector<vector<int>>>& boundary_map, vector<vector<int>>& heightMap, int i, int j, int h) {
        boundary_states[i][j] = 1;
        int w = 0;
        if (heightMap[i][j] < h) {
            w = h - heightMap[i][j];
            heightMap[i][j] = h; // filled water forms effective boundary
        }
        h = heightMap[i][j];
        boundary_heights.emplace(h);
        if (boundary_map.find(h) != boundary_map.end()) boundary_map[h].push_back(vector<int>{i,j});
        else boundary_map.emplace(h,vector<vector<int>>{vector<int>{i,j}});
        return w;
    }
public:
    int trapRainWater(vector<vector<int>>& heightMap) {
        if (heightMap.empty()) return 0;
        int m = heightMap.size();
        int n = heightMap[0].size();
        int h;
        
        // set up initial boundary
        priority_queue<int, vector<int>, greater<int>> boundary_heights;
        unordered_map<int, vector<vector<int>>> boundary_map;
        vector<vector<int>> boundary_states(m, vector<int>(n,0)); // 1 for boundary, 0 for inside, -1 for outside
        boundary_states[0][0] = -1;
        boundary_states[0][n-1] = -1;
        boundary_states[m-1][0] = -1;
        boundary_states[m-1][n-1] = -1;
        for (int j = 1; j < n-1; j++) {
            // top row
            boundary_states[0][j] = 1;
            h = heightMap[0][j];
            boundary_heights.emplace(h);
            if (boundary_map.find(h) != boundary_map.end()) boundary_map[h].push_back(vector<int>{0,j});
            else boundary_map.emplace(h,vector<vector<int>>{vector<int>{0,j}});
            // bottom row
            boundary_states[m-1][j] = 1;
            h = heightMap[m-1][j];
            boundary_heights.emplace(h);
            if (boundary_map.find(h) != boundary_map.end()) boundary_map[h].push_back(vector<int>{m-1,j});
            else boundary_map.emplace(h,vector<vector<int>>{vector<int>{m-1,j}});
        }
        for (int i = 1; i < m-1; i++) {
            // left column
            boundary_states[i][0] = 1;
            h = heightMap[i][0];
            boundary_heights.emplace(h);
            if (boundary_map.find(heightMap[i][0]) != boundary_map.end()) boundary_map[h].push_back(vector<int>{i,0});
            else boundary_map.emplace(heightMap[i][0],vector<vector<int>>{vector<int>{i,0}});
            // right column
            boundary_states[i][n-1] = 1;
            h = heightMap[i][n-1];
            boundary_heights.emplace(h);
            if (boundary_map.find(h) != boundary_map.end()) boundary_map[h].push_back(vector<int>{i,n-1});
            else boundary_map.emplace(h,vector<vector<int>>{vector<int>{i,n-1}});
        }
        
        // shrink boundary until to a point
        int water = 0;
        while (!boundary_heights.empty()) {
            // pop shortest boundary
            h = boundary_heights.top();
            boundary_heights.pop();
            int i = boundary_map[h].back()[0];
            int j = boundary_map[h].back()[1];
            boundary_map[h].pop_back();
            if (boundary_map[h].empty()) boundary_map.erase(h);
            boundary_states[i][j] = -1; // set state to be outside boundary
            // form new boundary
            if (i > 0 && !boundary_states[i-1][j]) water += setBoundary(boundary_states, boundary_heights, boundary_map, heightMap, i-1, j, h);
            if (i < m-1 && !boundary_states[i+1][j]) water += setBoundary(boundary_states, boundary_heights, boundary_map, heightMap, i+1, j, h);
            if (j > 0 && !boundary_states[i][j-1]) water += setBoundary(boundary_states, boundary_heights, boundary_map, heightMap, i, j-1, h);
            if (j < n-1 && !boundary_states[i][j+1]) water += setBoundary(boundary_states, boundary_heights, boundary_map, heightMap, i, j+1, h);
        }
        return water;
    }
};