class Solution {
public:
    int orderOfLargestPlusSign(int N, vector<vector<int>>& mines) {
        // create a hashmap of all the mines for easy look up
        unordered_map<int, unordered_set<int>> mine_map;
        for (int i = 0; i < mines.size(); i++) {
            int x = mines[i][0];
            int y = mines[i][1];
            if (mine_map.find(x) != mine_map.end()) mine_map[x].insert(y);
            else mine_map[x] = unordered_set<int>{y};
        }
        
        // initialize counts of contiguous 1's from all directions
        int L = 0;
        int R = 0;
        vector<int> U(N, 0);
        vector<int> D(N, 0);
        
        // go through each row to update counts while keeping track of largest plus sign
        int largest = 0;
        for (int i = 0; i < N; i++) {
            // infer grid value from mine_map
            vector<int> row(N, 1);
            for (int j = 0; j < N; j++) if (mine_map.find(i) != mine_map.end() && mine_map[i].find(j) != mine_map[i].end()) row[j] = 0;
            // update counts from all directions within the i'th row
            for (int j = 0; j < N; j++) {
                if (!row[j]) {
                    L = 0;
                    R = 0;
                    U[j] = 0;
                    D[j] = 0;
                }
                else { // row[j] == 1
                    if (!j) L = 1; // first count from left
                    else L++; // accumulate counts from left
                    
                    if (!R) { // need to find counts from right
                        int k = j;
                        while (k < N && row[k]) k++;
                        R = k - j;
                    }
                    else R--; // counts from right simply decumulates
                    
                    U[j]++; // accumulate counts from up
                    
                    if (!D[j]) { // need to find counts from down
                        int k = i;
                        int m = 1;
                        while (k < N && m) {
                            k++;
                            if (k < N && mine_map.find(k) != mine_map.end() && mine_map[k].find(j) != mine_map[k].end()) m = 0;
                        }
                        D[j] = k - i;
                    }
                    else D[j]--; // counts from down simply decumulates
                }
                
                // plus sign size is minimum over all directions
                int m1 = min(L, R);
                int m2 = min(U[j], D[j]);
                int m = min(m1, m2);
                if (largest < m) largest = m;
            }
        }
        return largest;
    }
};