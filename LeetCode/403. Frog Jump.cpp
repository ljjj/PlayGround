class Solution {
public:
    bool canCross(vector<int>& stones) {
        int i = 0;
        int jump = 1;
        int N = stones.size();
        vector<int> last(N,0); // last stone jumped from
        vector<vector<int>> choices(N,vector<int>{}); // possible choices of jump units available from this position
        vector<vector<bool>> tried(N,vector<bool>{}); // whether this jump unit has been explored
        choices[0].push_back(1);
        tried[0].push_back(true);
        while (i < stones.size()) {
            int j = i;
            while (j < stones.size() && stones[j]-stones[i] < jump) j++; // find next stone jump units away
            if (j < stones.size() && stones[j]-stones[i] == jump) { // jump from i to j
                if (j == N-1) return true;
                last[j] = i;
                i = j;
                bool kp = false, k0 = false, km = false; // only register new jump possibilities
                for (int k = 0; k < choices[i].size(); k++) {
                    if (choices[i][k] == jump + 1) kp = true;
                    if (choices[i][k] == jump) k0 = true;
                    if (choices[i][k] == jump - 1) km = true;
                }
                if (!kp) choices[i].push_back(jump+1);
                if (!k0) choices[i].push_back(jump);
                if (!km && jump > 1) choices[i].push_back(jump-1);
            }
            // find next jump possibility
            bool exhausted = true;
            while (exhausted) {
                exhausted = false;
                int k;
                for (k = 0; k < choices[i].size(); k++) {
                    if (tried[i][k]) continue;
                    jump = choices[i][k]; // find next jump choice
                    tried[i][k] = true;
                    break;
                }
                if (k == choices[i].size()) exhausted = true;
                if (exhausted) i = last[i]; // backtrack if jump choices exhausted
                if (i == 0) return false; // all possibilities exhausted
            }
        }
    }
};