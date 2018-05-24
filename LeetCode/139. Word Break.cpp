class Solution {
private:
    bool DPWordBreak(string s, const vector<string>& wordDict, vector<bool>& subResult, int index) {
        if (!subResult[index]) return false;
        
        // find bounds in wordDict where the first character matches
        int N = wordDict.size();
        char c = s[0];
        int lo = 0;
        int hi = N-1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (wordDict[mid][0] < c) lo = mid + 1;
            else hi = mid;
        }
        int left = hi;
        lo = 0;
        hi = N-1;
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (wordDict[mid][0] <= c) lo = mid;
            else hi = mid - 1;
        }
        int right = lo;

        // find word matches
        if (left <= right) {
            int i = left;
            string w = wordDict[i];
            int l = w.size();
            while (i < right && s.substr(0,l) != w) { // find the first word that matches
                ++i;
                w = wordDict[i];
                l = w.size();
            }

            while (i <= right) { // go through all possible word matches
                w = wordDict[i];
                l = w.size();
                if (s.substr(0,l) == w) {
                    if (s.size() == l) return true; // finished matching
                    if (DPWordBreak(s.substr(l), wordDict, subResult, index+l)) return true; // DP
                }
                ++i;
            }
        }

        // no word match
        subResult[index] = false;
        return false;
    }

public:
    bool wordBreak(string s, vector<string>& wordDict) {
        if (wordDict.empty()) return false;
        sort(wordDict.begin(), wordDict.end());
        vector<bool> subResult(s.size(), true);
        
        return DPWordBreak(s, wordDict, subResult, 0);
    }
};