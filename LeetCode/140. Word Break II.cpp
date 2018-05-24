class Solution {
private:
    bool DPWordBreak(string s, const vector<string>& wordDict, vector<bool>& subResult, int index, vector<int>& path, vector<string>& sentences) {
        if (!subResult[index]) return false; // this path doesn't work
        
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
        bool found = false;
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
                    path.push_back(i);
                    if (s.size() == l) {
                        found = true;
                        string sent="";
                        for (int j = 0; j < path.size(); j++) {
                            if (j > 0) sent += " ";
                            sent += wordDict[path[j]];
                        }
                        sentences.push_back(sent); // finished matching
                    }
                    else if (DPWordBreak(s.substr(l), wordDict, subResult, index+l, path, sentences)) found = true; // DP
                    path.pop_back();
                }
                ++i;
            }
        }
        
        if (!found) {
            subResult[index] = false; // no word match, no need to try again
            return false;
        }
        else return true;
    }

public:
    vector<string> wordBreak(string s, vector<string>& wordDict) {
        if (wordDict.empty()) return vector<string>{};
        sort(wordDict.begin(), wordDict.end());
        vector<bool> subResult(s.size(), true);
        vector<int> path{};
        vector<string> sentences{};
        
        DPWordBreak(s, wordDict, subResult, 0, path, sentences);
        return sentences;
    }
};