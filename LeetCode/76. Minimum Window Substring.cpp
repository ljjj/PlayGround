class Solution {
public:
    string minWindow(string s, string t) {
        // build hash map for last positions of characters in T, each character key maps to a value of deque
        unordered_map<char, deque<int>> last_pos;
        for (const auto & c : t) {
            if (last_pos.find(c) != last_pos.end()) last_pos[c].push_back(-1);
            else last_pos.emplace(c, deque<int>{-1});
        }
        
        // build hash map to indicate the remaining characters in T that have not been found
        unordered_map<char, int> found;
        for (auto it = last_pos.begin(); it != last_pos.end(); it++) found.emplace(it->first, it->second.size());
        
        // build a list of found positions in S of all characters in T
        vector<bool> found_pos(s.size(), false);
        
        string window = ""; // result window
        int j = -1; // position of the earliest character in S that matches T so that s[j:i+1] is a window
        bool all_found = false; // all characters in T have been found in S
        for (int i = 0; i < s.size(); i++) {
            if (last_pos.find(s[i]) != last_pos.end()) { // the character is in T
                // update found status
                if (j == -1) j = i;
                if (!all_found && found.find(s[i]) != found.end()) {
                    found[s[i]] -= 1;
                    if (!found[s[i]]) found.erase(s[i]);
                    if (found.empty()) {
                        all_found = true;
                        window = s.substr(j, i+1-j);
                    }
                }
                
                // update new last positions of s[i]
                int last = last_pos[s[i]].front();
                last_pos[s[i]].pop_front();
                last_pos[s[i]].push_back(i);
                
                // update found positions
                found_pos[i] = true;
                if (last >= 0) // guaranteed found[last] == True
                    found_pos[last] = false;
                    if (j == last) { // update window
                        while (!found_pos[j]) j++;
                        if (all_found && window.size() > i+1-j) window = s.substr(j, i+1-j);
                    }
            }
        }
        return window;
    }
};