from collections import deque

class Solution:
    def minWindow(self, s, t):
        """
        :type s: str
        :type t: str
        :rtype: str
        """
        # build hash map for last positions of characters in T, each character key maps to a value of deque
        last_pos = {}
        for c in t:
            if c in last_pos:
                last_pos[c].append(-1)
            else:
                last_pos[c] = deque([-1])
        
        # build hash map to indicate the remaining characters in T that have not been found
        found = {}
        for c in last_pos:
            found[c] = len(last_pos[c])
        
        # build a list of found positions in S of all characters in T
        found_pos = [False for _ in s]
        
        window = "" # result window
        j = -1 # position of the earliest character in S that matches T so that s[j:i+1] is a window
        all_found = False # all characters in T have been found in S
        for i in range(len(s)):
            if s[i] in last_pos: # the character is in T
                # update found status
                if j == -1: j = i
                if not all_found and s[i] in found:
                    found[s[i]] -= 1
                    if not found[s[i]]: del found[s[i]]
                    if not found:
                        all_found = True
                        window = s[j:i+1]
                
                # update new last positions of s[i]
                last = last_pos[s[i]][0]
                last_pos[s[i]].popleft()
                last_pos[s[i]].append(i)
                
                # update found positions
                found_pos[i] = True
                if last >= 0: # guaranteed found[last] == True
                    found_pos[last] = False
                    if j == last: # update window
                        while not found_pos[j]: j += 1
                        if all_found and len(window) > i-j+1: window = s[j:i+1]
        
        return window