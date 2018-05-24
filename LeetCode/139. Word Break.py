class Solution:
    def wordBreak(self, s, wordDict):
        """
        :type s: str
        :type wordDict: List[str]
        :rtype: bool
        """
        if not wordDict: return False
        wordDict.sort()
        N = len(wordDict)
        subResult = [True for _ in range(len(s))]
        
        def DPWordBreak(s, wordDict, subResult, index):
            if not subResult[index]: return False
            
            # find bounds in wordDict where the first character matches
            c = s[0]
            lo = 0
            hi = N-1
            while lo < hi:
                mid = (lo + hi) // 2
                if (wordDict[mid][0] < c): lo = mid + 1
                else: hi = mid
            left = hi
            lo = 0
            hi = N-1
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if (wordDict[mid][0] <= c): lo = mid
                else: hi = mid - 1
            right = lo

            # find word matches
            if left <= right:
                i = left
                w = wordDict[i]
                l = len(w)
                while i < right and s[:l] != w: # find the first word that matches
                    i += 1
                    w = wordDict[i]
                    l = len(w)

                while i <= right: # go through all possible word matches
                    w = wordDict[i]
                    l = len(w)
                    if s[:l] == w: 
                        if not s[l:]: return True # finished matching
                        if DPWordBreak(s[l:], wordDict, subResult, index+l): return True # DP
                    i += 1
            
            # no word match
            subResult[index] = False
            return False
        
        return DPWordBreak(s, wordDict, subResult, 0)