class Solution:
    def wordBreak(self, s, wordDict):
        """
        :type s: str
        :type wordDict: List[str]
        :rtype: List[str]
        """
        if not wordDict: return []
        wordDict.sort()
        N = len(wordDict)
        subResult = [True for _ in range(len(s))]
        path = []
        sentences = []
        
        def DPWordBreak(s, wordDict, subResult, index, path, sentences):
            if not subResult[index]: return False # this path doesn't work
            
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
            found = False
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
                        path.append(i)
                        if not s[l:]:
                            found = True
                            sentences.append(" ".join([wordDict[j] for j in path])) # finished matching
                        elif DPWordBreak(s[l:], wordDict, subResult, index+l, path, sentences): found = True # DP
                        path.pop()
                    i += 1
            
            if not found:                
                subResult[index] = False # no word match, no need to try again
                return False
            else: return True # some match is found
        
        DPWordBreak(s, wordDict, subResult, 0, path, sentences)
        return sentences