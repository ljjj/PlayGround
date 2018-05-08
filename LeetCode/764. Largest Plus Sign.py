class Solution:
    def orderOfLargestPlusSign(self, N, mines):
        """
        :type N: int
        :type mines: List[List[int]]
        :rtype: int
        """
        # create a hashmap of all the mines for easy look up
        mine_map = {}
        for i in range(len(mines)):
            x, y = mines[i]
            if x in mine_map: mine_map[x].add(y)
            else mine_map[x] = set([y])
        
        
        # initialize counts of contiguous 1's from all directions
        L = 0
        R = 0
        U = [0 for _ in range(N)]
        D = [0 for _ in range(N)]
        
        # go through each row to update counts while keeping track of largest plus sign
        largest = 0
        for i in range(N):
            # infer grid value from mine_map
            row = [1 for _ in range(N)]
            for j in range(N):
                if i in mine_map and j in mine_map[i]: row[j] = 0

            # update counts from all directions within the i'th row
            for j in range(N):
                if not row[j]:
                    L = 0
                    R = 0
                    U[j] = 0
                    D[j] = 0
                else:  # row[j] == 1
                    if not j: L = 1 # first count from left
                    else L += 1 # accumulate counts from left
                    
                    if not R:  # need to find counts from right
                        k = j
                        while k < N and row[k]: k += 1
                        R = k - j
                    else: R -= 1 # counts from right simply decumulates
                    
                    U[j] += 1 # accumulate counts from up
                    
                    if not D[j]:  # need to find counts from down
                        k = i
                        m = 1
                        while k < N and m:
                            k += 1
                            if k < N and k in mine_mape and j in mine_map[k]: m = 0
                        D[j] = k - i
                    else: D[j] -= 1 # counts from down simply decumulates
                
                
                # plus sign size is minimum over all directions
                m = min(L, R, U[j], D[j])
                if largest < m: largest = m
            
        return largest