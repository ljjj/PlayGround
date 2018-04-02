class Solution:
    def canCross(self, stones):
        """
        :type stones: List[int]
        :rtype: bool
        """
        i = 0
        jump = 1
        N = len(stones)
        last = [0] * N # last stone jumped from
        choices = [[] for _ in range(N)] # possible choices of jump units available from this position
        tried = [[] for _ in range(N)] # whether this jump unit has been explored
        choices[0].append(1)
        tried[0].append(True)
        while i < N:
            j = i
            while j < N and stones[j]-stones[i] < jump: j += 1 # find next stone jump units away
            if j < N and stones[j]-stones[i] == jump: # jump from i to j
                if j == N-1: return True
                last[j] = i
                i = j
                # only register new jump possibilities
                kp = False
                k0 = False
                km = False
                for k in range(len(choices[i])):
                    if choices[i][k] == jump + 1: kp = True
                    if choices[i][k] == jump: k0 = True
                    if choices[i][k] == jump - 1: km = True
                if not kp:
                    choices[i].append(jump+1)
                    tried[i].append(False)
                if not k0:
                    choices[i].append(jump)
                    tried[i].append(False)
                if not km and jump > 1:
                    choices[i].append(jump-1)
                    tried[i].append(False)
            # find next jump possibility
            exhausted = True
            while exhausted:
                exhausted = False
                k = 0
                while k < len(choices[i]):
                    if tried[i][k]:
                        k += 1
                        continue
                    jump = choices[i][k] # find next jump choice
                    tried[i][k] = True
                    break
                if k == len(choices[i]): exhausted = True
                if exhausted: i = last[i] # backtrack if jump choices exhausted
                if i == 0: return False # all possibilities exhausted