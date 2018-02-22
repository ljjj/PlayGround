class RangeModule:

    def __init__(self):
        self.range = []

    def addRange(self, left, right):
        """
        :type left: int
        :type right: int
        :rtype: void
        """
        i = 0
        while i < len(self.range) and left > self.range[i][0]: i += 1 # find where to insert or merge with left
        if i == 0 or left > self.range[i-1][1]: # insert because self.range[i-1][1] < left <= self.range[i][0]
            if i == len(self.range) or right < self.range[i][0]: self.range.insert(i, [left, right]) # pure insert
            else: # need to merge with right after insert
                j = i
                while j < len(self.range) and right > self.range[j][1]: j += 1 # find where to merge with right
                if j == len(self.range) or right < self.range[j][0]: # merge with right not including j
                    self.range[i] = [left, right] # self.range[j-1][1] < right < self.range[j][0]
                    self.range[i+1:j] = []
                else: # merge with right including j
                    self.range[i] = [left, self.range[j][1]] # self.range[j][0] <= right <= self.range[j][1]
                    self.range[i+1:j+1] = []
        else: # merge with left because self.range[i-1][0] < left <= self.range[i-1][1]
            # no need to merge with right
            if i == len(self.range) or right < self.range[i][0]: self.range[i-1][1] = max(right, self.range[i-1][1])
            else: # need to merge with right, similar to above
                j = i
                while j < len(self.range) and right > self.range[j][1]: j += 1
                if j == len(self.range) or right < self.range[j][0]:
                    self.range[i-1][1] = right # self.range[j-1][1] < right < self.range[j][0]
                    self.range[i:j] = []
                else:
                    self.range[i-1][1] = self.range[j][1] # self.range[j][0] <= right <= self.range[j][1]
                    self.range[i:j+1] = []
        
    def queryRange(self, left, right):
        """
        :type left: int
        :type right: int
        :rtype: bool
        """
        i = 0
        while i < len(self.range) and left >= self.range[i][0]: i += 1 # the query range could be in self.range[i-1]
        if i == 0 or left > self.range[i-1][1]: return False
        else:
            if right <= self.range[i-1][1]: return True
            else: return False

    def removeRange(self, left, right):
        """
        :type left: int
        :type right: int
        :rtype: void
        """
        i = 0
        while i < len(self.range) and left > self.range[i][0]: i += 1 # find where to start removal from left
        if i == 0 or left >= self.range[i-1][1]: # remove from i because self.range[i-1][1] < left <= self.range[i][0]
            if i == len(self.range) or right <= self.range[i][0]: pass # nothing to remove
            else:
                j = i
                while j < len(self.range) and right >= self.range[j][1]: j += 1 # find where to remove on the right
                if j < len(self.range) and right > self.range[j][0]: # remove part of j
                    self.range[j][0] = right # self.range[j][0] < right < self.range[j][1]
                self.range[i:j] = []
        else: # remove part of i-1 because self.range[i-1][0] < left < self.range[i-1][1]
            if i == len(self.range) or right <= self.range[i][0]:
                if right < self.range[i-1][1]: self.range.insert(i, [right, self.range[i-1][1]]) # i-1 is broken into two because self.range[i-1][0] < left < right < self.range[i-1][1]
            else:
                j = i
                while j < len(self.range) and right >= self.range[j][1]: j += 1 # find where to remove on the right
                if j < len(self.range) and right > self.range[j][0]: # remove part of j
                    self.range[j][0] = right # self.range[j][0] < right < self.range[j][1]
                self.range[i:j] = []
            self.range[i-1][1] = left # remove part of i-1
        
# Your RangeModule object will be instantiated and called as such:
# obj = RangeModule()
# obj.addRange(left,right)
# param_2 = obj.queryRange(left,right)
# obj.removeRange(left,right)