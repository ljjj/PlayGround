class Solution:
    def maxProfit(self, k, prices):
        """
        :type k: int
        :type prices: List[int]
        :rtype: int
        """
        if len(prices) <= 1: return 0 # no transactions possible
        
        # find all bottom to peak transactions
        transactions = []
        for i in range(1, len(prices)):
            if prices[i] > prices[i-1] and (not transactions or len(transactions[-1]) == 3): transactions.append([prices[i-1]])
            elif prices[i] < prices[i-1] and transactions and len(transactions[-1]) == 1: transactions[-1] = [prices[i-1]-transactions[-1][0], transactions[-1][0], prices[i-1]]
        # last price needs to sell if there is a buy (price guaranteed to be higher than sell)
        i = len(prices)-1
        if prices[i] >= prices[i-1] and transactions and len(transactions[-1]) == 1: transactions[-1] = [prices[i]-transactions[-1][0], transactions[-1][0], prices[i]]
        
        if not transactions: return 0 # no transaction happens
        
        # find all possible losses when a transaction is discarded or when two adjacent transactions merge
        N = t = len(transactions)
        losses = [x[0] for x in transactions] # keep track of loss
        actions = {} # record the action needed for each loss
        # when a transaction is discarded
        for i in range(t):
            loss = transactions[i][0]
            if loss in actions: actions[loss].append([i])
            else: actions[loss] = [[i]]
        # when two adjacent transactions merge
        for i in range(1, t):
            loss = transactions[i-1][2]-transactions[i][1]
            losses.append(loss)
            if loss in actions: actions[loss].append([i-1,i])
            else: actions[loss] = [[i-1,i]]
        heapq.heapify(losses) # keep track of minimal loss
        
        # discard or merge transactions until transaction numbers reach requirement
        while t > k:
            # get minimal loss and corresponding action
            min_loss = heapq.heappop(losses)
            action = actions[min_loss].pop()
            # take action
            if len(action) == 1: # discard a transaction
                d = action[0] # current transaction to discard
                if transactions[d][0] == 0 or transactions[d][0] != min_loss: continue # already discarded or merged
                else:
                    # find adjacent remaining transactions
                    m1 = d-1
                    m2 = d+1
                    if m1 >=0 and transactions[m1][0] == 0: m1 = transactions[m1][1]
                    if m2 < N and transactions[m2][0] == 0: m2 = transactions[m2][2]
                    transactions[d] = [0, m1, m2] # store adjacent remaining transactions for later use
                    t -= 1 # keep track of number of transactions
            else: # merge two adjacent transactions
                m1 = action[0]
                m2 = action[1]
                if transactions[m1][0] == 0 or transactions[m2][0] == 0: continue # already discarded or merged
                else:
                    # find buy and sell price and record to the earlier transaction
                    buy = transactions[m1][1]
                    sell = transactions[m2][2]
                    gain = sell-buy
                    transactions[m1] = [gain, buy, sell]
                    # record loss and action if this merged transaction is discarded
                    heapq.heappush(losses,gain)
                    if gain in actions: actions[gain].append([m1])
                    else: actions[gain] = [[m1]]
                    # find adjacent remaining transactions
                    d = m2
                    m2 = d+1
                    if m2 < N and transactions[m2][0] == 0: m2 = transactions[m2][2]
                    transactions[d] = [0, m1, m2] # store adjacent remaining transactions for later use
                    t -= 1 # keep track of number of transactions
            # record loss and action of if previously found adjacent remaining transactions are merged
            if m1 >= 0 and m2 < N:
                loss = transactions[m1][2]-transactions[m2][1]
                heapq.heappush(losses,loss)
                if loss in actions: actions[loss].append([m1,m2])
                else: actions[loss] = [[m1,m2]]
                # store adjacent remaining transactions for later use
                transactions[m1+1] = [0, m1, m2]
                transactions[m2-1] = [0, m1, m2]
        
        return sum([x[0] for x in transactions]) # sum total profit with required number of transactions