class Solution {
public:
    int maxProfit(int k, vector<int>& prices) {
        if (prices.size() <= 1) return 0; // no transactions possible
        
        // find all bottom to peak transactions
        vector<vector<int>> transactions; // each entry of transactions is a vector<int>{gain, buy, sell}
        int i;
        for (i=1; i<prices.size(); i++) {
            // first increase in price means price bottom
            if (prices[i] > prices[i-1] && (transactions.empty() || transactions.back().size() == 3)) transactions.push_back(vector<int>{prices[i-1]});
            // first decrease in price means price top
            else if (prices[i] < prices[i-1] && !transactions.empty() && transactions.back().size() == 1) transactions.back() = vector<int>{prices[i-1]-transactions.back()[0], transactions.back()[0], prices[i-1]};
        }
        // last price needs to sell if there is a buy (price guaranteed to be higher than sell)
        i = prices.size()-1;
        if (prices[i] >= prices[i-1] && !transactions.empty() && transactions.back().size() == 1) transactions.back() = vector<int>{prices[i]-transactions.back()[0], transactions.back()[0], prices[i]};
        
        if (transactions.empty()) return 0; // no transaction happens
        
        // find all possible losses when a transaction is discarded or when two adjacent transactions merge
        int t = transactions.size();
        priority_queue<int> losses; // keep track of minimal loss
        map<int, vector<vector<int>>> actions; // record the action needed for each loss
        // when a transaction is discarded
        for (i=0; i<t; i++) {
            int loss = -transactions[i][0];
            losses.push(loss);
            if (actions.find(loss) != actions.end()) actions[loss].push_back(vector<int>{i});
            else actions[loss] = vector<vector<int>>{vector<int>{i}};
        }
        // when two adjacent transactions merge
        for (i=1; i<t; i++) {
            int loss = -(transactions[i-1][2]-transactions[i][1]);
            losses.push(loss);
            if (actions.find(loss) != actions.end()) actions[loss].push_back(vector<int>{i-1,i});
            else actions[loss] = vector<vector<int>>{vector<int>{i-1,i}};
        }
        
        // discard or merge transactions until transaction numbers reach requirement
        while (t > k) {
            // get minimal loss and corresponding action
            int min_loss = losses.top();
            losses.pop();
            vector<int> action = actions[min_loss].back();
            actions[min_loss].pop_back();
            // take action
            int d, m1, m2;
            if (action.size() == 1) { // discard a transaction
                d = action[0]; // current transaction to discard
                if (transactions[d][0] == 0 || transactions[d][0] != -min_loss) continue; // already discarded or merged
                else {
                    // find adjacent remaining transactions
                    m1 = d-1;
                    m2 = d+1;
                    if (m1 >=0 && transactions[m1][0] == 0) m1 = transactions[m1][1];
                    if (m2 < transactions.size() && transactions[m2][0] == 0) m2 = transactions[m2][2];
                    transactions[d] = vector<int>{0, m1, m2}; // store adjacent remaining transactions for later use
                    t--; // keep track of number of transactions
                }
            }
            else { // merge two adjacent transactions
                m1 = action[0];
                m2 = action[1];
                if (transactions[m1][0] == 0 || transactions[m2][0] == 0) continue; // already discarded or merged
                else {
                    // find buy and sell price and record to the earlier transaction
                    int buy = transactions[m1][1];
                    int sell = transactions[m2][2];
                    int gain = sell-buy;
                    transactions[m1] = vector<int>{gain, buy, sell};
                    // record loss and action if this merged transaction is discarded
                    losses.push(-gain);
                    if (actions.find(-gain) != actions.end()) actions[-gain].push_back(vector<int>{m1});
                    else actions[-gain] = vector<vector<int>>{vector<int>{m1}};
                    // find adjacent remaining transactions
                    d = m2;
                    m2 = d+1;
                    if (m2 < transactions.size() && transactions[m2][0] == 0) m2 = transactions[m2][2];
                    transactions[d] = vector<int>{0, m1, m2}; // store adjacent remaining transactions for later use
                    t--; // keep track of number of transactions
                }
            }
            // record loss and action of if previously found adjacent remaining transactions are merged
            if (m1 >= 0 && m2 < transactions.size()) {
                int loss = -(transactions[m1][2]-transactions[m2][1]);
                losses.push(loss);
                if (actions.find(loss) != actions.end()) actions[loss].push_back(vector<int>{m1,m2});
                else actions[loss] = vector<vector<int>>{vector<int>{m1,m2}};
                // store adjacent remaining transactions for later use
                transactions[m1+1] = vector<int>{0, m1, m2};
                transactions[m2-1] = vector<int>{0, m1, m2};
            }
        }
        
        // sum total profit with required number of transactions
        int max_profit = 0;
        for (i=0; i<transactions.size(); i++) max_profit += transactions[i][0];
        return max_profit;
    }
};