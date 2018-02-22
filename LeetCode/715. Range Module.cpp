class RangeModule {
private:
    vector<vector<int>> range;
public:
    RangeModule() {
        range.clear();
    }
    
    void addRange(int left, int right) {
        int i = 0;
        while (i < range.size() && left > range[i][0]) i++; // find where to insert or merge with left
        if (i == 0 || left > range[i-1][1]) { // insert because range[i-1][1] < left <= range[i][0]
            if (i == range.size() || right < range[i][0]) range.insert(range.begin()+i, vector<int>{left, right}); // pure insert
            else { // need to merge with right after insert
                int j = i;
                while (j < range.size() && right > range[j][1]) j++; // find where to merge with right
                if (j == range.size() || right < range[j][0]) { // merge with right not including j
                    range[i] = vector<int>{left, right}; // range[j-1][1] < right < range[j][0]
                    range.erase(range.begin()+i+1, range.begin()+j);
                }
                else { // merge with right including j
                    range[i] = vector<int>{left, range[j][1]}; // range[j][0] <= right <= range[j][1]
                    range.erase(range.begin()+i+1, range.begin()+j+1);
                }
            }
        }
        else { // merge with left because range[i-1][0] < left <= range[i-1][1]
            // no need to merge with right
            if (i == range.size() || right < range[i][0]) range[i-1][1] = max(right, range[i-1][1]);
            else { // need to merge with right, similar to above
                int j = i;
                while (j < range.size() && right > range[j][1]) j++;
                if (j == range.size() || right < range[j][0]) {
                    range[i-1][1] = right; // range[j-1][1] < right < range[j][0]
                    range.erase(range.begin()+i, range.begin()+j);
                }
                else {
                    range[i-1][1] = range[j][1]; // range[j][0] <= right <= range[j][1]
                    range.erase(range.begin()+i, range.begin()+j+1);
                }
            }
        }
    }
    
    bool queryRange(int left, int right) {
        int i = 0;
        while (i < range.size() && left >= range[i][0]) i++; // the query range could be in range[i-1]
        if (i == 0 || left > range[i-1][1]) return false;
        else {
            if (right <= range[i-1][1]) return true;
            else return false;
        }
    }
    
    void removeRange(int left, int right) {
        int i = 0;
        while (i < range.size() && left > range[i][0]) i++; // find where to start removal from left
        if (i == 0 || left >= range[i-1][1]) { // remove from i because range[i-1][1] < left <= range[i][0]
            if (i < range.size() && right > range[i][0]) {
                int j = i;
                while (j < range.size() && right >= range[j][1]) j++; // find where to remove on the right
                if (j < range.size() && right > range[j][0]) range[j][0] = right; // remove part of j because range[j][0] < right < range[j][1]
                range.erase(range.begin()+i, range.begin()+j);
            }
        }
        else { // remove part of i-1 because range[i-1][0] < left < range[i-1][1]
            if (i == range.size() || right <= range[i][0]) {
                if (right < range[i-1][1]) range.insert(range.begin()+i, vector<int>{right, range[i-1][1]}); // i-1 is broken into two because range[i-1][0] < left < right < range[i-1][1]
            }
            else {
                int j = i;
                if (left == 10) cout << right << range[j][0] << range[j][1] << endl;
                while (j < range.size() && right >= range[j][1]) j++; // find where to remove on the right
                if (j < range.size() && right > range[j][0]) range[j][0] = right; // remove part of j because range[j][0] < right < range[j][1]
                range.erase(range.begin()+i, range.begin()+j);
            }
            range[i-1][1] = left; // remove part of i-1
        }
    }
};

/**
 * Your RangeModule object will be instantiated && called as such:
 * RangeModule obj = new RangeModule();
 * obj.addRange(left,right);
 * bool param_2 = obj.queryRange(left,right);
 * obj.removeRange(left,right);
 */