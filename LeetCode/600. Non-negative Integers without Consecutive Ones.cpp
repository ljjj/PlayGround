class Solution {
private:
    vector<vector<int>> oneFill {};
    int maxPossible(int msb) {
        int m = 1;
        for (int i = 1; i<=(msb+1)/2; i++) m += fillOnes(msb, i);
        return m;
    }
    int fillOnes(int x, int n) {
        if (n == 0) return 1;
        if (x < 0) return 0;
        if (n > (x+1)/2) return 0;
        while (oneFill.size() < x+1) oneFill.push_back(vector<int>{1, oneFill.size()});
        while (oneFill[x].size() < n+1) oneFill[x].push_back(0);
        int f = oneFill[x][n];
        if (f > 0) return f;
        int y = x;
        while (y > 0) {
            f += fillOnes(y-2, n-1);
            y--;
        }
        oneFill[x][n] = f;
        return f;
    }
public:
    int findIntegers(int num) {
        if (num == 0) return 1;
        int i, mask, fi = 0;
        for (i = 30; i>=0; i--) {
            mask = 1<<i;
            if ((mask & num) != 0) break;
        }
        fi += maxPossible(i);
        if (((mask>>1) & num) != 0) fi += maxPossible(i-1);
        else fi += findIntegers(num - mask);
        return fi;
    }
};