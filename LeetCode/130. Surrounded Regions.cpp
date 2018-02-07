class Solution {
private:
    // for debug
    void dump(vector<vector<char>>& board) {
        int M = board.size();
        int N = board[0].size();
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                cout << board[i][j];
            }
            cout << endl;
        }
        cout << endl;
    }
    
    int flip(vector<vector<char>>& board, int i, int j) {
        int M = board.size();
        int N = board[0].size();
        if (i < 0 || i >= M || j < 0 || j >=N) return 0;
        if (board[i][j] == 'D') return 0; // definitely not captured
        if (board[i][j] == 'X') return 1;
        
        if (board[i][j] == 'O') board[i][j] = '+';  // flip to pending
        else return 2; // know it's pending
        
        // short-circuit to definitely not captured when any neighbor isn't
        int a = flip(board, i, j-1);
        if (a == 0) {
            board[i][j] = 'D';
            return 0;
        }
        int b = flip(board, i, j+1);
        if (b == 0) {
            board[i][j] = 'D';
            return 0;
        }
        int c = flip(board, i-1, j);
        if (c == 0) {
            board[i][j] = 'D';
            return 0;
        }
        int d = flip(board, i+1, j);
        if (d == 0) {
            board[i][j] = 'D';
            return 0;
        }
        
        if (a==1 && b==1 && c==1 && d==1) { // definitely captured
            board[i][j] = 'X';
            return 1;
        }
        else { // pending
            board[i][j] = '+';
            return 2;
        }
    }

public:
    void solve(vector<vector<char>>& board) {
        int M = board.size();
        if (M < 2) return;
        int N = board[0].size();
        if (N < 2) return;
        
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (board[i][j] == 'X') continue;
                int f = flip(board, i, j);
                if (f == 0) { // any remaining pending block is definitely not captured
                    for (int m = 0; m < M; m++) {
                        for (int n = 0; n < N; n++) {
                            if (board[m][n] == '+') board[m][n] = 'D';
                        }
                    }
                }
                if (f == 2) { // a whole connected block are all pending means they're all captured
                    for (int m = 0; m < M; m++) {
                        for (int n = 0; n < N; n++) {
                            if (board[m][n] == '+') board[m][n] = 'X';
                        }
                    }
                }
            }
        }
        
        // remove definitely not captured flags
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (board[i][j] == 'D') board[i][j] = 'O';
            }
        }
    }
};