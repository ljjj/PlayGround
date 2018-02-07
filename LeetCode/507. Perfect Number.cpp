class Solution {
private:
    vector<int> primes{2};
    bool isNextPrime(int n) {
        for (int i = 0; i < primes.size(); i++) {
            if (n % primes[i] == 0) return false;
            if (primes[i] > sqrt(n)) break;
        }
        return true;
    }
public:
    bool checkPerfectNumber(int num) {
        if (num <= 5) return false;
        vector<int> prime_factors{};
        vector<int> prime_factor_powers{};
        int n = num;
        while (primes.back() <= sqrt(n)) {
            if (num % primes.back() == 0) {
                prime_factors.push_back(primes.back());
                int p = 0;
                while (n % primes.back() == 0) {
                    n /= primes.back();
                    p++;
                }
                prime_factor_powers.push_back(p);
            }
            
            int next = primes.back() + ((primes.back()==2)?1:2);
            while (!isNextPrime(next)) next += 2;
            primes.push_back(next);
        }
        if (n > 1) {
            prime_factors.push_back(n);
            prime_factor_powers.push_back(1);
        }
        
        vector<int> factors{};
        int N = prime_factors.size();
        vector<int> prime_factors_used(N,0);
        while (1) {
            int f = 1;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < prime_factors_used[i]; j++) {
                    f *= prime_factors[i];
                }
            }
            factors.push_back(f);
            
            prime_factors_used[N-1]++;
            int k = N-1;
            while (k > 0 && prime_factors_used[k] > prime_factor_powers[k]) {
                prime_factors_used[k] = 0;
                k--;
                prime_factors_used[k]++;
            }
            if (k == 0 && prime_factors_used[k] > prime_factor_powers[k]) break;
        }
        
        int sum = 0;
        for (int i=0; i<factors.size(); i++) {
            sum += factors[i];
        }
        
        if (sum == 2*num) return true;
        return false;
    }
};