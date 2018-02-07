class Solution {
private:
    vector<int> primes{2,3,5,7,11,13,17,19,23,29,31};
    set<int> getFactors(int n) {
        if (n<1) return set<int>{};
        if (n==1) return set<int>{1};
        int f = n;
        int d = 1;
        for (int i=0; i<primes.size(); i++) {
            if (primes[i] > sqrt(n)) break;
            if (n%primes[i] == 0) {
                d = primes[i];
                f = n/d;
                break;
            }
        }
        if (f==n) return set<int>{n,1};
        set<int> fac = getFactors(f);
        vector<int> new_fac;
        for (set<int>::iterator it=fac.begin(); it != fac.end(); it++) new_fac.push_back((*it)*d);
        for (int i=0; i<new_fac.size(); i++) fac.insert(new_fac[i]);
        return fac;
    }
public:
    int minSteps(int n) {
        set<int> fac = getFactors(n);
        vector<int> f;
        for (set<int>::iterator it=fac.begin(); it != fac.end(); it++) f.push_back(*it);
        sort(f.begin(), f.end());
        vector<int> s{0};
        for (int i=1; i<f.size(); i++) s.push_back(f[i]);
        for (int i=2; i<s.size(); i++) {
            for (int j=1; j<i; j++) {
                if (f[i] % f[j] == 0) {
                    int c = s[j] + f[i]/f[j];
                    if (c < s[i]) s[i] = c;
                }
            }
        }
        return s.back();
    }
};