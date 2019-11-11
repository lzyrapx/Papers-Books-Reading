#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const int mod = 1e9 + 7;

const int N = (1 << 25) + 12345678, P = 998244353, G = 3;

template <class T>
inline void read(T &x){
    int ch = 0, f = 0; x = 0;
    for(; !isdigit(ch); ch = getchar()) if(ch == '-') f = 1;
    for(; isdigit(ch); ch = getchar()) x = x * 10 + ch - 48;
    if(f) x = -x;
}

namespace poly
{
    int rev[N], W[N], invW[N], len, lg;
    inline int Pow(int a, int b)
    {
        int ans = 1;
        for (; b; b >>= 1, a = 1ll * a * a % P)
            if (b & 1)
                ans = 1ll * ans * a % P;
        return ans;
    }
    inline void init()
    {
        for (int k = 2; k < N; k <<= 1)
            W[k] = Pow(G, (P - 1) / k), invW[k] = Pow(W[k], P - 2);
    }
    inline void timesinit(int lenth)
    {
        for (len = 1, lg = 0; len <= lenth; len <<= 1, lg++)
            ;
        for (int i = 0; i < len; i++)
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (lg - 1));
    }
    inline void DFT(int *a, int sgn)
    {
        for (int i = 0; i < len; i++)
            if (i < rev[i])
                swap(a[i], a[rev[i]]);
        for (int k = 2; k <= len; k <<= 1)
        {
            int w = ~sgn ? W[k] : invW[k];
            for (int i = 0; i < len; i += k)
            {
                int now = 1;
                for (int j = i; j < i + (k >> 1); j++)
                {
                    int x = a[j], y = 1ll * a[j + (k >> 1)] * now % P;
                    a[j] = (x + y) % P, a[j + (k >> 1)] = (x - y + P) % P;
                    now = 1ll * now * w % P;
                }
            }
        }
        if (sgn == -1)
        {
            int Inv = Pow(len, P - 2);
            for (int i = 0; i < len; i++)
                a[i] = 1ll * a[i] * Inv % P;
        }
    }
    inline void getinv(int *a, int *b, int n)
    {
        static int tmp[N];
        if (n == 1)
            return (void)(b[0] = Pow(a[0], P - 2));
        getinv(a, b, (n + 1) / 2);
        timesinit(n * 2 - 1);
        for (int i = 0; i < len; i++)
            tmp[i] = i < n ? a[i] : 0;
        DFT(tmp, 1), DFT(b, 1);
        for (int i = 0; i < len; i++)
            b[i] = 1ll * (2 - 1ll * tmp[i] * b[i] % P + P) % P * b[i] % P;
        DFT(b, -1);
        for (int i = n; i < len; i++)
            b[i] = 0;
        for (int i = 0; i < len; i++)
            tmp[i] = 0;
    }
    inline void getsqrt(int *a, int *b, int n)
    {
        static int tmp1[N], tmp2[N];
        if (n == 1)
            return (void)(b[0] = 1);
        getsqrt(a, b, (n + 1) / 2);
        for (int i = 0; i < n; i++)
            tmp1[i] = a[i];
        getinv(b, tmp2, n), timesinit(n * 2 - 1);
        DFT(tmp1, 1), DFT(tmp2, 1);
        for (int i = 0; i < len; i++)
            tmp1[i] = 1ll * tmp1[i] * tmp2[i] % P;
        DFT(tmp1, -1);
        for (int i = 0; i < len; i++)
            b[i] = 1ll * (b[i] + tmp1[i]) % P * Pow(2, P - 2) % P;
        for (int i = n; i < len; i++)
            b[i] = 0;
        for (int i = 0; i < len; i++)
            tmp1[i] = tmp2[i] = 0;
    }
    inline void getln(int *a, int *b, int n)
    {
        static int tmp[N];
        getinv(a, b, n), timesinit(n * 2 - 1);
        for (int i = 1; i < n; i++)
            tmp[i - 1] = 1ll * a[i] * i % P;
        DFT(tmp, 1), DFT(b, 1);
        for (int i = 0; i < len; i++)
            b[i] = 1ll * tmp[i] * b[i] % P;
        DFT(b, -1);
        for (int i = len - 1; i; i--)
            b[i] = 1ll * b[i - 1] * Pow(i, P - 2) % P;
        b[0] = 0;
        for (int i = n; i < len; i++)
            b[i] = 0;
        for (int i = 0; i < len; i++)
            tmp[i] = 0;
    }
    inline void getexp(int *a, int *b, int n)
    {
        static int tmp[N];
        if (n == 1)
            return (void)(b[0] = 1);
        getexp(a, b, (n + 1) / 2);
        getln(b, tmp, n), timesinit(n * 2 - 1);
        for (int i = 0; i < n; i++)
            tmp[i] = (!i - tmp[i] + a[i] + P) % P;
        DFT(tmp, 1), DFT(b, 1);
        for (int i = 0; i < len; i++)
            b[i] = 1ll * b[i] * tmp[i] % P;
        DFT(b, -1);
        for (int i = n; i < len; i++)
            b[i] = 0;
        for (int i = 0; i < len; i++)
            tmp[i] = 0;
    }
    inline void power(int *a, int *b, int n, int m, ll k)
    {
        static int tmp[N];
        for (int i = 0; i < m; i++)
            b[i] = 0;
        int fir = -1;
        for (int i = 0; i < n; i++)
            if (a[i])
            {
                fir = i;
                break;
            }
        if (fir && k >= m)
            return;
        if (fir == -1 || 1ll * fir * k >= m)
            return;
        for (int i = fir; i < n; i++)
            b[i - fir] = a[i];
        for (int i = 0; i < n - fir; i++)
            b[i] = 1ll * b[i] * Pow(a[fir], P - 2) % P;
        getln(b, tmp, m);
        for (int i = 0; i < m; i++)
            b[i] = 1ll * tmp[i] * (k % P) % P, tmp[i] = 0;
        getexp(b, tmp, m);
        for (int i = m; i >= fir * k; i--)
            b[i] = 1ll * tmp[i - fir * k] * Pow(a[fir], k % (P - 1)) % P;
        for (int i = 0; i < fir * k; i++)
            b[i] = 0;
        for (int i = 0; i < m; i++)
            tmp[i] = 0;
    }
} // namespace poly

using poly::DFT;
using poly::Pow;
using poly::timesinit;
int a[N], b[N], A[N], f[N], g[N];
int n, m;
int main()
{
    // Given a set contains m kinds of numbers, the number of a[i] is b[i], (b[i] != 0)
    // subset sum is [1, n]
    // return the number of such subsets I modulo p, p is prime
    poly::init();
    cin >> n >> m;
    for (int i = 1; i <= m; i++)
    {
        cin >> a[i] >> b[i];
        if(!b[i]) b[i] = n / a[i];
        if(a[i] <= n) A[a[i]]++;
        if(1ll * a[i] * (b[i] + 1) <= n) A[a[i]*(b[i]+1)]--;
    }

    for(int i = 0; i <= n; i++) (A[i] += P) %= P;

    for(int j = 1; j <= n; j++) {
        int Inv = Pow(j, P - 2);
        for(int i = 1; i <= n / j; i++) {
            (f[i*j] += 1ll * A[i] * Inv % P) %= P;
        }
    }
    poly::getexp(f, g, n + 1);
    for(int i = 1; i <= n; i++) {
       cout << g[i] << " ";
    }
    cout << endl;
    return 0;
}
// g++ subset-sum-in-polynomial.cpp -std=c++11 -O2 -o main
// ./main