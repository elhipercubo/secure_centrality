#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iostream>
#include <map>
#include <thread>
#include <vector>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
using namespace std;

//#define DEBUG

#ifndef DEBUG
    #define BOOST_DISABLE_ASSERTS 
#endif

#include <boost/multi_array.hpp>
#include <NTL/mat_ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include "defs.h"


int num_hosts;
vector<Graph> graphs;               // Input graphs for the hosts
Graph A, B;                         // OR of all graphs, sum of all graphs

vector<double> weights;             // Weights for each distance (beta)
int D, n;                           // Number of iterations, total number of vertices

void initialize() {
    setbuf(stdout, NULL);           // no buffering for stdout

    // Initialize random generator
#ifndef DEBUG
    srand(time(NULL));
#else
    srand(1000);
#endif

    weights.resize(D + 1);
    for (int i = 1; i <= D; ++i)
        weights[i] = pow(0.5, i);
}

int randint(int mod) { return rand() % mod; }

void elapsed_time() {
    static bool first = true;
    static std::chrono::steady_clock::time_point last = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    if (!first)
        std::cout << "Elapsed time = " << (double)std::chrono::duration_cast<std::chrono::microseconds>(now - last).count() /
            1e6 << " seconds\n" << std::endl;
    first = false;
    last = now;
}

// assumes boolean graph
vector<double> inf1_vector(const Graph& G) {
    vector<double> inf(n);

    puts("Hola");
//    Matrix<bool> P(n, vector<bool>(n)), Q(n, vector<bool>(n));
    bool* P = (bool *)calloc(n * n, 1);
    bool* Q = (bool *)malloc(n * n);
    for (int i = 0; i < n; ++i)
        P[i * n + i] = true;
//        P[i][i] = true;

    for (int d = 1; d <= D; ++d) {
        printf("d = %i\n", d);
        // P = P * G
        for (auto p: G) {
            int j = p.first.first, k = p.first.second;
            if (d == 1) 
                Q[j * n + k] = 1;
//                Q[j][k] = 1;
            else {
                for (int i = 0; i < n; ++i)
                    if (P[i * n + j])
//                    if (P[i][j])
//                        Q[i][k] = 1;
                        Q[i * n + k] = 1;
            }
        }
        copy(&Q[0], &Q[n * n], &P[0]);
//        copy(&Q[0][0], &Q[n][0], &P[0][0]);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (P[i * n + j]) 
//                if (P[i][j])
                    inf[i] += weights[d];
    }

    return inf;
}


vector<double> inf1_vector_small_space(vector<unordered_set<int>>& adjA) {
    vector<double> inf(n);
    vector<double> cum_weights(D + 1);
    cum_weights[D] = weights[D];
    for (int i = D; i > 0; --i) 
        cum_weights[i] = cum_weights[i + 1] + weights[i];

    const int _infinity = 1000000000;
    vector<int> dist(n), queue(n);
    int each = max(1, n / 200);
    for (int i = 0; i < n; ++i) {
        if (i % each == 0)
            printf("  i = %i\n", i);
        
        // BFS from i
        fill(&dist[0], &dist[n], _infinity);
        dist[i] = 0;

        int start = 0, end = 0;
        queue[end++] = i;
        int edges = 0;
        while (start != end) {
            int v = queue[start++];
            if (dist[v] < D) 
                for (int w : adjA[v]) {
                    edges++;
                    if (dist[w] == _infinity) {
                        dist[w] = dist[v] + 1;
                        queue[end++] = w;
                    }
                }
        }

        for (int j = 0; j < n; ++j)
            if (dist[j] < _infinity)
                inf[i] += cum_weights[dist[j]];
    }
    return inf;
}

double reconstruct_result(const vector<int>& congr, const vector<int>& modulus) {
    NTL::ZZ a, p;
    a = 0;
    p = 1;
    for (int i = 0; i < congr.size(); ++i)
        NTL::CRT(a, p, congr[i], modulus[i]);
    if (a < 0) a += p;
    double x;
#ifdef DEBUG
    assert(isfinite(x));
#endif
    NTL::conv(x, a);
    return x;
}


vector<int> find_primes(int min_value, int how_many) {
    vector<int> ret;
    NTL::PrimeSeq s;
    s.reset(min_value);
    cout << "primes:";
    while (how_many-- > 0) {
        int p = s.next();
        cout << " " << p;
        ret.push_back(p);
    }
    cout << endl;
    return ret;
}

vector<double> inf23_vector(const Graph& G) {
    vector<double> inf(n);

    // We use doubles here to avoid overflow
    vector<double> v(n), w(n);                 // vector we disclose
    for (int i = 0; i < n; ++i) v[i] = 1;

    for (int d = 1; d <= D; ++d) {
        fill(w.begin(), w.end(), 0);
        for (auto x: G) {
            int i = x.first.first, j = x.first.second;
            int A_ij = x.second;
            w[i] += A_ij * v[j];
        }
//        cout << "w =" << endl; for (int i = 0; i < n; ++i) cout << w[i] << " "; cout << endl;

        v.swap(w);
        for (int i = 0; i < n; ++i)
            inf[i] += v[i] * weights[d];
    }

    return inf;
}

vector<double> exact_inf23_vector(const Graph& G, int prime_size) {
    vector<int> mod = find_primes(prime_size, D);
    int num_mod = mod.size();

    vector<double> inf(n);
//    vector<NTL::xdouble> inf(n);

    // We use doubles here to avoid overflow
    vector<vector<int>> v(num_mod, vector<int>(n));
    for (int md = 0; md < num_mod; ++md)
        for (int i = 0; i < n; ++i)
            v[md][i] = 1;

    for (int d = 1; d <= D; ++d) {
        printf("d = %i\n", d);
        vector<vector<int>> w(num_mod, vector<int>(n));
        for (auto x: G) {
            int i = x.first.first, j = x.first.second;
            int A_ij = x.second;
            if (A_ij != 0) 
                for (int md = 0; md < num_mod; ++md) {
                    w[md][i] = (w[md][i] + A_ij * v[md][j]) % mod[md];
                }
        }

        v = w;
        for (int i = 0; i < n; ++i) {
            vector<int> congr;
            for (int md = 0; md < num_mod; ++md) 
                congr.push_back(v[md][i]);
            double value_i = reconstruct_result(congr, mod);
            inf[i] += value_i * weights[d];
        }
    }

    return inf;
}

vector<double> secure_inf3_vector(int prime_size) {
    vector<int> mod = find_primes(prime_size, D);
    int num_mod = mod.size();
    vector<double> inf(n);

    vector<vector<int>> v(num_mod, vector<int>(n));                         // vector we disclose
    for (int md = 0; md < num_mod; ++md)
        for (int i = 0; i < n; ++i)
            v[md][i] = 1;

    int M[num_mod][num_hosts][n][num_hosts];
    int N[num_mod][num_hosts][n];
    for (int d = 1; d <= D; ++d) {
        printf("d = %i\n", d);
        for (int h1 = 0; h1 < num_hosts; ++h1) {
            printf("  h = %i\n", h1);
            auto G = graphs[h1];
            for (int md = 0; md < num_mod; ++md) {
                vector<int> host_w(n);
                for (auto x: G) {
                    int i = x.first.first, j = x.first.second;
                    int A_ij = x.second;
                    if (A_ij != 0)
                        host_w[i] = (host_w[i] + A_ij * v[md][j]) % mod[md];
                }
    //            cout << "host_w =" << endl; for (int i = 0; i < n; ++i) cout << host_w[i] << " "; cout << endl;

                for (int i = 0; i < n; ++i) {
                    int sum = 0;
                    for (int x = 1; x < num_hosts; ++x) {
                        int h2 = (x + h1) % num_hosts;
                        // Send message from h1 to h2
                        int value = randint(mod[md]);
                        M[md][h2][i][h1] = value;
                        sum = (sum + value) % mod[md];
                    }
                    M[md][h1][i][h1] = (mod[md] + host_w[i] - sum) % mod[md];
                }
            }
        }

        // Now all participants share the sums they get
        fill(&N[0][0][0], &N[num_mod][0][0], 0);
        for (int h2 = 0; h2 < num_hosts; ++h2)
            for (int md = 0; md < num_mod; ++md)
                for (int i = 0; i < n; ++i) {
                    for (int h1 = 0; h1 < num_hosts; ++h1) 
                        N[md][h2][i] = (N[md][h2][i] + M[md][h2][i][h1]) % mod[md];
                }

        for (int md = 0; md < num_mod; ++md) fill(v[md].begin(), v[md].end(), 0);

        for (int h = 0; h < num_hosts; ++h)
            for (int md = 0; md < num_mod; ++md)
                for (int i = 0; i < n; ++i) {
                    // Send message N[h][i] to all
                    v[md][i] = (v[md][i] + N[md][h][i]) % mod[md];
                }

        for (int i = 0; i < n; ++i) {
            vector<int> congr;
            for (int md = 0; md < num_mod; ++md) {
                congr.push_back(v[md][i]);
//                printf("d = %i; congr_%i[%i] = %i\n", d, i, md, v[md][i]);
            }
            double value_i = reconstruct_result(congr, mod);
            inf[i] += value_i * weights[d];
        }
    }

    return inf;
}

vector<int> get_ranking(const vector<double> v) {
    vector<pair<double, int>> w;
    vector<int> ret;
    for (int i = 0; i < v.size(); ++i)
        w.push_back(make_pair(-v[i], i));
    sort(w.begin(), w.end());
    for (int i = 0; i < w.size(); ++i)
        ret.push_back(w[i].second);
    return ret;
}

using namespace NTL;

vector<ZZ_p> sum_secrets(const vector<ZZ_p>& shares1, const vector<ZZ_p>& shares2) {
    vector<ZZ_p> ret;
    for (int i = 0; i < num_hosts; ++i)
        ret.push_back(shares1[i] + shares2[i]);
    return ret;
}

ZZ_pX random_poly(int degree, ZZ_p eval0) {
    ZZ_pX poly = random_ZZ_pX(degree + 1);
    SetCoeff(poly, 0, eval0);
    return poly;
}

bool check_secret(const vector<ZZ_p>& shares, int t) {
    if (shares.size() <= t) return false;
    static vec_ZZ_p v, w;
    v.SetLength(t + 1);
    for (int i = 0; i <= t; ++i)
        v[i] = i + 1;
    w.SetLength(t + 1);
    for (int i = 0; i <= t; ++i)
        w[i] = shares[i];
    ZZ_pX poly = interpolate(v, w);
    for (int i = 0; i < shares.size(); ++i)
        if (eval(poly, to_ZZ_p(i + 1)) != shares[i]) return false;
    return true;
}

vector<ZZ_p> make_secret(ZZ_p value, int t) {
    vector<ZZ_p> ret(num_hosts);
    ZZ_pX poly = random_poly(t, value);
//    cout << "degree = " << deg(poly) << endl; for (int i = 0; i <= deg(poly); ++i) cout << coeff(poly, i) << " "; cout << endl;
    for (int i = 0; i < num_hosts; ++i) {
        ret[i] = eval(poly, to_ZZ_p(i + 1));
//        cout << i << " " << ret[i] << endl;
    }
#ifdef DEBUG
    assert(check_secret(ret, t));
#endif
    return ret;
}

ZZ_p recover_secret(const vector<ZZ_p>& shares, int t) {
    static vec_ZZ_p v, w;
    if (v.length() == 0) {
        v.SetLength(t + 1);
        for (int i = 0; i <= t; ++i)
            v[i] = i + 1;
    }
    w.SetLength(t + 1);
    for (int i = 0; i < t + 1; ++i)
        w[i] = shares[i];
    ZZ_pX poly = interpolate(v, w);
    ZZ_p x;
    x = 0;
    return eval(poly, x);
}

vector<ZZ_p> mult_secrets(const vector<ZZ_p>& shares1, const vector<ZZ_p>& shares2, const vector<ZZ_p>& coef, int t) {
#ifdef DEBUG
    assert(2 * t + 1 <= num_hosts);
    assert(check_secret(shares1, t));
    assert(check_secret(shares2, t));
#endif

    vector<ZZ_p> ret(num_hosts, to_ZZ_p(0));
    ZZ_p g[num_hosts][num_hosts];
    for (int i = 0; i < num_hosts; ++i) {
        ZZ_pX poly = random_poly(t, shares1[i] * shares2[i]);
//        vector<ZZ_p> zorra;
        for (int j = 0; j < num_hosts; ++j) {
            g[i][j] = eval(poly, to_ZZ_p(j + 1));
//            zorra.push_back(g[i][j]);
        }
    }
    for (int i = 0; i < 2 * t + 1; ++i)
        for (int j = 0; j < num_hosts; ++j)
            ret[j] += coef[i + 1] * g[i][j];
//            ret[i] += coef[j + 1] * g[i][j];

#ifdef DEBUG
    if (!check_secret(ret, t)) {
        for (int i = 0; i < ret.size(); ++i) cout << shares1[i] << " "; cout << endl;
        for (int i = 0; i < ret.size(); ++i) cout << shares2[i] << " "; cout << endl;
        for (int i = 0; i < ret.size(); ++i) cout << ret[i] << " "; cout << endl;
    }
    assert(check_secret(ret, t));
    assert(recover_secret(ret, t) == recover_secret(shares1, t) * recover_secret(shares2, t));
#endif
    return ret;
}

vector<ZZ_p> expon_secret(const vector<ZZ_p>& shares, const ZZ& n, const vector<ZZ_p>& coef, int t) {
    if (n == 0) return vector<ZZ_p>(shares.size(), to_ZZ_p(1));
    if (n == 1) return shares;
    vector<ZZ_p> r = expon_secret(shares, n / 2, coef, t);
    r = mult_secrets(r, r, coef, t);
    if (IsOdd(n)) r = mult_secrets(r, shares, coef, t);
    return r;
}

vector<ZZ_p> find_coef(int t) {
    mat_ZZ_p A, InvA;
    A.SetDims(2 * t + 1, 2 * t + 1);

    ZZ_p det;
    for (int i = 1; i <= 2 * t + 1; ++i)
        for (int j = 1; j <= 2 * t + 1; ++j) {
            ZZ_p r;
            r = i;
            A[i - 1][j - 1] = power(r, j - 1);
        }
    inv(det, InvA, A);
//    cout << A << endl; cout << InvA << endl;

    vector<ZZ_p> coef(2 * t + 2);
    for (int j = 1; j <= 2 * t + 1; ++j)
        coef[j] = InvA[0][j - 1];
    return coef;
}

vector<double> secure_inf1_vector(int prime_size) {
    assert(num_hosts >= 3);
    int t = (num_hosts - 1) / 2;

    // Find a big enough prime
    ZZ bigprime;
    bigprime = find_primes(prime_size, 1)[0];
    bigprime = power(bigprime, D);
    cout << bigprime << endl;
    ZZ_p::init(bigprime);

    // Find coefficients for multiplication method
    vector<ZZ_p> coef = find_coef(t);

    // Start process
    vector<ZZ_p> A[n][n], P[n][n], Q[n][n];
    vector<double> inf(n);
    for (int d = 1; d <= D; ++d) {
        printf("d = %i\n", d);
        if (d == 1) {
            // Share union of matrices
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    A[i][j].assign(num_hosts, to_ZZ_p(1));
                    P[i][j].resize(num_hosts);
                }

            for (int h = 0; h < num_hosts; ++h) {
                printf("  h = %i\n", h);

                int mat[n][n];
                fill(&mat[0][0], &mat[n][0], 0);
                for (auto x: graphs[h])
                    mat[x.first.first][x.first.second] = 1;

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j) {
                        auto share_h_ij = make_secret(to_ZZ_p(1 - mat[i][j]), t);
                        A[i][j] = mult_secrets(A[i][j], share_h_ij, coef, t);
                    }
            }

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < num_hosts; ++k)
                        A[i][j][k] = 1 - A[i][j][k];
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    P[i][j] = A[i][j];
        } else {
            // Q = 0
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    Q[i][j].assign(num_hosts, to_ZZ_p(0));

            // Q = P * A
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < n; ++k)
                        Q[i][j] = sum_secrets(Q[i][j], mult_secrets(P[i][k], A[k][j], coef, t));

            // Q = P
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    P[i][j] = Q[i][j];
        }

        if (d > 1) {
            // We need to truncate P
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    P[i][j] = expon_secret(P[i][j], bigprime - 1, coef, t);
        }

        // Reconstruct vector of influences
        vector<ZZ_p> v(n);
        for (int i = 0; i < n; ++i) {
            vector<ZZ_p> share(n, to_ZZ_p(0));
            for (int j = 0; j < n; ++j)
                share = sum_secrets(share, P[i][j]);
            v[i] = recover_secret(share, t);
        }

        for (int i = 0; i < n; ++i) {
            ZZ warra = v[i].LoopHole();
            double value_i;
            conv(value_i, warra);
            inf[i] += value_i * weights[d];
        }
    }

    return inf;
}

vector<double> secure_inf2_vector(int prime_size) {
    assert(num_hosts >= 3);
    int t = (num_hosts - 1) / 2;

    // Find a big enough prime
    ZZ bigprime;
    bigprime = find_primes(prime_size, 1)[0];
    bigprime = power(bigprime, D);
    cout << bigprime << endl;
    ZZ_p::init(bigprime);

    // Find coefficients for multiplication method
    vector<ZZ_p> coef = find_coef(t);

    // Start process
    vector<ZZ_p> A[n][n];
    vector<double> inf(n);
    vector<ZZ_p> v(n, to_ZZ_p(1));
    for (int d = 1; d <= D; ++d) {
        printf("d = %i\n", d);
        if (d == 1) {
            // Share union of matrices
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    A[i][j].assign(num_hosts, to_ZZ_p(1));

            for (int h = 0; h < num_hosts; ++h) {
                printf("  h = %i\n", h);
                int mat[n][n];
                fill(&mat[0][0], &mat[n][0], 0);
                for (auto x: graphs[h])
                    mat[x.first.first][x.first.second] = 1;

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j) {
                        auto share_h_ij = make_secret(to_ZZ_p(1 - mat[i][j]), t);
                        ZZ_p Aij = recover_secret(A[i][j], t), b = recover_secret(share_h_ij, t);
                        vector<ZZ_p> old = A[i][j];
                        A[i][j] = mult_secrets(A[i][j], share_h_ij, coef, t);
                        ZZ_p C = recover_secret(A[i][j], t);
                        if (Aij *  b != C) {
                            cout << Aij << " * " << b << " = " << C << endl;
                            cout << "prime = " << bigprime;
                            check_secret(old, t);
                            check_secret(share_h_ij, t);
                            check_secret(A[i][j], t);

                            for (int k = 0; k < num_hosts; ++k) cout << old[k] << " "; cout << endl;
                            for (int k = 0; k < num_hosts; ++k) cout << share_h_ij[k] << " "; cout << endl;
                            for (int k = 0; k < num_hosts; ++k) cout << A[i][j][k] << " "; cout << endl;
                        }
                        assert(Aij * b == C);
                    }
            }

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < num_hosts; ++k)
                        A[i][j][k] = 1 - A[i][j][k];

            /*
            int edges = 0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    if (recover_secret(A[i][j], t) == 1)
                        ++edges;
            printf("edges = %i\n", edges);
            exit(1);
            */
        }

        // Reconstruct vector of influences
        vector<vector<ZZ_p>> w(n, vector<ZZ_p>(num_hosts, to_ZZ_p(0)));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < num_hosts; ++k)
                    w[i][k] = w[i][k] + A[i][j][k] * v[j];

        for (int i = 0; i < n; ++i) {
            v[i] = recover_secret(w[i], t);
            cout << i << ":" << v[i] << endl;
        }

        for (int i = 0; i < n; ++i) {
            ZZ warra = v[i].LoopHole();
            double value_i;
            conv(value_i, warra);
            inf[i] += value_i * weights[d];
        }
    }

    return inf;
}

void warra(ZZ mod) {
    assert(num_hosts >= 3);
    int t = (num_hosts - 1) / 2;

    cout << mod << endl;
    ZZ_p::init(mod);
    mat_ZZ_p A, InvA;
    A.SetDims(2 * t + 1, 2 * t + 1);

    ZZ_p det;
    for (int i = 1; i <= 2 * t + 1; ++i)
        for (int j = 1; j <= 2 * t + 1; ++j) {
            ZZ_p r;
            r = i;
            A[i - 1][j - 1] = power(r, j - 1);
        }
    inv(det, InvA, A);
    cout << A << endl;
    cout << InvA << endl;

    vector<ZZ_p> coef(2 * t + 2);
    for (int j = 1; j <= 2 * t + 1; ++j)
        coef[j] = InvA[0][j - 1];

    ZZ_p secret[2];
    secret[0] = 47;
    secret[1] = 84;
    vector<ZZ_p> shares[2];

    for (int i = 0; i < 2; ++i) {
        shares[i] = make_secret(secret[i], t);
        cout << recover_secret(shares[i], t) << endl;
    }
    auto x = sum_secrets(shares[0], shares[1]);
    cout << recover_secret(x, t) << endl;
    x = mult_secrets(shares[0], shares[1], coef, t);
    cout << recover_secret(x, t) << endl;
    return;
}

void compare(const vector<double>& inf, const vector<double>& test) {
    printf("\n");
    double max_ratio = 1.0;
    for (int i = 0; i < min(15, n); ++i) {
        double r = inf[i] / test[i];
        if (fabs(inf[i] - test[i]) > 1e-6)
            max_ratio = max(max_ratio, max(r, 1.0 / r));
        printf("%i: %lf %lf ratio = %lf\n", i, inf[i], test[i], r);
    }
    auto r = get_ranking(test);
        for (int k = 0; k < min(10, n); ++k)
            printf("  v%i: rank=%i\tinfl=%lf\n", k, r[k], test[r[k]]);
    printf("max_ratio = %lf\n\n", max_ratio);
}

// copy graphs.resize(argc - 1)
int main(int argc, char* argv[]) {
    const int prime_size = 100000000;
//    num_hosts = 5;ZZ_p::init(to_ZZ(1039)); make_secret(to_ZZ_p(43), 2); return -1;
    for (int i = 0; i < argc; ++i)
        printf("%s ", argv[i]);
    puts("\n");

    int opt;
    bool duplicate_edges = false;
    bool comparison = false;
    bool all = false;
    while ((opt = getopt(argc, argv, "uca")) != -1) {
        switch (opt) {
        case 'u':
            duplicate_edges = true;
            printf("Duplicating edges (undirected graph).\n");
            break;
        case 'c':
            comparison = true;
            break;
        case 'a':
            all = true;
            break;
        default:
            fprintf(stderr, "Wrong option: %c.\n", opt);
            exit(1);
        }
    }

    if (argc - optind < 2) {
        fprintf(stderr, "Usage: %s D EDGE_FILE_1 ... EDGE_FILE_N\n", argv[0]);
        exit(-1);
    }

    D = atoi(argv[optind++]);
    printf("D = %i\n", D);
    graphs.resize(num_hosts = argc - optind);
    printf("num_hosts = %i\n", num_hosts);
    printf("Reading graphs...\n");

    printf("optind = %i\n", optind);
    elapsed_time();
    for (int i = optind; i < argc; ++i)
//        read_graph(argv[i], graphs[i - 2]);
        read_graph(argv[i], graphs[i - optind], duplicate_edges);
    n = vertices.size();
    assert(num_hosts * n * 10 < prime_size);
    int edgesA = 0, edgesB = 0;
    elapsed_time();

    printf("Computing A and  B...\n");
    vector<unordered_map<int, int>> adjB(n);
    vector<unordered_set<int>> adjA(n);
    for (int i = 0; i < num_hosts; ++i) {
        printf("h = %i\n", i);
        for (auto x: graphs[i]) {
            adjB[x.first.first][x.first.second] += x.second;
            adjA[x.first.first].insert(x.first.second);
//            B[x.first] += x.second;
            edgesB += x.second;
        }
    }
    for (int i = 0; i < n; ++i)
        for (auto x : adjB[i])
            B.push_back(make_pair(make_pair(i, x.first), x.second));
    adjB.clear();
    for (auto x: B) {
        A.push_back(make_pair(make_pair(x.first.first, x.first.second), 1));
        ++edgesA;
    }
    printf("total number of vertices = %i\n", n);
    printf("edgesA = %i\nedgesB = %i\n", edgesA, edgesB);
    initialize();
    elapsed_time();

    printf("Computing secure inf2...\n");
    elapsed_time();

    vector<double> inf[4];

    printf("Computing inf2...\n");
//    inf[2] = inf23_vector(A);
//    inf[2] = secure_inf2_vector(prime_size);
    inf[2] = exact_inf23_vector(A, prime_size);

    printf("Computing inf3...\n");
//    inf[3] = inf23_vector(B);
    inf[3] = exact_inf23_vector(B, prime_size);
    printf("Inf3 computed\n");
    elapsed_time();

    if (comparison) {
        printf("Computing inf1...\n");
        inf[1] = inf1_vector_small_space(adjA);
        elapsed_time();
        adjA.clear();

        vector<int> ranking[4];
        for (int i = 1; i <= 3; ++i) {
            printf("inf_%i:\n", i);
            ranking[i] = get_ranking(inf[i]);
            for (int k = 0; k < min(10, n); ++k)
                printf("  v%i: rank=%i\tinfl=%lf\n", k, ranking[i][k], inf[i][ranking[i][k]]);
            puts("");
        }

        for (int i = 1; i <= 3; ++i) {
            fprintf(stderr, "inf%i = [", i);
            for (int j = 0; j < inf[i].size(); ++j) {
                if (j > 0) fprintf(stderr, ",");
                fprintf(stderr, " %lf", inf[i][j]);
            }
            fprintf(stderr, " ]\n");
        }
    } else {
        vector<double> test;
        printf("Computing inf3 securely...\n");
        test = secure_inf3_vector(prime_size);
        printf("Secure Inf3 computed\n");
        elapsed_time();
        compare(inf[3], test);

        if (all)  {
            printf("Computing inf2 securely...\n");
            test = secure_inf2_vector(prime_size);
            printf("Secure Inf2 computed\n");
            elapsed_time();
            compare(inf[2], test);

            printf("Computing inf1...\n");
            inf[1] = inf1_vector_small_space(adjA);
            elapsed_time();
            test = secure_inf1_vector(prime_size);
            printf("Secure Inf1 computed\n");
            elapsed_time();
            compare(inf[1], test);
        }
    }

    return 0;
}
