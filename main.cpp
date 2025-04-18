#include <iostream>
#include <cmath>
#include <thread>
#include <chrono>
#include <vector>
#include <omp.h>

using namespace std;
using namespace std::chrono;


float random_int(int low, int high);
float random_data(float low, float hi);
float random_float_step(float low, float high, float step);
float pi = 3.14159;
float twosqrtpi = sqrt(pi*2.0f);

struct OptionDetails {
    vector<float> S_list;
    float S_1;
    float K_1;
    float T_1;
    float r_1;
    float v_1;
    int N = 1024;
    vector<float> S, r, v, K, T;
    int iter;

    OptionDetails(int iterations = 10, bool is_test = true) {
        iter = iterations;
        if (is_test) {
            S_list = {90, 95, 100, 105, 110};
            K_1 = 100;
            T_1 = 1;
            r_1 = 0.03;
            v_1 = 0.3;
        } else {
            S.resize(iter);
            r.resize(iter);
            v.resize(iter);
            K.resize(iter);
            T.resize(iter);
        }
    }
};

float jarrowRuddCall(float S, float K, float T, float r, float v, int N);
float jarrowRuddPut(float S, float K, float T, float r, float v, int N);
float Normal(float d);
float d1(float S, float K, float T, float r, float v);
float d2(float S, float K, float T, float r, float v);
float Normalprime(float d);
float calc_rho(float S, float K, float T, float r, float v, int type);
float calc_delta(float S, float K, float T, float r, float v, int Type);
float calc_theta(float S, float K, float T, float r, float v, int type);
float calc_vega(float S, float K, float T, float r, float v);
float calc_gamma(float S, float K, float T, float r, float v);

int main() {

    OptionDetails params(true);

    for (int i = 0; i < params.S_list.size(); i++) {
        float call_prices = jarrowRuddCall(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1, params.N);
        float put_prices = jarrowRuddPut(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1, params.N);
        float delta_call = calc_delta(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1, 0);
        float delta_put = calc_delta(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1, 1);
        float vega = calc_vega(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1);
        float theta_call = calc_theta(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1, 0);
        float theta_put = calc_theta(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1, 1);
        float rho_call = calc_rho(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1, 0);
        float rho_put = calc_rho(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1, 1);
        float gamma = calc_gamma(params.S_list[i], params.K_1, params.T_1, params.r_1, params.v_1);

        cout << "Stock Price: " << params.S_list[i] << ", Call Price: " << call_prices << ", Delta: " << delta_call << ", Gamma: " << gamma <<
            ", Vega: " << vega << ", Theta: " << theta_call << ", Rho: " << rho_call << endl;

        cout << "Stock Price: " << params.S_list[i] << ", Put Price: " << put_prices << ", Delta: " << delta_put << ", Gamma: " << gamma <<
            ", Vega: " << vega << ", Theta: " << theta_put << ", Rho: " << rho_put << endl;
    }

    // part 2
    int iter = 10000;

    OptionDetails param_timed(iter, false);

    for (int i = 0; i < iter; i++) {
        param_timed.S[i]  = random_int(20,100);
        param_timed.r[i] = random_data(0.01,0.05);
        param_timed.v[i] = random_float_step(0.1, 1, 0.1);
        param_timed.K[i] =  param_timed.S[i] + random_int(-10, 10);
        param_timed.T[i] = random_float_step(0.5, 3, 0.5);

    }


    high_resolution_clock::time_point t3;
    t3 = high_resolution_clock::now();

#pragma omp parallel for
    for (int i = 0; i < iter; i++) {
        float S = param_timed.S[i];
        float K = param_timed.K[i];
        float T = param_timed.T[i];
        float r = param_timed.r[i];
        float v = param_timed.v[i];
        jarrowRuddCall(S, K, T, r, v, param_timed.N);
        calc_delta(S, K, T, r, v, 0);
        calc_vega(S, K, T, r, v);
        calc_theta(S, K, T, r, v, 0);
        calc_rho(S, K, T, r, v, 0);
        calc_gamma(S, K, T, r, v);
        jarrowRuddPut(S, K, T, r, v, param_timed.N);
        calc_delta(S, K, T, r, v, 1);
        calc_theta(S, K, T, r, v, 1);
        calc_rho(S, K, T, r, v, 1);
    }


    high_resolution_clock::time_point t4 =
                high_resolution_clock::now();


    cout<<"Elapsed Time:" << duration_cast<milliseconds>(t4 - t3).count()<< " ms" <<endl;
}

struct Node {
    float S = 0.0f;
    float C = 0.0f;
    float P = 0.0f;
};

float jarrowRuddCall(float S, float K, float T, float r, float v, int N) {
    const float dt = T / N;
    const float nu = r - 0.5 * v * v;
    const float disc = exp(- r * dt);
    const float sqrt_dt = sqrt(dt);

    const float p = 0.5f;
    const float u = exp(nu * dt + v * sqrt_dt);
    const float d = exp(nu * dt - v * sqrt_dt);

    vector<vector<Node>> tree(N + 1);

    for (int i = 0; i <= N; i++) {
        tree[i].resize(i + 1);
    }

    tree[0][0].S = S;

    for (int i = 1; i <= N; i++) {
#pragma omp parallel for schedule(dynamic)
            for (int j = 0; j <= i; j++) {
                if (j == 0) {
                    tree[i][j].S = tree[i - 1][j].S * d;
                } else {
                    tree[i][j].S = tree[i - 1][j - 1].S * u;
                }
            }
        }

#pragma omp parallel for
        for (int j = 0; j <= N; j++) {
            tree[N][j].C = max(0.0f, tree[N][j].S - K);
        }

        for (int ir = N - 1; ir >= 0; --ir) {
#pragma omp simd
            for (int j = 0; j <= ir; j++) {
                tree[ir][j].C = disc * (p * tree[ir + 1][j + 1].C + (1 - p) * tree[ir + 1][j].C);
            }
        }
        return tree[0][0].C;
    }

float jarrowRuddPut(float S, float K, float T, float r, float v, int N) {

    const float sigmasq = v * v;
    const float dt = T / N;
    const float nu = r - 0.5 * sigmasq;
    const float disc = exp(-r * dt);
    const float sqrt_dt = sqrt(dt);

    const float p = 0.5f;
    const float u = exp(nu * dt + v * sqrt_dt);
    const float d = exp(nu * dt - v * sqrt_dt);

    vector<vector<Node>> tree(N + 1);

    for (long i = 0; i <= N; i++) {
        tree[i].resize(i + 1);
    }

    tree[0][0].S = S;

    for (int i = 1; i <= N; i++) {
        for (int j = 0; j <= i; j++) {
            if (j == 0) {
                tree[i][j].S = tree[i - 1][j].S * d;
            } else {
                tree[i][j].S = tree[i - 1][j - 1].S * u;
            }
        }
    }

    for (int j = 0; j <= N; j++) {
        tree[N][j].P = max(0.0f, K - tree[N][j].S);
    }

    for (int ir = N - 1; ir >= 0; --ir) {
        for (int j = 0; j <= ir; j++) {
            tree[ir][j].P = disc * (p * tree[ir + 1][j + 1].P + (1 - p) * tree[ir + 1][j].P);
        }
    }

    return tree[0][0].P;
}

float random_data(float low, float hi) {
    return low + ((float)rand() / (float)RAND_MAX) * (hi - low);
}

float random_int(int low, int high) {
    return low + (rand() % (high - low));
}

float random_float_step(float low, float high, float step) {
    return low + (rand() / (float)RAND_MAX) * (high - low);
}

float d1(float S, float K, float T, float r, float v) {
    return (log(S / K) + (r + 0.5f * v * v) * T) / (v * sqrt(T));
}

float d2(float S, float K, float T, float r, float v) {
    return d1(S, K, T, r, v) - v * sqrt(T);
}

float Normal(float d) {
    return 0.5f * (1.0f + erf(d / twosqrtpi));
}

float calc_gamma(float S, float K, float T, float r, float v) {
    float d1_x = d1(S, K, T, r, v);
    return ((1.0f / twosqrtpi) * exp((-pow(d1_x, 2.0f)) / 2.0f)) / (S * v * sqrt(T));
}

float calc_delta(float S, float K, float T, float r, float v, int Type) {
    float d1_x = d1(S, K, T, r, v);
    if (Type == 0) {
        return Normal(d1_x);
    }
    if (Type == 1) {
        return Normal(d1_x) - 1;
    }
}

float Normalprime(float d) {
    return (1.0f / twosqrtpi) * exp(-0.5f * d * d);
}

float calc_vega(float S, float K, float T, float r, float v) {
    float d1_x = d1(S, K, T, r, v);
    return 0.01f * S * sqrt(T) * Normalprime(d1_x);
}

float calc_theta(float S, float K, float T, float r, float v, int type) {

    float d1_x = d1(S, K, T, r, v);
    float d2_x = d2(S, K, T, r, v);
    if (type == 0) {
        return -((S * Normalprime(d1_x) * v) / (2.0f * sqrt(T)) - r * K * exp(-r * T) * Normal(d2_x)) / 365.0f;
    }
    if (type == 1) {
        return -((S * Normalprime(d1_x) * v) / (2.0f * sqrt(T)) + r * K * exp(-r * T) * Normal(d2_x)) / 365.0f;
    }
}

float calc_rho(float S, float K, float T, float r, float v, int type) {
    if (type == 0) {
        return 0.01f * K * T * exp(-r * T) * Normal(d2(S, K, T, r, v));
    }
    if (type == 1) {
        return 0.01f * K * T * exp(-r * T) * Normal(-d2(S, K, T, r, v));
    }
}