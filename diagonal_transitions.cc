#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std;

namespace EditDistance {

constexpr int MAX_LOG_N = 20;
constexpr int MAX_N = 1 << MAX_LOG_N;

int n;
int len_pattern, len_text;
int suffix_array_pos[MAX_LOG_N][MAX_N];
int suffix_array[MAX_N];

// lcp_maxes[k][i] = max (lcp[i], lcp[i+1}, ... lcp[i-1+2^^k])
//      where lcp[i] = Longest common prefix of the suffixes i and i+1
int lcp_mins[MAX_LOG_N][MAX_N];


int max(int a, int b) {
    return (a>b)?a:b;
}
int max(int a, int b, int c) {
    return max(a, max(b, c));
}

int FloorLog(int x) {
    int res = 0;
    while ((1 << res) <= x) res++;
    res--;
    return res;
}

int CeilLog(int x) {
    int res = 0;
    while ((1<<res) < x) res++;
    return res;
}


// The complexity of this function is O(NlgN) where
// N is the length of the concatenation pattern + text
void ComputeSuffixArray(const string& p, const string& t) {
    len_pattern = p.length();
    len_text = t.length();

    const string str = p + t;
    cout << str << endl;
    n = str.length();

    for (int i = 0; i < n; i++) {
        suffix_array_pos[0][i] = str[i] - 'a';
        cout << suffix_array_pos[0][i] << " ";
    }
    cout << endl;

    vector< tuple<int, int, int> > l (n);
    const int log_n = CeilLog(n);
    for (int k = 1, len = 1; k <= log_n; k++, len *= 2) {
        for (int i = 0; i < n; i++) {
            l[i] = make_tuple(suffix_array_pos[k-1][i], suffix_array_pos[k-1][i+len], i);
        }
        sort(l.begin(), l.end());
        for (int i = 0; i < n; i++) {
            const bool equal_to_last = i > 0 &&
                get<0>(l[i-1]) == get<0>(l[i]) &&
                get<1>(l[i-1]) == get<1>(l[i]);
            const int pos = get<2>(l[i]);
            if (equal_to_last) {
                const int pos_last = get<2>(l[i-1]);
                suffix_array_pos[k][pos] = suffix_array_pos[k][pos_last];
            } else {
                suffix_array_pos[k][pos] = i;
            }
        }
        for (int i = 0; i < n; i++) {
            cout << suffix_array_pos[k][i] << " ";
        }
        cout << endl;
    }

    for (int i = 0; i < n; i++) {
        suffix_array[suffix_array_pos[log_n][i]] = i;
    }
}

// pi and ti are indexes in the concatenation pattern + text
// The complexity of this function is O(lgN) where
// N is the length of the concatenation pattern + text.
int ComputeLongestCommonPrefixInConcat (int pi, int ti) {
    int result = 0;
    const int log_n = CeilLog(n);
    for (int k = log_n; k >= 0 && pi < n && ti < n; k--) {
        if (suffix_array_pos[k][pi] == suffix_array_pos[k][ti]) {
            result += (1 << k);
            pi += (1 << k);
            ti += (1 << k);
        }
    }
    return result;
}

// pi is index in the the pattern.
// ti is index in the text.
int ComputeLongestCommonPrefix (int pi, int ti) {
    return ComputeLongestCommonPrefixInConcat(pi, ti+len_pattern);
}

// The complexity of this function is O(NlgN) where
// N is the length of the concatenation pattern + text
void PreprocessLCP() {
    const int log_n = CeilLog(n);

    /*
    cout << "\n\nLCP\n";
    for (int i = 0; i+1 < n; i++) {
        cout << i%10 << " ";
    }
    cout << endl;
    */

    for (int i = 0; i+1 < n; i++) {
        lcp_mins[0][i] =
            ComputeLongestCommonPrefixInConcat(suffix_array[i], suffix_array[i+1]);
        cout << lcp_mins[0][i] << " ";
    }
    cout << endl;

    for (int k = 1; k <= log_n; k++) {
        for (int i = 0; i+1 < n; i++) {
            const int next_half_start = i + (1<<(k-1));
            lcp_mins[k][i] = lcp_mins[k-1][i];

            if (next_half_start+1 < n) {
                lcp_mins[k][i] =
                    min(lcp_mins[k-1][i], lcp_mins[k-1][next_half_start]);
            }
            //cout << lcp_mins[k][i] << " ";
        }
        //cout << endl;
    }
}

// i and j are indexes in the concatenation pattern+text
// The complexity of this function is O(1)
int GetLongestCommonPrefixInConcat(int i, int j) {
    const int log_n = CeilLog(n);

    int suff_i = suffix_array_pos[log_n][i];
    int suff_j = suffix_array_pos[log_n][j];
    if (suff_j < suff_i) swap(suff_i, suff_j);

    const int len = suff_j-suff_i;
    const int log_len = FloorLog(len);
    /*
    cout << "len = " << len << endl;
    cout << "suff_j = " << suff_j << endl;
    cout << "suff_i = " << suff_i << endl;
    cout << "log_len  = " << log_len << endl;
    */
    const int lcp = min(lcp_mins[log_len][suff_i], lcp_mins[log_len][suff_j-(1<<log_len)]);
    return lcp;
}

// pi is index in the the pattern.
// ti is index in the text.
// The complexity of this function is O(1)
int GetLongestCommonPrefix(int pi, int ti) {
    return GetLongestCommonPrefixInConcat(pi, ti+len_pattern);
}

int ComputeEditDistance(const string& p,
                    const string& t, const int k) {

    ComputeSuffixArray(p, t);
    PreprocessLCP();

    const int m = p.length(), n = t.length();
    if (k >= p.length()) return 0;

    const int table_w = n - m + k + 3;
    const int table_h = k + 2;
    int col, d;
    vector<int> preprevious(table_h, -2);
    vector<int> previous(table_h, -1);
    vector<int> current(table_h);

    int matches = 0;

    for (int j = 2; j < table_w; j++)
    {
        current[0] = j - 2;
        for (int i = 1; i < table_h; i++)
        {
            col = max(preprevious[i - 1] + 1,
                      previous[i - 1] + 1,
                      current[i - 1]);
            d = j - i - 1;

            col += GetLongestCommonPrefix( col - d, col);
            //while (col < n && col - d < m && p[col - d] == t[col]) ++col;
            current[i] = min(col, m + d);
        }
        if (current[k + 1] == m + j - k - 2) {
                matches++; // matching substring found at position current[k + 1]
                cout << "pos: " << current[k + 1] << endl;
        }
        swap(preprevious, previous);
        swap(previous, current);
    }
    return matches;
}

} // namespace EditDistance

int main() {
    cout << "res= " << EditDistance::ComputeEditDistance("pesho", "peshxxxxpeshox", 1) << endl;
    // cout << EditDistance::GetLongestCommonPrefix(0, 9) << endl;
    // cout << EditDistance::ComputeLongestCommonPrefix(0, 9) << endl;
}
