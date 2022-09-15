#include <time.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <queue>
#include <stack>
#include <string>
#include <typeinfo>
#include <vector>

// To save std::
using namespace std;
// OR
using std::cout;
using std::endl;

using std::queue;
using std::stack;

//===================================================================
// macro
int power(int base, int exponent) {
    int result = 1;
    if (exponent == 0) {
        return result;
    }
    for (int i = 1; i <= exponent; i++) {
        result = result * base;
    }
    return result;
}

//==================================================================
// Basics
// tail recursion
int mod(int a, int b) {
    while (a >= b) {
        int d = a - b;
        a = d;
    }
    return a;
}
// structrual recursion
// int mod(int a, int b) {
//     if (a < b) {
//         return a;
//     } else {
//         return mod(a - b, b);
//     }
// }

// tail recursion
int gcd(int a, int b) {
    while (b != 0) {
        int r = a % b;
        a = b;
        b = r;
    }
    return a;
}
// structural recursion
// int gcd(int a, int b) {
//     if (b == 0) {
//         return a;
//     } else {
//         return gcd(b, a % b);
//     }
// }

// tail recursion
int gcdSteps(int a, int b) {
    int step = 0;
    while (b != 0) {
        int r = a % b;
        a = b;
        b = r;
        step = step + 1;
    }
    return step;
}
// structrual recursion
// int gcdSteps(int a, int b) {
//     if (b == 0) {
//         return 0;
//     } else {
//         return 1 + gcdSteps(b, a % b);
//     }
// }

int is_congruent(int a, int b, int n) {
    return (((a - b) % n) == 0);
}

// compute the residual class of a mod n
int residue(int a, int n) {
    int r = a % n;
    if (r >= 0) {
        return r;
    } else {
        return (n + r);
    }
}

pair<int, int> bezoutCoef(int m, int n) {
    int g = gcd(m, n);
    m = m / g;
    n = n / g;
    pair<int, int> v1(1, 0);
    pair<int, int> v2(0, 1);
    while (n > 1) {
        int r = m % n;
        int q = (int)m / n;
        pair<int, int> v3 = make_pair(v1.first - q * v2.first, v1.second - q * v2.second);
        v1 = v2;
        v2 = v3;
        m = n;
        n = r;
    }
    return v2;
}

// find inverse of a modulo n
int inverse(int a, int n) {
    if (gcd(a, n) != 1) {
        cout << "No inverse" << endl;
        return -1;
    } else {
        pair<int, int> coeff = bezoutCoef(a, n);
        int inverse = coeff.first;
        if (inverse < 0) {
            inverse = n + inverse;
        }
        return inverse;
    }
}

// n must be positive integer
// find highest power of base b in n
int highestPower(int n, int b) {
    int degree = 0;
    while (n >= (power(b, degree))) {
        degree = degree + 1;
    }
    return (degree - 1);
}

// comupte the number of digits of an integer
// digitNum = function(b) {
//  return (floor(log10(b) + 1))
// }

int digitNum(int a) {
    return (highestPower(a, 10) + 1);
}

int highestPowerDivisible(int n, int b) {
    int degree = 0;
    while (n % power(b, degree) == 0) {
        degree = degree + 1;
    }
    return (degree - 1);
}

// n must be postive integer
vector<int> binaryPowerRep(int n) {
    vector<int> powerList;
    while (n > 0) {
        int p = highestPower(n, 2);
        powerList.push_back(p);
        n = n - power(2, p);
    }
    return powerList;
}

// compute a^e modulo n
int modExp(int a, int e, int n) {
    int result = 1;
    if (e == 0) {
        return result;
    } else {
        vector<int> powerList = binaryPowerRep(e);
        // ascending order
        sort(powerList.begin(), powerList.end());

        int startPwr = 0;
        int temp = 1;
        for (int i = 0; i < powerList.size(); i++) {
            int pwr = powerList[i];
            for (int j = startPwr; j <= pwr; j++) {
                if (j == 0) {
                    temp = a;
                } else {
                    temp = power(temp, 2) % n;
                }
            }
            startPwr = pwr + 1;
            result = (result * temp) % n;
        }
        return (result);
    }
}

vector<int> discreteLog(int b, int r, int n) {
    vector<int> logs;
    for (int i = 1; i <= 50; i++) {
        if ((modExp(b, i, n) - r) % n == 0) {
            logs.push_back(i);
        }
    }
    return logs;
}

// reduced residue system
vector<int> coprimeClass(int n) {
    vector<int> classes;
    for (int i = 0; i <= n - 1; i++) {
        if (gcd(i, n) == 1) {
            classes.push_back(i);
        }
    }
    return classes;
}

// compute phi(n) number
int eulerPhi(int n) {
    return coprimeClass(n).size();
}

// compute order of an integer modulo n
// ord=0 => ord=phi(n)
int ord(int a, int n) {
    int i = 1;
    while (modExp(a, i, n) != 1) {
        i = i + 1;
    }
    return i;
}

vector<int> primitiveRoot(int n) {
    vector<int> classes = coprimeClass(n);
    vector<int> roots;
    int phi = eulerPhi(n);
    for (int i = 0; i < classes.size(); i++) {
        int a = classes[i];
        int d = ord(a, n);
        if (d == phi) {
            roots.push_back(a);
        }
    }
    return roots;
}

// compute the index i such that r^i = a (mod n) where a is from integers with gcd(a, n) = 1
int ind(int r, int a, int n) {
    int i = 0;
    while (modExp(r, i, n) != a) {
        i = i + 1;
    }
    return i;
}

// display the index table for every a in {a | gcd(a, n) = 1}
vector<int> indexTable(int r, int n) {
    vector<int> indices;
    vector<int> reducedResidues = coprimeClass(n);
    for (int a : reducedResidues) {
        int i = ind(r, a, n);
        indices.push_back(i);
    }
    return indices;
}

// determine given integer is a perfect square
// Note: To determine a double is an integer, we can use:
// if Math.round(n) == n;
// if (n % 1 == 0); i.e (312.15 % 1 = 0.15)
bool isSquare(int n) {
    if (n < 1) {
        return false;
    }
    double root = sqrt(n);
    return root == (int)root;
}

// Fermat's factorization, produce a pair of factors in array
pair<int, int> fermatFactor(int n) {
    double root = sqrt(n);
    // Math.ceil() returns double even though it is an integer
    int f = (int)ceil(root);
    int r = power(f, 2) - n;
    bool stop = false;
    // increasing f until r is a perfect square or until f - sqrt(r) < 1
    // worst case: n is prime => n = n * 1 => f - sqrt(r) must be at least 1
    // if f - sqrt(r) < 1, then n must be a prime, then no need to continue
    while (!stop) {
        if (isSquare(r)) {
            stop = true;
        } else {
            double d = f - sqrt(r);
            if (d < 1) {
                stop = true;
            }
            f = f + 1;
            r = power(f, 2) - n;
        }
    }
    // In worst case when n is prime, we can also try f until f = (n+1)/2
    // since n = n * 1 = (f + (f-1)) * (f - (f-1)) => f = (n+1)/2
    //        while (!isSquare(r) && f <= (n+1)/2) {
    //            f = f + 1;
    //            r = power(f, 2) - n;
    //        }
    pair<int, int> factors;
    // if casting to int, result would be rounded (4 down 5 up)
    int rt = (int)sqrt(r);
    int factor1 = f - rt;
    int factor2 = f + rt;
    factors.first = factor1;
    factors.second = factor2;
    return factors;
}

// base b representation
// b = 2
int toBinary(int n) {
    if (n == 0) {
        return 0;
    }
    if (n == 1) {
        return 1;
    }
    int multiplier = n / 2;
    int remainder = n % 2;
    return 10 * toBinary(multiplier) + remainder;
}

// b = 10
int toDecimal(int n) {
    if (n == 0) {
        return 0;
    }
    if (n == 1) {
        return 1;
    }
    int multiplier = n / 10;
    int remainder = n % 10;
    return 2 * toBinary(multiplier) + remainder;
}

bool is_prime(int n) {
    if (n <= 1) {
        return false;
    }
    if (n == 2) {
        return true;
    } else {
        bool check = true;
        for (int i = 2; i < n; i++) {
            if (n % i == 0) {
                return false;
            }
        }
        return check;
    }
}

// compute the max distance between all consecutive two primes in [start, end]
int prime_dist(int start, int end) {
    int first = start;
    while (!is_prime(first) && first <= end) {
        first++;
    }
    if (first > end) {
        cout << "no prime!" << endl;
        return -1;
    }
    int maxDist = 0;

    while (first <= end) {
        int second = first + 1;
        while (second <= end && !is_prime(second)) {
            second++;
        }
        if (second <= end) {
            int currDist = second - first;
            if (currDist > maxDist) {
                maxDist = currDist;
            }
        }
        first = second;
    }
    return maxDist;
}

// row sum
int rowSum(vector<vector<int>> table, int row) {
    int sum = 0;
    for (int i = 0; i < table[row].size(); i++) {
        int curr = table[row][i];
        sum = sum + curr;
    }
    return sum;
}

// col sum
int colSum(vector<vector<int>> table, int col) {
    int sum = 0;
    for (int i = 0; i < table.size(); i++) {
        int curr = table[i][col];
        sum = sum + curr;
    }
    return sum;
}

// smallest digit sensitive prime (last 2 digits)
int digitsensitive() {
    int base = 0;
    int ret = 0;
    while (ret == 0) {
        // create table (10 x 10)
        base = base + 100;
        vector<vector<int>> table;
        for (int i = 0; i <= 9; i++) {
            int row_start = base + 10 * i;
            vector<int> row;
            for (int j = 0; j <= 9; j++) {
                int curr = row_start + j;
                if (is_prime(curr)) {
                    row.push_back(curr);
                } else {
                    row.push_back(0);
                }
            }
            table.push_back(row);
        }

        // search table row & col
        for (int r = 0; r < table.size(); r++) {
            for (int c = 0; c < table[r].size(); c++) {
                int entry = table[r][c];
                if (is_prime(entry) &&
                    rowSum(table, r) == entry &&
                    colSum(table, c) == entry) {
                    ret = entry;
                }
            }
        }
    }

    return ret;
}

//=================================================================
// Coding theory
int hammingDistance(string s1, string s2) {
    int len1 = s1.length();
    int len2 = s2.length();
    if (len1 != len2) {
        cout << "Invalid input" << endl;
        return -1;
    } else {
        int count = 0;
        for (int i = 0; i < len1; i++) {
            if (s1[i] != s2[i]) {
                count = count + 1;
            }
        }
        return count;
    }
}

//=============================================================
// Cryptography

// define alphabet
const vector<string> alphabet = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};

template <class T>
int getIndex(vector<T> v, T e) {
    for (int i = 0; i < v.size(); i++) {
        if (v[i] == e) {
            return i;
        }
    }
    return -1;
}

// shift cipher
// encode with key
string shiftEncode(string p, int key) {
    string c = "";
    for (int i = 0; i < p.length(); i++) {
        string sp = p.substr(i, 1);
        int indexP = getIndex(alphabet, sp);
        int indexC = residue((indexP + key) % 26, 26);
        string sc = alphabet[indexC];
        c = c + sc;
    }
    return c;
}

// decode with key
string shiftDecode(string c, int key) {
    string p = "";
    for (int i = 0; i < c.length(); i++) {
        string sc = c.substr(i, 1);
        int indexC = getIndex(alphabet, sc);
        int indexP = residue((indexC - key) % 26, 26);
        string sp = alphabet[indexP];
        p = p + sp;
    }
    return p;
}

//=============================================================
int main() {
    // clock_t start = clock();
    // auto start = chrono::system_clock::now();

    // wall time
    auto start = chrono::high_resolution_clock::now();

    // cout << mod(24, 2) << endl;
    // cout << mod(100, 13) << endl;
    // cout << gcd(23, 17) << endl;
    // cout << gcdSteps(23, 17) << endl;
    // cout << gcdSteps(15673, 2532) << endl;

    // cout << is_congruent(5, 17, 13) << endl;

    // cout << is_prime(17) << endl;
    // cout << is_prime(2) << endl;
    // cout << is_prime(3) << endl;
    // cout << is_prime(33) << endl;

    // cout << prime_dist(2, 101) << endl;

    // cout << bezoutCoef(173, 167).first << " " << bezoutCoef(173, 167).second << endl;
    // cout << bezoutCoef(603, 626).first << " " << bezoutCoef(603, 626).second << endl;

    // cout << modExp(504, 1109, 2881) << endl;
    // cout << modExp(1874, 1109, 2881) << endl;
    // cout << modExp(347, 1109, 2881) << endl;
    // cout << modExp(515, 1109, 2881) << endl;
    // cout << modExp(2088, 1109, 2881) << endl;
    // cout << modExp(2356, 1109, 2881) << endl;
    // cout << modExp(736, 1109, 2881) << endl;
    // cout << modExp(468, 1109, 2881) << endl;

    // cout << digitsensitive() << endl;

    // clock_t end = clock();
    // auto end = chrono::system_clock::now();
    auto end = chrono::high_resolution_clock::now();

    // double timeTaken = (end - start) / CLOCKS_PER_SEC;
    auto timeTaken = chrono::duration_cast<chrono::milliseconds>(end - start);

    cout << "Time taken: " << timeTaken.count() << " milliseconds" << endl;
    return 0;
}