#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define boolean int
#define true 1
#define false 0

// enum "type" {...} "variable name"
enum BOOLEAN { TRUE = 1,
               FALSE = 0 } bool;

// enum BOOLEAN bool;

// MATH 312
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

int congruent(int a, int b, int n) {
    return (((a - b) % n) == 0);
}

// coding theory
int hammingDistance(char s1[], char s2[]) {
    int len1 = strlen(s1);
    int len2 = strlen(s2);
    if (len1 != len2) {
        printf("Invalid input");
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

// boolean millerTest(int n, int b) {
//     boolean composite = false;
//     int p = n - 1;
//     int r = (b^p) % n;
//     if (!congruent(r, 1, n)) {
//         composite = true;
//     } else {
//         r = b^(p/2) % n;
//         while (congruent(r, 1, n)) {
//             if ((p/2) % 2 == 0) {
//                 r = b^(p/2) % n;
//             } else {
//                 return millerTest(n, b+1);
//             }

//         }
//         if (congruent(r, -1, n)) {
//              return millerTest(n, b+1);
//         }
//         if (!congruent(r, 1, n)) {
//             composite = true;
//         }
//     }
// }

void main() {
    clock_t start = clock();

    printf("%d \n", mod(24, 2));
    printf("%d \n", mod(100, 13));
    printf("%d \n", gcd(23, 17));
    printf("%d \n", gcdSteps(23, 17));
    printf("%d \n", gcdSteps(15673, 2532));

    printf("%d \n", congruent(5, 17, 13));

    clock_t end = clock();
    double timeTaken = (end - start) / CLOCKS_PER_SEC;
    printf("Time taken: %.2fs\n", timeTaken);
}