/* Exhaustive search over all subspaces K of codimension c in row(R) = F_2^11 (order-24 Paley matrix).
   K = ker M for a c x 11 matrix M of rank c in reduced row echelon form; the slices of K are the fibres of
   r -> M r on the 24 rows. Counts the K whose slice weights are all >= 3, and those whose weights are all
   divisible by 3 (Proposition 9.1), and prints the weight patterns found.
   Usage: h24_min c < rows.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define D 11
static unsigned X[24];
int main(int argc, char **argv) {
    int c = atoi(argv[1]);
    for (int i = 0; i < 24; i++) if (scanf("%u", &X[i]) != 1) return 1;
    long long total = 0, ge3 = 0, div3 = 0;
    int npiv = 0, pivs[2048][4];
    for (int a = 0; a < (1 << D); a++) if (__builtin_popcount(a) == c) {
        int k = 0; for (int j = 0; j < D; j++) if ((a >> j) & 1) pivs[npiv][k++] = j;
        npiv++;
    }
    #pragma omp parallel for schedule(dynamic,1) reduction(+:total,ge3,div3)
    for (int t = 0; t < npiv; t++) {
        int *P = pivs[t];
        int isp[D] = {0}; for (int k = 0; k < c; k++) isp[P[k]] = 1;
        int fr[64], fc[64], nf = 0;      /* free positions (row k, column j > P[k], j not a pivot) */
        for (int k = 0; k < c; k++) for (int j = P[k] + 1; j < D; j++) if (!isp[j]) { fr[nf] = k; fc[nf] = j; nf++; }
        for (long long mask = 0; mask < (1LL << nf); mask++) {
            unsigned m[4] = {0, 0, 0, 0};
            for (int k = 0; k < c; k++) m[k] = 1u << P[k];
            for (int f = 0; f < nf; f++) if ((mask >> f) & 1) m[fr[f]] |= 1u << fc[f];
            int cnt[16] = {0};
            for (int i = 0; i < 24; i++) {
                int s = 0;
                for (int k = 0; k < c; k++) s |= (__builtin_parity(m[k] & X[i])) << k;
                cnt[s]++;
            }
            total++;
            int ok = 1, ok3 = 1;
            for (int s = 0; s < (1 << c); s++) {
                if (cnt[s] == 1 || cnt[s] == 2) ok = 0;
                if (cnt[s] % 3) ok3 = 0;
            }
            if (ok) {
                ge3++;
                if (ok3) div3++;
                #pragma omp critical
                {
                    if (ge3 <= 5 || ok3) {
                        printf("all weights >= 3:");
                        for (int s = 0; s < (1 << c); s++) if (cnt[s]) printf(" %d", cnt[s]);
                        printf("%s\n", ok3 ? "  [all divisible by 3]" : "");
                    }
                }
            }
        }
    }
    printf("codim %d: %lld subspaces, %lld with all slice weights >= 3, %lld with all weights divisible by 3\n",
           c, total, ge3, div3);
    return 0;
}
