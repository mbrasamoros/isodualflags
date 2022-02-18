/*
 *
 * Copyright (c) 2022 
 * Maria Bras-Amoros, Julio Fernandez-Gonzalez, Cesar Marin Rodriguez
 *
 * Distributed under the terms of the GNU General Public License
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * The full text is available at http:
 *
 * Last update: February 16, 2022
 *
 */
#include <atomic>
#include <bitset>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
#define MAXNUMBITS 128
#define MAXGAMMA (MAXNUMBITS >> 1) + 4
typedef unsigned __int128 bbint;
const bbint constant_one = 1;
bbint array_1shiftn[MAXNUMBITS] __attribute__((aligned(16)));
bbint array_3shiftn[MAXNUMBITS] __attribute__((aligned(16)));
bbint array_7shiftn[MAXNUMBITS] __attribute__((aligned(16)));
bbint array_1shiftn_minusone[MAXNUMBITS] __attribute__((aligned(16)));
bbint array_1shiftn_minuseight[MAXNUMBITS] __attribute__((aligned(16)));
typedef struct
{
  bbint G;
  bbint S;
  uint_fast64_t c;
  uint_fast64_t m;
  uint_fast64_t u;
  uint_fast64_t v;
  uint_fast64_t r;
  uint_fast64_t gd;
  uint64_t ng;
} Tentry;
Tentry t[1024 * 100];
uint_fast64_t indexi = 0;
uint64_t millis() {
  struct timespec spec;
  clock_gettime(CLOCK_MONOTONIC, &spec);
  return spec.tv_sec * 1000 + spec.tv_nsec / 1e6;
}
int_fast64_t Sdesc(const bbint G, bbint S, const uint_fast64_t c, const uint_fast64_t m, const uint_fast64_t u, const uint_fast64_t v, uint_fast64_t r, uint_fast64_t gd) {
  uint_fast64_t newc, ss = 1, sold = 0;
  int_fast64_t ng = 0;
  bbint A, AA;
  if (gd) {
    --gd;
    A = G;
    if (S & constant_one) {
      if (S & array_1shiftn[m])
        ng = Sdesc(G | array_1shiftn[c - 1], (S >> 1) | array_3shiftn[c - 1], c + 1, m, u, v, r--, gd);
      else
        ng = Sdesc(G | array_1shiftn[c - 1], (S >> 1) | array_3shiftn[c - 1], c + 1, m, u, v, --r, gd);
    }
    while (r) {
      if (S & array_1shiftn[ss]) {
        newc = c + ss + 1;
        for (uint_fast64_t i = sold; i < ss; i++) {
          A <<= 1;
          S &= A;
        }
        if (S & array_1shiftn[m + ss])
          ng += Sdesc(G | array_1shiftn[newc - 2], S >> (ss + 1) | array_7shiftn[newc - 3], newc, m, u, v, r--, gd);
        else
          ng += Sdesc(G | array_1shiftn[newc - 2], S >> (ss + 1) | array_7shiftn[newc - 3], newc, m, u, v, --r, gd);
        sold = ss;
      }
      ss++;
    }
  } else {
    A = S & (S >> m) & array_1shiftn_minusone[u];
    ss = 0;
    while (A) {
      ss++;
      A &= A - 1;
    }
    ng += ss * (r - 1);
    ng += r * (r - 1) * (r - 2) / 6;
    if (c > m + u) {
      AA = array_1shiftn_minusone[v] & S & (S >> u) & (S >> (m + u));
      while (AA) {
        ng++;
        AA &= AA - 1;
      }
    }
  }
  return ng;
}
void fillentry(uint_fast64_t i, const bbint G, bbint S, const uint_fast64_t c, const uint_fast64_t m, const uint_fast64_t u, const uint_fast64_t v, uint_fast64_t r, uint_fast64_t gd) {
  t[i].G = G;
  t[i].S = S;
  t[i].c = c;
  t[i].m = m;
  t[i].u = u;
  t[i].v = v;
  t[i].r = r;
  t[i].gd = gd;
  t[i].ng = 0;
}
void lowrank(const uint_fast64_t m, const uint_fast64_t u, const uint_fast64_t gamma) {
  uint_fast64_t gd, min, mu, c;
  bbint G, S;
  mu = m + u;
  if (u == m) {
    gd = gamma - mu;
    if (m < gd)
      min = m;
    else
      min = gd;
    c = mu;
    for (uint_fast64_t v = 1; v < min; v++) {
      gd--;
      G = array_1shiftn_minusone[c];
      c++;
      G -= array_1shiftn[m - 1];
      G -= array_1shiftn[mu - 1];
      S = array_1shiftn_minusone[c];
      S -= array_1shiftn[m - v];
      S -= array_1shiftn[mu - v];
      fillentry(indexi++, G, S - 1, c, m, u, v, m - 2, gd);
    }
    gd--;
    G = array_1shiftn_minusone[c];
    c++;
    G -= array_1shiftn[m - 1];
    G -= array_1shiftn[mu - 1];
    S = array_1shiftn_minusone[c];
    S -= array_1shiftn[m - min];
    S -= array_1shiftn[mu - min];
    fillentry(indexi++, G, S, c, m, u, min, m - 1, gd);
  } else {
    if (m - u < gamma - mu)
      min = m - u;
    else
      min = gamma - mu;
    c = mu;
    gd = gamma - m - u;
    for (uint_fast64_t v = 1; v < min; v++) {
      gd--;
      G = array_1shiftn_minusone[c];
      c++;
      G -= array_1shiftn[m - 1];
      G -= array_1shiftn[mu - 1];
      S = array_1shiftn_minusone[c];
      S -= array_1shiftn[m - u - v];
      S -= array_1shiftn[m - v];
      S -= array_1shiftn[mu - v];
      if (u < v)
        fillentry(indexi++, G, S - 1, c, m, u, v, m - 4, gd);
      else
        fillentry(indexi++, G, S - 1, c, m, u, v, m - 3, gd);
    }
    gd--;
    G = array_1shiftn_minusone[c];
    c++;
    G -= array_1shiftn[m - 1];
    G -= array_1shiftn[mu - 1];
    S = array_1shiftn_minusone[c];
    S -= array_1shiftn[m - u - min];
    S -= array_1shiftn[m - min];
    S -= array_1shiftn[mu - min];
    if (u < min)
      fillentry(indexi++, G, S, c, m, u, min, m - 3, gd);
    else
      fillentry(indexi++, G, S, c, m, u, min, m - 2, gd);
  }
}
int main(int numvars, char **vars) {
  int_fast64_t ng;
  uint_fast64_t gamma;
  long seconds, secondsafter;
  __cilkrts_set_param("nworkers", "24");
  gamma = (uint_fast64_t)atoi(vars[1]);
  if (gamma > MAXGAMMA || gamma < 7) {
    cout << "target genus must be >= 7 and <= " << MAXGAMMA << endl;
    exit(1);
  }
  seconds = millis();
  cilk_for(uint_fast64_t u = 0; u < MAXNUMBITS; u++) {
    array_1shiftn[u] = (bbint)1 << u;
    array_3shiftn[u] = (bbint)3 << u;
    array_7shiftn[u] = (bbint)7 << u;
    array_1shiftn_minusone[u] = ((bbint)1 << u) - 1;
    array_1shiftn_minuseight[u] = ((bbint)1 << u) - 8;
  }
  ng = (gamma - 4) * (gamma - 5) * (gamma - 6) / 6 + (gamma - 2) * (gamma - 3) + 8;
  for (uint_fast64_t m = 3; m < gamma - 2; m++) {
    uint_fast64_t x, min;
    bbint G, S;
    x = gamma - m - 1;
    if (m - 1 < x)
      min = m - 1;
    else
      min = x;
    fillentry(indexi++, array_1shiftn_minusone[m - 1], array_1shiftn_minuseight[m], m, m, 1, 1, m - 3, x - 1);
    for (uint_fast64_t v = 2; v < min; v++) {
      bbint GG, SS;
      GG = array_1shiftn_minusone[m + v];
      GG -= array_3shiftn[m - 1];
      SS = array_1shiftn_minusone[m + v + 1];
      SS -= array_7shiftn[m - v - 1];
      fillentry(indexi++, GG, SS - 1, m + v + 1, m, 1, v, m - 4, x - v);
    }
    G = array_1shiftn_minusone[m + min];
    G -= array_3shiftn[m - 1];
    S = array_1shiftn_minusone[m + min + 1];
    S -= array_7shiftn[m - min - 1];
    fillentry(indexi++, G, S, m + min + 1, m, 1, min, m - 3, x - min);
    if (x > m) {
      for (uint_fast64_t u = 2; u <= m; u++)
        lowrank(m, u, gamma);
    } else {
      for (uint_fast64_t u = 2; u < x; u++)
        lowrank(m, u, gamma);
      int_fast64_t tempi;
      tempi = (m - 1) * (m - 2) * (m - 3) / 6;
      if (m < 2 * x)
        tempi += (x - 1) * (m - 2);
      else
        tempi += x * (m - 2);
      if (x == m || 2 * x != m)
        tempi++;
      if (x == m - 1)
        tempi++;
      if (x < m - 1 && 2 * x != m - 1)
        tempi += 2;
      ng += tempi;
    }
  }
#pragma cilk grainsize = 1;
  cilk_for(uint_fast64_t a = 0; a < indexi; a++) {
    t[a].ng = Sdesc(t[a].G, t[a].S, t[a].c, t[a].m, t[a].u, t[a].v, t[a].r, t[a].gd);
  }
  for (uint_fast64_t a = 0; a < indexi; a++) {
    ng += t[a].ng;
  }
  secondsafter = millis();
  cout << "n" << (int)gamma << "=" << ng << " (" << (int)(__cilkrts_get_nworkers()) << " workers) time taken " << (int)(secondsafter - seconds) << " miliseconds" << endl;
  return 0;
}
