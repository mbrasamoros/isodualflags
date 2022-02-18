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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
typedef unsigned __int128 bbint;
int gamma;
long long int Sdesc(const bbint G, bbint S, const int c, const int m, const int u, const int v, int r, int gd) {
  int newc, s = 1, sold = 0;
  bbint A;
  long long int ng = 0;
  if (gd) {
    --gd;
    A = G;
    if (S & (bbint)1) {
      if (S & ((bbint)1 << m))
        ng += Sdesc(G | ((bbint)1 << (c - 1)), (S >> 1) | ((bbint)3 << (c - 1)), c + 1, m, u, v, r--, gd);
      else
        ng += Sdesc(G | ((bbint)1 << (c - 1)), (S >> 1) | ((bbint)3 << (c - 1)), c + 1, m, u, v, --r, gd);
    }
    while (r) {
      if (S & (bbint)1 << s) {
        newc = c + s + 1;
        for (int i = sold; i < s; i++) {
          A <<= 1;
          S &= A;
        }
        if (S & ((bbint)1 << (m + s)))
          ng += Sdesc(G | ((bbint)1 << (newc - 2)), S >> (s + 1) | ((bbint)7 << (newc - 3)), newc, m, u, v, r--, gd);
        else
          ng += Sdesc(G | ((bbint)1 << (newc - 2)), S >> (s + 1) | ((bbint)7 << (newc - 3)), newc, m, u, v, --r, gd);
        sold = s;
      }
      s++;
    }
  } else {
    A = S & (S >> m) & (((bbint)1 << u) - 1);
    s = 0;
    while (A) {
      s++;
      A &= A - 1;
    }
    ng += s * (r - 1);
    ng += r * (r - 1) * (r - 2) / 6;
    if (c > m + u) {
      A = (((bbint)1 << v) - 1) & S & (S >> u) & (S >> (m + u));
      while (A) {
        ng++;
        A &= A - 1;
      }
    }
  }
  return ng;
}
long long int lowrank(const int m, const int u, const int gamma) {
  int gd, min, mu, c;
  bbint G, S;
  long long int ng = 0;
  mu = m + u;
  if (u == m) {
    gd = gamma - mu;
    if (m < gd)
      min = m;
    else
      min = gd;
    c = mu;
    for (int v = 1; v < min; v++) {
      gd--;
      G = ((bbint)1 << c) - 1;
      c++;
      G -= ((bbint)1 << (m - 1));
      G -= ((bbint)1 << (mu - 1));
      S = ((bbint)1 << c) - 1;
      S -= ((bbint)1 << (m - v));
      S -= ((bbint)1 << (mu - v));
      ng += Sdesc(G, S - 1, c, m, u, v, m - 2, gd);
    }
    gd--;
    G = ((bbint)1 << c) - 1;
    c++;
    G -= ((bbint)1 << (m - 1));
    G -= ((bbint)1 << (mu - 1));
    S = ((bbint)1 << c) - 1;
    S -= ((bbint)1 << (m - min));
    S -= ((bbint)1 << (mu - min));
    ng += Sdesc(G, S, c, m, u, min, m - 1, gd);
  } else {
    if (m - u < gamma - mu)
      min = m - u;
    else
      min = gamma - mu;
    c = mu;
    gd = gamma - m - u;
    for (int v = 1; v < min; v++) {
      gd--;
      G = ((bbint)1 << c) - 1;
      c++;
      G -= ((bbint)1 << (m - 1));
      G -= ((bbint)1 << (mu - 1));
      S = ((bbint)1 << c) - 1;
      S -= ((bbint)1 << (m - u - v));
      S -= ((bbint)1 << (m - v));
      S -= ((bbint)1 << (mu - v));
      if (u < v)
        ng += Sdesc(G, S - 1, c, m, u, v, m - 4, gd);
      else
        ng += Sdesc(G, S - 1, c, m, u, v, m - 3, gd);
    }
    gd--;
    G = ((bbint)1 << c) - 1;
    c++;
    G -= ((bbint)1 << (m - 1));
    G -= ((bbint)1 << (mu - 1));
    S = ((bbint)1 << c) - 1;
    S -= ((bbint)1 << (m - u - min));
    S -= ((bbint)1 << (m - min));
    S -= ((bbint)1 << (mu - min));
    if (u < min)
      ng += Sdesc(G, S, c, m, u, min, m - 3, gd);
    else
      ng += Sdesc(G, S, c, m, u, min, m - 2, gd);
  }
  return ng;
}
int main(int numvars, char **vars) {
  long long int ng = 0;
  bbint G, S;
  int x, min;
  time_t seconds, secondsafter;
  gamma = (int)atoi(vars[1]);
  seconds = time(NULL);
  ng = (gamma - 4) * (gamma - 5) * (gamma - 6) / 6 + (gamma - 2) * (gamma - 3) + 8; // case m=gamma-2 and the hyperelliptic;
  for (int m = 3; m < gamma - 2; m++) {
    x = gamma - m - 1;
    if (m - 1 < x)
      min = m - 1;
    else
      min = x;
    ng += Sdesc(((bbint)1 << (m - 1)) - 1, ((bbint)1 << m) - 8, m, m, 1, 1, m - 3, x - 1);
    for (int v = 2; v < min; v++) {
      G = ((bbint)1 << (m + v)) - 1;
      G -= ((bbint)3 << (m - 1));
      S = ((bbint)1 << (m + v + 1)) - 1;
      S -= ((bbint)7 << (m - v - 1));
      ng += Sdesc(G, S - 1, m + v + 1, m, 1, v, m - 4, x - v);
    }
    G = ((bbint)1 << (m + min)) - 1;
    G -= ((bbint)3 << (m - 1));
    S = ((bbint)1 << (m + min + 1)) - 1;
    S -= ((bbint)7 << (m - min - 1));
    ng += Sdesc(G, S, m + min + 1, m, 1, min, m - 3, x - min);
    if (x > m) {
      for (int u = 2; u <= m; u++)
        ng += lowrank(m, u, gamma);
    } else {
      for (int u = 2; u < x; u++) {
        ng += lowrank(m, u, gamma);
      }
      ng += (m - 1) * (m - 2) * (m - 3) / 6;
      if (m < 2 * x)
        ng += (x - 1) * (m - 2);
      else
        ng += x * (m - 2);
      if (x == m || 2 * x != m)
        ng++;
      if (x == m - 1)
        ng++;
      if (x < m - 1 && 2 * x != m - 1)
        ng += 2;
    }
  }
  secondsafter = time(NULL);
  printf("n%d=%lld (without parallelization) time taken %d\n", (int)gamma, ng, (int)(secondsafter - seconds));
  return 0;
}
