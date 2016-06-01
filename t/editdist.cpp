#include <iostream>
#include <deque>
#include <vector>
#include <algorithm>

using namespace std;

int main (int argc, char** argv) {
  if (argc != 3 && argc != 4) {
    cout << "Usage: " << argv[0] << " <string1> <string2> [<diag_band_width>]" << endl;
    exit (EXIT_FAILURE);
  }

  const string x (argv[1]);
  const string y (argv[2]);
  
  const int xlen = x.length();
  const int ylen = y.length();

  const int band = argc == 4 ? atoi (argv[3]) : max(xlen*2,ylen*2);

  const int diff = ylen - xlen;
  const int bmin = max (band/2, -diff);
  const int bmax = max (band/2, diff);

  deque<vector<int> > cell;
  int prev_jmin, prev_jmax;
  const int inf = xlen + ylen;  // max possible edit distance
  for (int i = 0; i <= xlen; ++i) {
    const int jmin = max (0, i - bmin);
    const int jmax = min (ylen, i + bmax);
    if (i >= 2)
      cell.pop_front();
    cell.push_back (vector<int> (jmax + 1 - jmin, inf));
    const char xi = i > 0 ? x[i-1] : '?';
    for (int j = jmin; j <= jmax; ++j) {
      const char yj = j > 0 ? y[j-1] : '?';
      int sc = inf;
      if (i == 0 && j == 0)
	sc = 0;
      if (i > 0 && j <= prev_jmax) {
	sc = min (sc, cell.front()[j-prev_jmin] + 1);
	if (j > prev_jmin)
	  sc = min (sc, cell.front()[j-prev_jmin-1] + (xi == yj ? 0 : 1));
      }
      if (j > jmin)
	sc = min (sc, cell.back()[j-jmin-1] + 1);
      cell.back()[j-jmin] = sc;
    }
    prev_jmin = jmin;
    prev_jmax = jmax;
  }

  cout << ((cell.size() > 1 && cell.back().size() > 0) ? cell.back().back() : inf) << endl;

  exit (EXIT_SUCCESS);
}
