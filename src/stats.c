#include "lise.h"

void initstats(STATS *restrict st) {
    st->max = MAX;
    st->min = MIN;
    st->range = 0;
}
void setstats(int n, STATS *restrict st) {
    if (n > st->max) {
        st->max = n;
        return;
    } else if ((n < st->min)) {
        st->min = n;
        return;
    }
}
void finstats(STATS *restrict st) {
    st->max = floor(st->max + ast) + 1;
	st->min = floor(st->min - ast);
	st->range = st->max - st->min+1;
}