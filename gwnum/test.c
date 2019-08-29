#include "gwnum.h"

int main () {
	gwhandle gwdata;
	gwnum	x;

	gwinit (&gwdata);
	gwsetup (&gwdata, 1.0, 2, 640, -1);

	x = gwalloc (&gwdata);
	dbltogw (&gwdata, 2.0, x);
	gwsetnormroutine (&gwdata, 0, 1, 0); /* Enable error checking */
	gwsquare (&gwdata, x);
}
