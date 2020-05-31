

#ifndef mylbfgsb_h
#define mylbfgsb_h

void optimLBFGSB(int p,
                 void (*)(int, double*, double*),
                 void (*)(int, double*, double*),
                 double*,
                 double*,
                 double*,
                 int,
                 double,
                 double*,
                 double*);

#endif /* mylbfgsb_h*/
