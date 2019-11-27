
long int ipow(long int base, unsigned int exp);
double complex dot_product(double complex *v, double complex *w, long int N, unsigned int D);
void scalar_vec(double complex *out, double complex *vec, double complex *s, long int N, unsigned int D);
void add_vec(double complex *out, double complex *w double complex *v, long int N, unsigned int D);
void sub_vec(double complex *out, double complex *w double complex *v, long int N, unsigned int D);
void assign_vec(double complex *out, double complex *in, long int N, unsigned int D);
void cg(double complex *out, void (*f)(double complex *, double complex *), double complex *in, int max_iter, double tol, long int N, unsigned int D);
