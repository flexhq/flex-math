// Include flex header 
#include "flex-core/flex-core.h"

#define FLEX_TRUE 1u
#define FLEX_FALSE 0u
#define FLEX_BOOL uint8_t

// const double M_PI = 3.141592653589793
// const double M_E = 2.718281828459045

fObject* fMath_acos(fVM* vm, double x) {
	return fVM_PushDouble(vm, acos(x));
}

fObject* fMath_acosh(fVM* vm, double x) {
	return fVM_PushDouble(vm, acosh(x));
}

fObject* fMath_asin(fVM* vm, double x) {
	return fVM_PushDouble(vm, asin(x));
}

fObject* fMath_asinh(fVM* vm, double x) {
	return fVM_PushDouble(vm, asinh(x));
}

fObject* fMath_atan(fVM* vm, double x) {
	return fVM_PushDouble(vm, atan(x));
}

fObject* fMath_atan2(fVM* vm, double dx, double dy) {
	// return fVM_PushDouble(vm, atan2(dx, dy));

	// printf("dx: %f dy: %f\n", dx, dy);
	if (dx == 0.0 && dy == 0.0) {
		if (!signbit(dx))
			return fVM_PushDouble(vm, dy);
		if (!signbit(dy))
			return fVM_PushDouble(vm, M_PI);
		return fVM_PushDouble(vm, -M_PI);
	}
	return fVM_PushDouble(vm, atan2(dx, dy));
}

fObject* fMath_atanh(fVM* vm, double x) {
	/* check for domain error */
	if (x <  -1.0 || +1.0 < x) throw(RuntimeException, "Domain error: atanh");
	/* check for pole error */
	if (x == -1.0) return fVM_PushDouble(vm, -INFINITY);
	if (x == +1.0) return fVM_PushDouble(vm, +INFINITY);
	return fVM_PushDouble(vm, atanh(x));
}

fObject* fMath_cbrt(fVM* vm, double f) {
	double r = cbrt(f);
#if defined __GLIBC__
	if (isfinite(r)) {
		r = (2.0 * r + (f / r / r)) / 3.0;
	}
#endif
	return fVM_PushDouble(vm, r);
}

fObject* fMath_cos(fVM* vm, double x) {
	return fVM_PushDouble(vm, cos(x));
}

fObject* fMath_cosh(fVM* vm, double x) {
	return fVM_PushDouble(vm, cosh(x));
}

fObject* fMath_erf(fVM* vm, double x) {
	return fVM_PushDouble(vm, erf(x));
}

fObject* fMath_erfc(fVM* vm, double x) {
	return fVM_PushDouble(vm, erfc(x));
}

fObject* fMath_exp(fVM* vm, double x) {
	return fVM_PushDouble(vm, exp(x));
}

fObject* fMath_frexp(fVM* vm, double x) {
	int exp;

	double d = frexp(x, &exp);

	return fVM_PushVector(vm, 2, fVM_PushDouble(vm, d), fVM_PushInt32(vm, exp));

	// return fVM_PushDouble(vm, exp(x));
}

#define numberof(ary) (int)(sizeof(ary)/sizeof((ary)[0]))
static const double fact_table[] = {
	/* fact(0) */ 1.0,
	/* fact(1) */ 1.0,
	/* fact(2) */ 2.0,
	/* fact(3) */ 6.0,
	/* fact(4) */ 24.0,
	/* fact(5) */ 120.0,
	/* fact(6) */ 720.0,
	/* fact(7) */ 5040.0,
	/* fact(8) */ 40320.0,
	/* fact(9) */ 362880.0,
	/* fact(10) */ 3628800.0,
	/* fact(11) */ 39916800.0,
	/* fact(12) */ 479001600.0,
	/* fact(13) */ 6227020800.0,
	/* fact(14) */ 87178291200.0,
	/* fact(15) */ 1307674368000.0,
	/* fact(16) */ 20922789888000.0,
	/* fact(17) */ 355687428096000.0,
	/* fact(18) */ 6402373705728000.0,
	/* fact(19) */ 121645100408832000.0,
	/* fact(20) */ 2432902008176640000.0,
	/* fact(21) */ 51090942171709440000.0,
	/* fact(22) */ 1124000727777607680000.0,
	/* fact(23)=25852016738884976640000 needs 56bit mantissa which is
	 * impossible to represent exactly in IEEE 754 double which have
	 * 53bit mantissa. */
};
enum {NFACT_TABLE = numberof(fact_table)};

fObject* fMath_gamma(fVM* vm, double d) {
	/* check for domain error */
	if (isinf(d)) {
		if (signbit(d)) throw(RuntimeException, "Domain error: gamma");
		return fVM_PushDouble(vm, HUGE_VAL);
	}
	if (d == 0.0) {
		return signbit(d) ? fVM_PushDouble(vm, -HUGE_VAL) : fVM_PushDouble(vm, HUGE_VAL);
	}
	if (d == floor(d)) {
		if (d < 0.0) throw(RuntimeException, "Domain error: gamma");
		if (1.0 <= d && d <= (double)NFACT_TABLE) {
			return fVM_PushDouble(vm, fact_table[(int)d - 1]);
		}
	}
	return fVM_PushDouble(vm, tgamma(d));
}

fObject* fMath_hypot(fVM* vm, double dx, double dy) {
	return fVM_PushDouble(vm, hypot(dx, dy));
}

fObject* fMath_ldexp(fVM* vm, double x, int32_t n) {
	return fVM_PushDouble(vm, ldexp(x, n));
}

fObject* fMath_lgamma(fVM* vm, double d) {
	int sign=1;
	fObject* v;

	/* check for domain error */
	if (isinf(d)) {
		if (signbit(d)) throw(RuntimeException, "Domain error: lgamma");
		return fVM_PushVector(vm, 2, fVM_PushDouble(vm, HUGE_VAL), fVM_PushInt32(vm, 1));
	}
	if (d == 0.0) {
		fObject* vsign = signbit(d) ? fVM_PushInt32(vm, -1) : fVM_PushInt32(vm, +1);
		return fVM_PushVector(vm, 2, fVM_PushDouble(vm, HUGE_VAL), vsign);
	}
	v = fVM_PushDouble(vm, lgamma_r(d, &sign));
	return fVM_PushVector(vm, 2, v, fVM_PushInt32(vm, sign));
}

fObject* fMath_log10(fVM* vm, double d) {
	size_t numbits = 0;

	/* check for domain error */
	if (d < 0.0) throw(RuntimeException, "Domain error: log10");
	/* check for pole error */
	if (d == 0.0) return fVM_PushDouble(vm, -HUGE_VAL);

	return fVM_PushDouble(vm, log10(d) + numbits * log10(2)); /* log10(d * 2 ** numbits) */
}

fObject* fMath_log2(fVM* vm, double d) {
	size_t numbits = 0;

	/* check for domain error */
    if (d < 0.0) throw(RuntimeException, "Domain error: log2");
    /* check for pole error */
    if (d == 0.0) return fVM_PushDouble(vm, -HUGE_VAL);

    return fVM_PushDouble(vm, log2(d) + numbits); /* log2(d * 2 ** numbits) */
}

fObject* fMath_sin(fVM* vm, double x) {
	return fVM_PushDouble(vm, sin(x));
}

fObject* fMath_sinh(fVM* vm, double x) {
	return fVM_PushDouble(vm, sinh(x));
}

fObject* fMath_sqrt(fVM* vm, double x) {
	return fVM_PushDouble(vm, sqrt(x));
}

fObject* fMath_tan(fVM* vm, double x) {
	return fVM_PushDouble(vm, tan(x));
}

fObject* fMath_tanh(fVM* vm, double x) {
	return fVM_PushDouble(vm, tanh(x));
}