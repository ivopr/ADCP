/* =============================================================================
 * cyclic_energy() must restrain the macrocycle to the peptide-bond geometry.
 *
 * This is the check behind the decision recorded in docs/cyclic-port.md: the
 * closure model was taken from upstream's `cyclic` branch, but with two
 * arithmetic defects repaired. Ported verbatim, upstream's version reduces to
 *
 *     ans += CaDistance;        // CaDistance is d^2
 *
 * for every geometry a peptide actually visits -- a harmonic centred on ZERO,
 * which drags the two termini together until van der Waals repulsion stops
 * them. Measured over the 18 backbone-cyclic targets, the shipped ADFRsuite
 * binary closes the ring at a median CA1-CAn of 2.52 A against 3.83 A for the
 * crystallographic ligands.
 *
 * The assertions below fail on that behaviour and pass on the repaired one.
 * =============================================================================
 */
#include <cmath>
#include <cstdio>

#include <vector>

#include "error.h"
#include "params.h"
#include "vector.h"
#include "rotation.h"
#include "aadict.h"
#include "peptide.h"
#include "energy.h"

/* Release builds define NDEBUG, which compiles assert() away entirely -- so
 * this check does its own reporting rather than relying on it. */
static int failures = 0;
#define CHECK(cond, msg)                                                      \
	do {                                                                  \
		if (!(cond)) {                                                \
			fprintf(stderr, "FAIL: %s\n      (%s)\n", msg, #cond); \
			failures++;                                           \
		}                                                             \
	} while (0)

/* CA1-CAn distance the model restrains to, and the N-C peptide bond length. */
static const double CA_EQUILIBRIUM = 3.819;
static const double NC_EQUILIBRIUM = 1.345;

/* Two residues whose CA-CA separation is `d`, with the N-C pair parked at its
 * own ideal length far away, so only the CA term varies across the scan. */
static double closure_energy(double d)
{
	AA a{}, b{};

	b.ca[0] = 0.0;      b.ca[1] = 0.0; b.ca[2] = 0.0;
	a.ca[0] = d;        a.ca[1] = 0.0; a.ca[2] = 0.0;

	b.c[0]  = 0.0;            b.c[1] = 0.0; b.c[2] = 100.0;
	a.n[0]  = NC_EQUILIBRIUM; a.n[1] = 0.0; a.n[2] = 100.0;

	return cyclic_energy(&a, &b, 0);
}

int main(void)
{
	/* 1. The minimum sits on the peptide-bond CA-CA distance, not on zero.
	 *    Upstream's verbatim version is monotone increasing from d = 0, so its
	 *    minimum over this scan lands at the first sample. */
	double best_d = 0.0, best_E = 1e300;
	for (int i = 0; i <= 11000; i++) {
		double d = 1.0 + i * 0.001;          /* 1.0 .. 12.0 A */
		double E = closure_energy(d);
		if (E < best_E) { best_E = E; best_d = d; }
	}
	printf("minimum at CA1-CAn = %.3f A, E = %.6g\n", best_d, best_E);
	CHECK(fabs(best_d - CA_EQUILIBRIUM) < 0.005,
	      "cyclic_energy() is not minimised at the peptide-bond CA-CA distance");
	CHECK(fabs(best_E) < 1e-9,
	      "the restraint does not vanish at its own equilibrium");

	/* 2. The collapsed ring the shipped binary produces must cost more than
	 *    the crystallographic one. */
	double E_collapsed = closure_energy(2.52);   /* ADFRsuite 1.0 median */
	double E_crystal   = closure_energy(3.83);   /* crystallographic median */
	printf("E(2.52 A) = %.6g   E(3.83 A) = %.6g\n", E_collapsed, E_crystal);
	CHECK(E_collapsed > E_crystal,
	      "a 2.52 A ring is not penalised against a 3.83 A one");

	/* 3. Continuous where the far-field d^2 pull switches on. Upstream's own
	 *    cutoff steps by ~23 there; ours subtracts the boundary value. */
	double below = closure_energy(5.0 - 1e-6);
	double above = closure_energy(5.0 + 1e-6);
	printf("E(5-) = %.6g   E(5+) = %.6g   step = %.3g\n",
	       below, above, fabs(above - below));
	CHECK(fabs(above - below) < 1e-3,
	      "the closure restraint steps at the 5 A boundary");

	/* 4. Strictly increasing on both sides of the equilibrium, so the search
	 *    always has a gradient pointing back to it. */
	int up_ok = 1, down_ok = 1;
	for (double d = CA_EQUILIBRIUM; d < 11.99; d += 0.01)
		if (closure_energy(d + 0.01) <= closure_energy(d)) up_ok = 0;
	for (double d = CA_EQUILIBRIUM; d > 1.01; d -= 0.01)
		if (closure_energy(d - 0.01) <= closure_energy(d)) down_ok = 0;
	CHECK(up_ok, "restraint is not monotone above the equilibrium");
	CHECK(down_ok, "restraint is not monotone below the equilibrium");

	if (failures) {
		fprintf(stderr, "%d check(s) failed\n", failures);
		return 1;
	}
	printf("PASS: ring closure restrained to %.3f A, continuous, monotone\n",
	       CA_EQUILIBRIUM);
	return 0;
}
