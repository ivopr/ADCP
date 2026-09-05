/* =============================================================================
 * cyclic_energy() must reproduce upstream `cyclic`'s closure model EXACTLY --
 * defects included.
 *
 * This is not the physically correct model, and the assertions below are not
 * describing a peptide. This repo tracks the behaviour of the reference
 * implementation rather than a corrected version of it, and upstream's
 * cyclic_energy() reads:
 *
 *     //NCDistance = distance(a->n, b->c);   // assignment commented out
 *     if (CaDistance < 5) {
 *             ans += (sqrt(CaDistance) - 3.819)^2;
 *             ans += (sqrt(NCDistance) - 1.345)^2;
 *     } else ans += CaDistance;
 *
 * NCDistance is never assigned, so the second term is the constant 1.809
 * whatever the geometry; and distance() returns the SQUARE of the distance, so
 * `CaDistance < 5` means d < 2.24 A and the harmonic branch essentially never
 * runs. What executes is `ans += CaDistance`, i.e. d^2 -- a harmonic centred on
 * ZERO. Measured over the 18 backbone-cyclic targets it closes the ring at a
 * median CA1-CAn of ~2.5 A against 3.83 A for the crystallographic ligands.
 *
 * The point of this test is that the behaviour cannot drift silently: it fails
 * if someone repairs the arithmetic, and it fails if someone changes the
 * weights. See docs/cyclic-port.md for the decision and docs/compares/7.md for
 * what it costs.
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

/* Where upstream's `CaDistance < 5` test actually switches, in real angstroms:
 * distance() returns d^2, so the branch flips at d = sqrt(5). */
static const double UPSTREAM_SWITCH = 2.2360679;

/* Two residues whose CA-CA separation is `d`. `nc` is the N-to-C distance,
 * which upstream never reads -- the third argument exists so a test can prove
 * that. */
static double closure_energy(double d, double nc)
{
	AA a{}, b{};

	b.ca[0] = 0.0; b.ca[1] = 0.0; b.ca[2] = 0.0;
	a.ca[0] = d;   a.ca[1] = 0.0; a.ca[2] = 0.0;

	b.c[0]  = 0.0; b.c[1]  = 0.0; b.c[2]  = 100.0;
	a.n[0]  = nc;  a.n[1]  = 0.0; a.n[2]  = 100.0;

	return cyclic_energy(&a, &b, 0);
}

int main(void)
{
	/* 1. The N-C term is geometry-independent -- upstream's commented-out
	 *    assignment. Moving the nitrogen must change nothing at all. */
	double nc_near = closure_energy(3.819, 1.345);
	double nc_far  = closure_energy(3.819, 40.0);
	printf("E(d=3.819, nc=1.345) = %.6g   E(d=3.819, nc=40) = %.6g\n",
	       nc_near, nc_far);
	CHECK(nc_near == nc_far,
	      "the N-C term responds to geometry; upstream's is a constant");

	/* 2. The minimum sits at upstream's branch switch, not at the
	 *    peptide-bond CA-CA distance. A repaired model puts it at 3.819. */
	double best_d = 0.0, best_E = 1e300;
	for (int i = 0; i <= 11000; i++) {
		double d = 1.0 + i * 0.001;          /* 1.0 .. 12.0 A */
		double E = closure_energy(d, 1.345);
		if (E < best_E) { best_E = E; best_d = d; }
	}
	printf("minimum at CA1-CAn = %.3f A, E = %.6g\n", best_d, best_E);
	CHECK(fabs(best_d - UPSTREAM_SWITCH) < 0.005,
	      "cyclic_energy() is not minimised at upstream's sqrt(5) switch");

	/* 3. The collapse itself: the ring the shipped binary produces must cost
	 *    LESS than the crystallographic one. This is the assertion that
	 *    inverts if the model is repaired. */
	double E_collapsed = closure_energy(2.52, 1.345);   /* ADFRsuite median */
	double E_crystal   = closure_energy(3.83, 1.345);   /* crystal median   */
	printf("E(2.52 A) = %.6g   E(3.83 A) = %.6g\n", E_collapsed, E_crystal);
	CHECK(E_collapsed < E_crystal,
	      "a 2.52 A ring is not preferred over a 3.83 A one");

	/* 4. Above the switch the restraint is d^2 about zero, so it rises
	 *    monotonically all the way out with no equilibrium anywhere. */
	int mono = 1;
	for (double d = UPSTREAM_SWITCH + 0.01; d < 11.99; d += 0.01)
		if (closure_energy(d + 0.01, 1.345) <= closure_energy(d, 1.345))
			mono = 0;
	CHECK(mono, "restraint is not monotone above the switch");

	/* 5. And it is exactly d^2 there -- no equilibrium offset survived. */
	double e6 = closure_energy(6.0, 1.345);
	printf("E(6.0 A) = %.6g   (d^2 = 36)\n", e6);
	CHECK(fabs(e6 - 36.0) < 1e-9,
	      "the far branch is not the bare d^2 upstream computes");

	if (failures) {
		fprintf(stderr, "%d check(s) failed\n", failures);
		return 1;
	}
	printf("PASS: upstream `cyclic` closure model reproduced, collapse included\n");
	return 0;
}
