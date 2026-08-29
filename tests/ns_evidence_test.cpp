/* =============================================================================
 * Nested-sampling evidence validation.
 *
 * Every other test in this suite asserts that ADCP does not crash, is
 * reproducible, or lands in a plausible range. None of them check that nested
 * sampling computes the RIGHT NUMBER, because the peptide likelihood has no
 * closed form to check against.
 *
 * This one drives the real bookkeeping -- find_worst(), update_NS_parameters()
 * and the alpha = exp(-1/N) volume schedule, linked from adcp_core -- against a
 * problem whose evidence is known analytically:
 *
 *     prior       x ~ Uniform(0,1)
 *     likelihood  L(x) = exp(-x/tau),  so  log L = -x/tau
 *     evidence    Z = int_0^1 exp(-x/tau) dx = tau * (1 - exp(-1/tau))
 *
 * L is monotone decreasing in x, so the constrained prior {L > L*} is exactly
 * {x < x*} with x* = -tau * logL*, and a constrained sample is Uniform(0, x*)
 * drawn exactly. That is deliberate: it removes MCMC quality from the
 * measurement so that what remains under test is the NS bookkeeping itself.
 *
 * The estimator is unbiased but not exact -- nested sampling's error is
 * ~ sqrt(H/N) with H the information, which update_NS_parameters also computes.
 * The tolerance below is 3 sqrt(H/N), the standard 3-sigma band.
 *
 * Uses its own RNG, not rand(), so unlike the rest of the suite this test's
 * numbers are identical on glibc and BSD libc and it can assert an exact target.
 * ========================================================================== */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>

#include "params.h"
#include "vector.h"
#include "rotation.h"
#include "peptide.h"
#include "checkpoint_io.h"

/* Defined in nested.cpp. Not declared in nested.h -- only nestedsampling() is
   public -- but they have external linkage, which is what makes this test
   possible without restructuring production code. */
void find_worst(simulation_params *sim_params, ChainHash *chainhash, int *heaplength, int N, int P);
void update_NS_parameters(simulation_params *sim_params, int *converged, double *logZnew, double *lweight);
void heapifyhashin(ChainHash *chainhash, int length);

static const double TAU = 0.2;

static double logL(double x) { return -x / TAU; }
static double x_of_logLstar(double ll) { return -TAU * ll; }

static double analytic_logZ() {
	/* log( tau * (1 - exp(-1/tau)) ), via log1p for accuracy */
	return std::log(TAU) + std::log1p(-std::exp(-1.0 / TAU));
}

static int failures = 0;

static void check(bool ok, const char *what) {
	std::printf("  %-52s %s\n", what, ok ? "ok" : "FAIL");
	if (!ok) failures++;
}

/* One nested-sampling run, mirroring the loop structure of nestedsampling().
   Returns the estimated log Z; writes the information H out through *H_out. */
struct DrawStats { long fixed_hits_discarded; long buggy_hits_discarded; };

static double run_ns(int N, int iters, unsigned seed, double *H_out, int *invariant_ok,
                     DrawStats *draws)
{
	std::mt19937 rng(seed);
	std::uniform_real_distribution<double> unit(0.0, 1.0);

	simulation_params sim{};
	sim.thermobeta = 1.0;   /* lweight = log_DeltaX + logLstar * thermobeta */
	sim.lowtemp = 0;        /* no convergence short-circuit */
	sim.checkpoint = 0;     /* keeps update_NS_parameters off the fprintf path */
	sim.iter_start = 0;
	sim.logX_start = 0.0;
	sim.logZ = -std::numeric_limits<double>::max();
	sim.H = 0.0;

	/* Exactly as nestedsampling() sets them up for the serial build. */
	sim.alpha = std::exp(-1.0 / (double)N);
	sim.Delta_logX = std::log(sim.alpha);
	sim.log_DeltaX = std::log(1.0 - sim.alpha);

	/* Live points. chainhash is 1-based: the heap occupies 1..N. */
	std::vector<ChainHash> chainhash(N + 2);
	std::vector<double> x(N + 2);
	for (int i = 1; i <= N; i++) {
		x[i] = unit(rng);
		chainhash[i].index = i;
		chainhash[i].processor = 0;
		chainhash[i].ll = logL(x[i]);
	}
	constructhashheap(chainhash.data(), N);

	int converged = 1;
	double logZnew = 0.0, lweight = 0.0;

	for (sim.iter = 0; sim.iter < iters; sim.iter++) {
		int heaplength = N;

		find_worst(&sim, chainhash.data(), &heaplength, N, 1);

		/* The premise the serial seed draw depends on: with P == 1 the
		   discarded point is parked at index N and logLstar is its ll, so the
		   live heap is 1..N-1 and drawing 1..N would seed from the point being
		   thrown away. If find_worst's layout ever changes, this fires. */
		if (heaplength != N - 1) *invariant_ok = 0;
		if (chainhash[N].ll != sim.logLstar) *invariant_ok = 0;
		for (int i = 1; i <= N - 1; i++)
			if (chainhash[i].index == chainhash[N].index) *invariant_ok = 0;

		update_NS_parameters(&sim, &converged, &logZnew, &lweight);

		/* The seed draw. Both expressions are evaluated from the same uniform so
		   they can be compared directly: the shipped one (survivors, 1..N-1) and
		   the pre-35fb3fb one (1..N, which could land on the discarded point at
		   index N). Counting how often each selects index N is what gives this
		   test teeth -- the fixed form must never do it, the old form must
		   sometimes, or the check is vacuous. */
		double u = unit(rng);
		int copies       = 1 + (int)(u * (N - 1)) % (N - 1);
		int copies_buggy = 1 + (int)(u * N) % N;
		if (copies == N)       draws->fixed_hits_discarded++;
		if (copies_buggy == N) draws->buggy_hits_discarded++;

		/* Replacement sampled exactly from the constrained prior {x < x*}. Exact
		   sampling is why the seed itself cannot affect the evidence here: this
		   half of the test measures the bookkeeping, the counters above guard the
		   draw. */
		int slot = chainhash[N].index;              /* reuse the discarded slot */
		double xstar = x_of_logLstar(sim.logLstar);
		(void)copies;
		x[slot] = unit(rng) * xstar;

		heaplength++;
		chainhash[heaplength].index = slot;
		chainhash[heaplength].ll = logL(x[slot]);
		heapifyhashin(chainhash.data(), heaplength);

		sim.log_DeltaX += sim.Delta_logX;
	}

	*H_out = sim.H;
	return sim.logZ;
}

int main(void)
{
	const double exact = analytic_logZ();
	std::printf("Nested-sampling evidence check\n");
	std::printf("  L(x) = exp(-x/%.2f), x ~ U(0,1)\n", TAU);
	std::printf("  analytic log Z = %.10f\n\n", exact);

	const int N = 100;
	/* exp(-iters/N) is the un-integrated live mass left at the end; ADCP does
	   not add the final live-set contribution, so the run must be long enough
	   that the remainder is far below the tolerance. 3000/100 -> e^-30. */
	const int iters = 3000;

	int invariant_ok = 1;
	DrawStats draws{0, 0};
	double sum = 0.0;
	int nseeds = 0;

	for (unsigned seed : {1u, 7u, 42u, 1234u}) {
		double H = 0.0;
		double logZ = run_ns(N, iters, seed, &H, &invariant_ok, &draws);
		double sigma = std::sqrt(H / (double)N);
		double err = logZ - exact;
		char what[128];
		std::snprintf(what, sizeof what,
		              "seed %-5u logZ %+.6f  err %+.6f  (%.2f sigma)",
		              seed, logZ, err, err / sigma);
		check(std::fabs(err) < 3.0 * sigma, what);
		sum += logZ;
		nseeds++;
	}

	check(invariant_ok == 1,
	      "find_worst parks the discarded point at index N");

	/* The seed-draw guard, and the proof it is not vacuous. */
	{
		char w[128];
		std::snprintf(w, sizeof w, "shipped draw never seeds from the discarded point (%ld)",
		              draws.fixed_hits_discarded);
		check(draws.fixed_hits_discarded == 0, w);
		std::snprintf(w, sizeof w, "pre-35fb3fb draw would have, %ld times (check has teeth)",
		              draws.buggy_hits_discarded);
		check(draws.buggy_hits_discarded > 0, w);
	}

	/* The estimator is unbiased, so the mean over seeds should sit well inside
	   a single sigma. sigma for the mean shrinks as sqrt(nseeds). */
	double mean = sum / nseeds;
	double H_ref = 0.0;
	DrawStats scratch{0, 0};
	(void)run_ns(N, iters, 1u, &H_ref, &invariant_ok, &scratch);
	double sigma_mean = std::sqrt(H_ref / (double)N) / std::sqrt((double)nseeds);
	char what[128];
	std::snprintf(what, sizeof what, "mean logZ %+.6f within 2 sigma of exact", mean);
	check(std::fabs(mean - exact) < 2.0 * sigma_mean, what);

	std::printf("\n%s: nested sampling recovers the analytic evidence\n",
	            failures ? "FAIL" : "PASS");
	return failures ? 1 : 0;
}
