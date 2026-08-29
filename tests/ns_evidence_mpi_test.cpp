/* =============================================================================
 * Nested-sampling evidence check, driven across real MPI ranks.
 *
 * tests/ns_evidence_test.cpp already checks find_worst()/update_NS_parameters()
 * against a closed-form evidence, but it hardcodes P=1 and never touches a
 * network. This is the P>1 companion: it drives the actual cross-rank
 * functions that "Phase 2 progress -- MPI unblocked" in MIGRATION.md newly
 * makes real -- collect_chains(), return_and_reheap_chains(), and through them
 * mpi_send_chain()/mpi_rec_chain() -- under a real mpirun, against the same
 * problem with a known answer:
 *
 *     prior       x ~ Uniform(0,1)
 *     likelihood  L(x) = exp(-x/tau),  so  log L = -x/tau
 *     evidence    Z = tau * (1 - exp(-1/tau))
 *
 * As in the serial test, L is monotone in x, so a constrained sample is
 * Uniform(0, x*) drawn exactly -- no MCMC, no peptide, no energy function.
 * That is what makes it possible to build valid Chain objects without going
 * through build_peptide_from_sequence/biasmap_initialise/aat_init at all: a
 * Chain allocated by allocmem_chain(&c, 1, 1) with only .ll set is everything
 * collect_chains/return_and_reheap_chains/mpi_send_chain/mpi_rec_chain touch
 * for this problem, since none of them interpret Chain's geometry fields --
 * they just move whatever is there.
 *
 * The initial population is built identically on every rank from a shared
 * seed (same std::mt19937 sequence run on every process), then each rank
 * keeps only the points store_chain's real round-robin scheme would have
 * given it: point i (1-indexed) belongs to processor (i-1) % P, at local
 * index (i-1) / P. That sidesteps needing an initial-population broadcast of
 * our own, while still handing collect_chains/return_and_reheap_chains
 * exactly the chainhash/cpoints shape nestedsampling() itself builds.
 *
 * If any of the newly-unblocked MPI plumbing silently drops or misroutes a
 * point across ranks, the evidence estimate is the thing that would show it:
 * find_worst's boundary would be wrong on some ranks, or a stale/duplicated
 * .ll would leak in, and the estimate would drift well outside the standard
 * nested-sampling error band.
 * ========================================================================== */
#include <mpi.h>

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

/* Defined in nested.cpp, PARALLEL-only. Not declared in nested.h -- only
   nestedsampling() is public -- but they have external linkage. */
void find_worst(simulation_params *sim_params, ChainHash *chainhash, int *heaplength, int N, int P);
void update_NS_parameters(simulation_params *sim_params, int *converged, double *logZnew, double *lweight);
void collect_chains(ChainHash *chainhash, Chain *cpoints, Chain *chaincopies, simulation_params *sim_params,
                     int rank, int P, int N, MPI_Comm *NS_WORLD);
void return_and_reheap_chains(ChainHash *chainhash, Chain *cpoints, Chain *chaincopies, simulation_params *sim_params,
                               int rank, int P, int *heaplength, MPI_Comm *NS_WORLD);

static const double TAU = 0.2;
static const int NAA = 1, NCHAINS = 1; /* geometry is irrelevant here -- see file header */

static double logL(double x) { return -x / TAU; }
static double x_of_logLstar(double ll) { return -TAU * ll; }
static double analytic_logZ() { return std::log(TAU) + std::log1p(-std::exp(-1.0 / TAU)); }

/* One NS run across all P ranks. Returns this rank's view: only rank 0's
   *logZ_out/*H_out are meaningful -- update_NS_parameters only ever runs
   there, exactly as in nestedsampling(). */
static void run_ns_mpi(int N, int iters, unsigned seed, int rank, int P,
                        double *logZ_out, double *H_out, int *invariant_failures)
{
	MPI_Comm ns_world = MPI_COMM_WORLD;

	/* Same sequence on every rank -- see file header. Only used to build the
	   values each rank needs for the indices it owns. */
	std::mt19937 shared_rng(seed);
	std::uniform_real_distribution<double> unit(0.0, 1.0);
	std::vector<double> x(N + 1);
	for (int i = 1; i <= N; i++) x[i] = unit(shared_rng);

	std::vector<ChainHash> chainhash(N + 2);
	if (rank == 0) {
		for (int i = 1; i <= N; i++) {
			chainhash[i].processor = (i - 1) % P;
			chainhash[i].index = (i - 1) / P;
			chainhash[i].ll = logL(x[i]);
		}
		constructhashheap(chainhash.data(), N);
	}

	std::vector<Chain> cpoints;
	int current_stored = 0;
	for (int i = 1; i <= N; i++) {
		if ((i - 1) % P == rank) {
			cpoints.resize(current_stored + 1);
			cpoints[current_stored].NAA = NAA;
			cpoints[current_stored].Nchains = NCHAINS;
			allocmem_chain(&cpoints[current_stored], NAA, NCHAINS);
			cpoints[current_stored].ll = logL(x[i]);
			current_stored++;
		}
	}

	int size_of_chaincopies = (rank == 0) ? P : 1;
	std::vector<Chain> chaincopies(size_of_chaincopies);
	for (int k = 0; k < size_of_chaincopies; k++) {
		chaincopies[k].NAA = NAA;
		chaincopies[k].Nchains = NCHAINS;
		allocmem_chain(&chaincopies[k], NAA, NCHAINS);
	}

	simulation_params sim{};
	sim.thermobeta = 1.0;
	sim.lowtemp = 0;
	sim.checkpoint = 0;
	sim.iter_start = 0;
	sim.N = N;
	sim.logX_start = 0.0;
	sim.logZ = -std::numeric_limits<double>::max();
	sim.H = 0.0;
	/* nestedsampling()'s own PARALLEL setup -- a linear approximation to
	   exp(-1/N), not the serial exp(-1/N) itself. */
	sim.alpha = 1.0 - (double)P / (double)(N + 1);
	sim.Delta_logX = std::log(sim.alpha);
	sim.log_DeltaX = std::log(1.0 - sim.alpha);

	std::mt19937 draw_rng(seed * 1000u + (unsigned)rank + 1u);
	int converged = 1;

	for (sim.iter = 0; sim.iter < iters && converged; sim.iter++) {
		int heaplength = N;
		if (rank == 0) {
			double logZnew = 0.0, lweight = 0.0;
			find_worst(&sim, chainhash.data(), &heaplength, N, P);
			update_NS_parameters(&sim, &converged, &logZnew, &lweight);
		}

		/* The function under test: routes the P discarded points'
		   replacements across ranks via mpi_send_chain/mpi_rec_chain, and
		   (from==0 there) is also how every non-zero rank learns this
		   iteration's logLstar. */
		collect_chains(chainhash.data(), cpoints.data(), chaincopies.data(), &sim, rank, P, N, &ns_world);

		/* Exact constrained draw, replacing move()'s MCMC -- see file header.
		   Every rank does its own, exactly as nestedsampling() does. */
		double xstar = x_of_logLstar(sim.logLstar);
		double new_x = unit(draw_rng) * xstar;
		chaincopies[0].ll = logL(new_x);
		if (!(chaincopies[0].ll > sim.logLstar)) (*invariant_failures)++;

		return_and_reheap_chains(chainhash.data(), cpoints.data(), chaincopies.data(), &sim, rank, P, &heaplength, &ns_world);

		if (rank == 0) sim.log_DeltaX += sim.Delta_logX;
		MPI_Bcast(&converged, 1, MPI_INT, 0, ns_world);
	}

	for (auto &c : cpoints) freemem_chain(&c);
	for (auto &c : chaincopies) freemem_chain(&c);

	*logZ_out = sim.logZ;
	*H_out = sim.H;
}

int main(void)
{
	MPI_Init(NULL, NULL);
	int rank = 0, P = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &P);

	int failures = 0;

	if (rank == 0) {
		std::printf("Nested-sampling evidence check (MPI, P=%d)\n", P);
		std::printf("  L(x) = exp(-x/%.2f), x ~ U(0,1)\n", TAU);
		std::printf("  analytic log Z = %.10f\n\n", analytic_logZ());
	}

	const double exact = analytic_logZ();
	const int N = 100;
	/* Matches the serial test's total shrinkage depth (~3000 steps) however
	   many ranks split the work: each iteration here removes P points, not 1. */
	const int iters = (3000 + P - 1) / P;

	double sum = 0.0;
	int nseeds = 0;
	for (unsigned seed : {11u, 2027u}) {
		if (rank == 0) srand(seed);
		double logZ = 0.0, H = 0.0;
		int invariant_failures = 0;
		run_ns_mpi(N, iters, seed, rank, P, &logZ, &H, &invariant_failures);

		int total_invariant_failures = 0;
		MPI_Reduce(&invariant_failures, &total_invariant_failures, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			double sigma = std::sqrt(H / (double)N);
			double err = logZ - exact;
			bool ok = std::fabs(err) < 3.0 * sigma;
			std::printf("  seed %-5u logZ %+.6f  err %+.6f  (%.2f sigma)  %s\n",
			            seed, logZ, err, err / sigma, ok ? "ok" : "FAIL");
			if (!ok) failures++;
			if (total_invariant_failures != 0) {
				std::printf("  seed %-5u %d draws violated ll > logLstar (FAIL)\n",
				            seed, total_invariant_failures);
				failures++;
			}
			sum += logZ;
			nseeds++;
		}
	}

	if (rank == 0) {
		double mean = sum / nseeds;
		bool ok = std::fabs(mean - exact) < 1.0; /* loose: 2 seeds only, this just catches gross breakage */
		std::printf("  mean logZ %+.6f  %s\n", mean, ok ? "ok" : "FAIL");
		if (!ok) failures++;

		std::printf("\n%s: nested sampling (MPI, P=%d) recovers the analytic evidence\n",
		            failures ? "FAIL" : "PASS", P);
	}

	int global_failures = 0;
	MPI_Bcast(&failures, 1, MPI_INT, 0, MPI_COMM_WORLD);
	global_failures = failures;

	MPI_Finalize();
	return global_failures ? 1 : 0;
}
