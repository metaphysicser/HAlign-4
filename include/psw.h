#ifndef PSW_H_
#define PSW_H_

#include <stdint.h>

#define PSW_NEG_INF   (-0x40000000)
#define PSW_NEG_INF_F (-1e30f)

#define PSW_FLAG_SCORE_ONLY   0x01
#define PSW_FLAG_GLOBAL       0x02
#define PSW_FLAG_SEMIGLOBAL   0x04
#define PSW_FLAG_LOCAL        0x08
#define PSW_FLAG_REV_CIGAR    0x10
#define PSW_FLAG_EQX          0x20

#define PSW_CIGAR_MATCH  0
#define PSW_CIGAR_INS    1
#define PSW_CIGAR_DEL    2
#define PSW_CIGAR_N_SKIP 3
#define PSW_CIGAR_SOFTCLIP 4
#define PSW_CIGAR_HARDCLIP 5
#define PSW_CIGAR_EQ     7
#define PSW_CIGAR_X      8

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint32_t max : 31, zdropped : 1;
	int max_q, max_t;
	int score;
	int mqe, mte;
	int mqe_t, mte_q;
	int reach_end;
	int m_cigar, n_cigar;
	uint32_t *cigar;
} psw_extz_t;

typedef struct {
	int len;
	int dim;
	int depth;              /* profile 总序列数 */
	const uint32_t *prof;   /* 每列 dim 个计数；前 m 个通常是 residue/base 计数 */
} psw_prof_t;

/**
 * NW-like extension
 *
 * @param km        memory pool, when used with kalloc
 * @param q_len     query length
 * @param query     query profile
 * @param t_len     target length
 * @param target    target profile
 * @param m         number of residue types
 * @param mat       m*m scoring matrix in 1D array
 * @param q         gap open penalty
 * @param e         gap extension penalty
 * @param w         band width (<0 to disable)
 * @param zdrop     off-diagonal drop-off to stop extension (positive; <0 to disable)
 * @param flag      flags (see PSW_FLAG_* macros)
 * @param ez        (out) scores and cigar
 */
void psw_extz_pp(void *km, int q_len, const psw_prof_t *query,
              int t_len, const psw_prof_t *target,
              int8_t m, const int8_t *mat,
              int8_t q, int8_t e, int w, int zdrop, int flag, psw_extz_t *ez);

void psw_extz_ps(void *km, int qlen, const uint8_t *query,
              int t_len, const psw_prof_t *target,
              int8_t m, const int8_t *mat,
              int8_t q, int8_t e, int w, int zdrop, int flag, psw_extz_t *ez);

void psw_extz_sse_pp(void *km, int q_len, const psw_prof_t *query,
                     int t_len, const psw_prof_t *target,
                     int8_t m, const int8_t *mat,
                     int8_t q, int8_t e, int w, int zdrop, int flag, psw_extz_t *ez);

void psw_extz_sse_ps(void *km, int qlen, const uint8_t *query,
                     int t_len, const psw_prof_t *target,
                     int8_t m, const int8_t *mat,
                     int8_t q, int8_t e, int w, int zdrop, int flag, psw_extz_t *ez);

/**
 * Global alignment for profile-profile alignment
 *
 * @param km        memory pool
 * @param qlen      query length
 * @param query     query profile
 * @param tlen      target length
 * @param target    target profile
 * @param m         number of residue types
 * @param mat       m*m scoring matrix in 1D array
 * @param gapo      gap open penalty; a gap of length l costs -(gapo + l*gape)
 * @param gape      gap extension penalty
 * @param w         band width (<0 to disable)
 * @param m_cigar_  (modified) max CIGAR length; feed 0 if *cigar_==0
 * @param n_cigar_  (out) number of CIGAR elements
 * @param cigar_    (out) BAM-encoded CIGAR; caller should deallocate with kfree(km, ...)
 *
 * @return          alignment score
 */
float psw_gg_pp(void *km, int qlen, const psw_prof_t *query,
                int tlen, const psw_prof_t *target,
                int8_t m, const int8_t *mat,
                int8_t gapo, int8_t gape, int w,
                int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

float psw_gg_ps(void *km, int qlen, const uint8_t *query,
				int tlen, const psw_prof_t *target,
				int8_t m, const int8_t *mat,
				int8_t gapo, int8_t gape, int w,
				int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

float psw_gg2_pp(void *km, int qlen, const psw_prof_t *query,
                 int tlen, const psw_prof_t *target,
                 int8_t m, const int8_t *mat,
                 int8_t gapo, int8_t gape, int w,
                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

float psw_gg2_ps(void *km, int qlen, const uint8_t *query,
                 int tlen, const psw_prof_t *target,
                 int8_t m, const int8_t *mat,
                 int8_t gapo, int8_t gape, int w,
                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

float psw_gg3_pp(void *km, int qlen, const psw_prof_t *query,
                 int tlen, const psw_prof_t *target,
                 int8_t m, const int8_t *mat,
                 int8_t gapo, int8_t gape, int w,
                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

float psw_gg3_ps(void *km, int qlen, const uint8_t *query,
                 int tlen, const psw_prof_t *target,
                 int8_t m, const int8_t *mat,
                 int8_t gapo, int8_t gape, int w,
                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

float psw_gg3_sse_pp(void *km, int qlen, const psw_prof_t *query,
                     int tlen, const psw_prof_t *target,
                     int8_t m, const int8_t *mat,
                     int8_t gapo, int8_t gape, int w,
                     int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

float psw_gg3_sse_ps(void *km, int qlen, const uint8_t *query,
                     int tlen, const psw_prof_t *target,
                     int8_t m, const int8_t *mat,
                     int8_t gapo, int8_t gape, int w,
                     int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

	float psw_sw_pp(void *km, int qlen, const psw_prof_t *query,
				 int tlen, const psw_prof_t *target,
				 int8_t m, const int8_t *mat,
				 int8_t gapo, int8_t gape, int w,
				 int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

	float psw_sw_ps(void *km, int qlen, const uint8_t *query,
					 int tlen, const psw_prof_t *target,
					 int8_t m, const int8_t *mat,
					 int8_t gapo, int8_t gape, int w,
					 int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

#ifdef __cplusplus
}
#endif

/* -------------------------------------------------------------------------- */
/* allocator abstraction                                                      */
/* -------------------------------------------------------------------------- */

#ifdef HAVE_KALLOC
#include "kalloc.h"
#else
#include <stdlib.h>
#define kmalloc(km, size)          malloc((size))
#define kcalloc(km, count, size)   calloc((count), (size))
#define krealloc(km, ptr, size)    realloc((ptr), (size))
#define kfree(km, ptr)             free((ptr))
#endif

/* -------------------------------------------------------------------------- */
/* cigar helpers                                                              */
/* -------------------------------------------------------------------------- */

static inline uint32_t *psw_push_cigar(void *km, int *n_cigar, int *m_cigar,
                                       uint32_t *cigar, uint32_t op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1] & 0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar ? (*m_cigar) << 1 : 4;
			cigar = (uint32_t*)krealloc(km, cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = (len << 4) | op;
	} else {
		cigar[(*n_cigar) - 1] += (len << 4);
	}
	return cigar;
}

/*
 * In the backtrack matrix, value p[] has the following structure:
 *   bit 0-2: which type gets the max
 *            0 for H, 1 for E, 2 for F, 3 for \tilde{E}, 4 for \tilde{F}
 *   bit 3/0x08: 1 if continuation on E   (bit 5/0x20 for \tilde{E})
 *   bit 4/0x10: 1 if continuation on F   (bit 6/0x40 for \tilde{F})
 */
static inline void psw_backtrack(void *km, int is_rot, int is_rev, int min_intron_len,
                                 const uint8_t *p, const int *off, const int *off_end,
                                 int n_col, int i0, int j0,
                                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
	uint32_t *cigar = *cigar_, tmp;

	while (i >= 0 && j >= 0) {
		int force_state = -1;
		if (is_rot) {
			r = i + j;
			if (i < off[r]) force_state = 2;
			if (off_end && i > off_end[r]) force_state = 1;
			tmp = force_state < 0 ? p[(size_t)r * n_col + i - off[r]] : 0;
		} else {
			if (j < off[i]) force_state = 2;
			if (off_end && j > off_end[i]) force_state = 1;
			tmp = force_state < 0 ? p[(size_t)i * n_col + j - off[i]] : 0;
		}

		if (state == 0) state = tmp & 7;
		else if (!((tmp >> (state + 2)) & 1)) state = 0;
		if (state == 0) state = tmp & 7;
		if (force_state >= 0) state = force_state;

		if (state == 0) {
			cigar = psw_push_cigar(km, &n_cigar, &m_cigar, cigar, PSW_CIGAR_MATCH, 1);
			--i; --j;
		} else if (state == 1 || (state == 3 && min_intron_len <= 0)) {
			cigar = psw_push_cigar(km, &n_cigar, &m_cigar, cigar, PSW_CIGAR_DEL, 1);
			--i;
		} else if (state == 3 && min_intron_len > 0) {
			cigar = psw_push_cigar(km, &n_cigar, &m_cigar, cigar, PSW_CIGAR_N_SKIP, 1);
			--i;
		} else {
			cigar = psw_push_cigar(km, &n_cigar, &m_cigar, cigar, PSW_CIGAR_INS, 1);
			--j;
		}
	}

	if (i >= 0) {
		cigar = psw_push_cigar(km, &n_cigar, &m_cigar, cigar,
		                       min_intron_len > 0 && i >= min_intron_len ? PSW_CIGAR_N_SKIP : PSW_CIGAR_DEL,
		                       i + 1);
	}
	if (j >= 0) {
		cigar = psw_push_cigar(km, &n_cigar, &m_cigar, cigar, PSW_CIGAR_INS, j + 1);
	}

	if (!is_rev) {
		for (i = 0; i < (n_cigar >> 1); ++i) {
			tmp = cigar[i];
			cigar[i] = cigar[n_cigar - 1 - i];
			cigar[n_cigar - 1 - i] = tmp;
		}
	}

	*m_cigar_ = m_cigar;
	*n_cigar_ = n_cigar;
	*cigar_ = cigar;
}

/* convert MATCH cigar to EQ/X cigar */
static inline void psw_cigar2eqx(void *km,
                                 const uint8_t *target, const uint8_t *query,
                                 int nc0, const uint32_t *ci0,
                                 int *mc1, int *nc1, uint32_t **ci1)
{
	int i, k, x = 0, y = 0;
	*nc1 = 0;

	for (k = 0; k < nc0; ++k) {
		int op = ci0[k] & 0xf;
		int len = ci0[k] >> 4;

		if (op == PSW_CIGAR_MATCH) {
			for (i = 0; i < len; ++i) {
				*ci1 = psw_push_cigar(km, nc1, mc1, *ci1,
				                      target[x + i] == query[y + i] ? PSW_CIGAR_EQ : PSW_CIGAR_X,
				                      1);
			}
			x += len;
			y += len;
		} else {
			*ci1 = psw_push_cigar(km, nc1, mc1, *ci1, op, len);
			if (op == PSW_CIGAR_DEL || op == PSW_CIGAR_N_SKIP) x += len;
			else if (op == PSW_CIGAR_INS) y += len;
			else if (op == PSW_CIGAR_EQ || op == PSW_CIGAR_X) {
				x += len;
				y += len;
			}
		}
	}
}

static inline void psw_reset_extz(psw_extz_t *ez)
{
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0;
	ez->score = ez->mqe = ez->mte = PSW_NEG_INF;
	ez->n_cigar = 0;
	ez->m_cigar = 0;
	ez->zdropped = 0;
	ez->reach_end = 0;
	ez->cigar = 0;
}

static inline int psw_apply_zdrop(psw_extz_t *ez, int is_rot, int32_t H,
                                  int a, int b, int zdrop, int8_t e)
{
	int r, t;
	if (is_rot) r = a, t = b;
	else r = a + b, t = a;

	if (H > (int32_t)ez->max) {
		ez->max = H;
		ez->max_t = t;
		ez->max_q = r - t;
	} else if (t >= ez->max_t && r - t >= ez->max_q) {
		int tl = t - ez->max_t;
		int ql = (r - t) - ez->max_q;
		int l = tl > ql ? tl - ql : ql - tl;
		if (zdrop >= 0 && ez->max - H > (uint32_t)(zdrop + l * e)) {
			ez->zdropped = 1;
			return 1;
		}
	}
	return 0;
}

#endif /* PSW_H_ */

