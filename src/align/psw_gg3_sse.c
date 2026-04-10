#include <limits.h>
#include <stdint.h>
#include <emmintrin.h>
#include "psw.h"

#ifndef PSW_NEG_INF_F
#define PSW_NEG_INF_F (-1e30f)
#endif

#ifndef PSW_GG3_SCALE
#define PSW_GG3_SCALE 100
#endif

typedef struct { int16_t u, v, x, y; } uvxy16_t;

static inline int16_t psw_sat16(int32_t x)
{
	if (x > INT16_MAX) return INT16_MAX;
	if (x < INT16_MIN) return INT16_MIN;
	return (int16_t)x;
}

static inline int16_t psw_pick_safe_scale_pow2(int8_t m, const int8_t *mat, int8_t gapo, int8_t gape)
{
	int a, b;
	int32_t smax = 0;
	int32_t coeff, max_scale;
	int32_t pow2 = 1;

	for (a = 0; a < m; ++a)
		for (b = 0; b < m; ++b)
			if ((int32_t)mat[a * m + b] > smax)
				smax = (int32_t)mat[a * m + b];

	coeff = smax + (int32_t)2 * gapo + (int32_t)2 * gape;
	if (coeff <= 0) coeff = 1;
	max_scale = INT16_MAX / coeff;
	while ((pow2 << 1) > 0 && (pow2 << 1) <= max_scale)
		pow2 <<= 1;
	if (pow2 < 2) return 0;
	return (int16_t)pow2;
}

static inline int8_t psw_scale_to_shift(int16_t scale)
{
	int8_t sh = 0;
	if (scale <= 0) return -1;
	while ((scale & 1) == 0) {
		scale >>= 1;
		++sh;
	}
	return scale == 1 ? sh : -1;
}

static inline int16_t psw_to_scaled_ratio(int32_t num, int32_t den, int16_t scale)
{
	int64_t den64, scaled, rounded;

	if (den <= 0) den = 1;
	if (scale <= 0) scale = 1;
	den64 = (int64_t)den;
	scaled = (int64_t)num * (int64_t)scale;
	if (scaled >= 0) rounded = (scaled + den64 / 2) / den64;
	else             rounded = (scaled - den64 / 2) / den64;

	if (rounded > INT16_MAX) return INT16_MAX;
	if (rounded < INT16_MIN) return INT16_MIN;
	return (int16_t)rounded;
}

static inline int psw_make_norm_prof(void *km, const psw_prof_t *src, psw_prof_t *dst, int8_t m, int16_t scale)
{
	int i, a;
	int32_t depth;
	uint32_t *prof;
	if (src == 0 || dst == 0 || src->prof == 0 || src->len < 0 || src->dim <= 0 || m <= 0 || scale <= 0) return 0;
	if (src->dim < m) return 0;

	prof = (uint32_t*)kmalloc(km, (size_t)src->len * src->dim * sizeof(uint32_t));
	if (prof == 0) return 0;

	depth = src->depth > 0 ? (int32_t)src->depth : 1;
	for (i = 0; i < src->len; ++i) {
		const uint32_t *in = src->prof + (size_t)i * src->dim;
		uint32_t *out = prof + (size_t)i * src->dim;
		for (a = 0; a < src->dim; ++a) {
			if (a < m) {
				int16_t v = psw_to_scaled_ratio((int32_t)in[a], depth, scale);
				if (v < 0) v = 0;
				if (v > scale) v = scale;
				out[a] = (uint32_t)v;
			} else out[a] = in[a];
		}
	}

	dst->len = src->len;
	dst->dim = src->dim;
	dst->depth = scale;
	dst->prof = prof;
	return 1;
}

static inline void psw_free_norm_prof(void *km, psw_prof_t *p)
{
	if (p && p->prof) {
		kfree(km, (void*)p->prof);
		p->prof = 0;
	}
}

static inline int16_t psw_dot_scaled_sse_m5(const int16_t *x, const int16_t *y, int8_t scale_shift)
{
	__m128i xv = _mm_loadl_epi64((const __m128i*)x);
	__m128i yv = _mm_loadl_epi64((const __m128i*)y);
	__m128i madd = _mm_madd_epi16(xv, yv);

	// Move lane1 to lane0 then accumulate.
	__m128i hi = _mm_srli_si128(madd, 4);
	madd = _mm_add_epi32(madd, hi);

	int32_t acc = _mm_cvtsi128_si32(madd);
	acc += (int32_t)x[4] * y[4];
	acc >>= scale_shift;
	return psw_sat16(acc);
}

static inline int32_t psw_hsum_epi32_sse2(__m128i v)
{
	__m128i hi = _mm_srli_si128(v, 8);
	__m128i s = _mm_add_epi32(v, hi);
	hi = _mm_srli_si128(s, 4);
	s = _mm_add_epi32(s, hi);
	return _mm_cvtsi128_si32(s);
}

static inline int32_t psw_dot_acc_sse_m20(const int16_t *x, const int16_t *y)
{
	__m128i s = _mm_setzero_si128();
	__m128i xv, yv;

	xv = _mm_loadu_si128((const __m128i*)x);
	yv = _mm_loadu_si128((const __m128i*)y);
	s = _mm_add_epi32(s, _mm_madd_epi16(xv, yv));

	xv = _mm_loadu_si128((const __m128i*)(x + 8));
	yv = _mm_loadu_si128((const __m128i*)(y + 8));
	s = _mm_add_epi32(s, _mm_madd_epi16(xv, yv));

	xv = _mm_loadl_epi64((const __m128i*)(x + 16));
	yv = _mm_loadl_epi64((const __m128i*)(y + 16));
	s = _mm_add_epi32(s, _mm_madd_epi16(xv, yv));

	return psw_hsum_epi32_sse2(s);
}

static inline int16_t psw_dot_scaled_sse_m20(const int16_t *x, const int16_t *y, int8_t scale_shift)
{
	int32_t acc = psw_dot_acc_sse_m20(x, y);
	acc >>= scale_shift;
	return psw_sat16(acc);
}

static inline int16_t psw_dot_scaled_sse_m21(const int16_t *x, const int16_t *y, int8_t scale_shift)
{
	int32_t acc = psw_dot_acc_sse_m20(x, y) + (int32_t)x[20] * y[20];
	acc >>= scale_shift;
	return psw_sat16(acc);
}

static inline int16_t psw_dot_scaled_sse(const int16_t *x, const int16_t *y, int m, int8_t scale_shift)
{
	__m128i sum = _mm_setzero_si128();
	int i = 0;
	int32_t acc;
	int32_t tmp[4];

	if (scale_shift <= 0) return 0;
	if (m == 5) return psw_dot_scaled_sse_m5(x, y, scale_shift);
	if (m == 20) return psw_dot_scaled_sse_m20(x, y, scale_shift);
	if (m == 21) return psw_dot_scaled_sse_m21(x, y, scale_shift);
	for (; i + 8 <= m; i += 8) {
		__m128i xv = _mm_loadu_si128((const __m128i*)(x + i));
		__m128i yv = _mm_loadu_si128((const __m128i*)(y + i));
		sum = _mm_add_epi32(sum, _mm_madd_epi16(xv, yv));
	}
	_mm_storeu_si128((__m128i*)tmp, sum);
	acc = tmp[0] + tmp[1] + tmp[2] + tmp[3];
	for (; i < m; ++i)
		acc += (int32_t)x[i] * y[i];
	acc >>= scale_shift;
	return psw_sat16(acc);
}

static inline int16_t *psw_gen_base_freq_i16(void *km, int len, const psw_prof_t *p, int8_t m, int16_t scale)
{
	int i, a;
	int16_t *bf;

	bf = (int16_t*)kmalloc(km, (size_t)len * sizeof(int16_t));
	if (bf == 0) return 0;

	for (i = 0; i < len; ++i) {
		const uint32_t *col = p->prof + (size_t)i * p->dim;
		int32_t sum = 0;
		for (a = 0; a < m; ++a)
			sum += (int32_t)col[a];
		bf[i] = psw_sat16(sum);
		if (bf[i] < 0) bf[i] = 0;
		if (bf[i] > scale) bf[i] = scale;
	}
	return bf;
}

static inline int16_t *psw_gen_qp_i16(void *km, int qlen, const psw_prof_t *query, int8_t m, const int8_t *mat)
{
	int a, b, j;
	int16_t *qp;

	qp = (int16_t*)kmalloc(km, (size_t)qlen * m * sizeof(int16_t));
	if (qp == 0) return 0;

	for (j = 0; j < qlen; ++j) {
		const uint32_t *qcol = query->prof + (size_t)j * query->dim;
		int16_t *dst = qp + (size_t)j * m;
		for (b = 0; b < m; ++b) {
			int32_t s = 0;
			for (a = 0; a < m; ++a)
				s += (int32_t)qcol[a] * (int32_t)mat[a * m + b];
			dst[b] = psw_sat16(s);
		}
	}
	return qp;
}

static inline int16_t *psw_gen_tf_i16(void *km, int tlen, const psw_prof_t *target, int8_t m)
{
	int i, b;
	int16_t *tf = (int16_t*)kmalloc(km, (size_t)tlen * m * sizeof(int16_t));
	if (tf == 0) return 0;

	for (i = 0; i < tlen; ++i) {
		const uint32_t *tcol = target->prof + (size_t)i * target->dim;
		int16_t *dst = tf + (size_t)i * m;
		for (b = 0; b < m; ++b)
			dst[b] = psw_sat16((int32_t)tcol[b]);
	}
	return tf;
}

static inline int16_t *psw_gen_tp_i16(void *km, int tlen, const psw_prof_t *target,
                                      int8_t m, const int8_t *mat)
{
	int a, b, i;
	int16_t *tp;

	tp = (int16_t*)kmalloc(km, (size_t)m * tlen * sizeof(int16_t));
	if (tp == 0) return 0;

	for (i = 0; i < tlen; ++i) {
		const uint32_t *tcol = target->prof + (size_t)i * target->dim;
		for (a = 0; a < m; ++a) {
			int32_t s = 0;
			for (b = 0; b < m; ++b)
				s += (int32_t)mat[a * m + b] * (int32_t)tcol[b];
			tp[(size_t)a * tlen + i] = psw_sat16(s);
		}
	}
	return tp;
}

static inline int32_t psw_gap_only_target_i16(const int16_t *go_t, const int16_t *ge_t, int tlen)
{
	int i;
	int32_t s = 0;
	for (i = 0; i < tlen; ++i)
		s -= (i == 0 ? (int32_t)go_t[i] + ge_t[i] : (int32_t)ge_t[i]);
	return s;
}

static inline int32_t psw_gap_only_query_i16(const int16_t *go_q, const int16_t *ge_q, int qlen)
{
	int j;
	int32_t s = 0;
	for (j = 0; j < qlen; ++j)
		s -= (j == 0 ? (int32_t)go_q[j] + ge_q[j] : (int32_t)ge_q[j]);
	return s;
}

static inline int32_t psw_gap_only_query_scalar_i16(int16_t go_q, int16_t ge_q, int qlen)
{
	if (qlen <= 0) return 0;
	return -((int32_t)go_q + (int32_t)qlen * ge_q);
}

#ifdef __SSE2__

typedef struct {
	__m128i *u, *v, *x, *y, *s;
	int16_t *u16, *v16, *x16, *y16;
	int n_vec;
	uint8_t *mem;
} psw_sse_state_t;

static inline int psw_sse_state_init(void *km, int tlen, psw_sse_state_t *st)
{
	size_t n_vec = (size_t)((tlen + 7) / 8);
	uint8_t *mem;
	if (st == 0) return 0;
	mem = (uint8_t*)kcalloc(km, n_vec * 5 + 1, 16);
	if (mem == 0) return 0;
	st->mem = mem;
	st->u = (__m128i*)(((size_t)mem + 15) >> 4 << 4);
	st->v = st->u + n_vec;
	st->x = st->v + n_vec;
	st->y = st->x + n_vec;
	st->s = st->y + n_vec;
	st->u16 = (int16_t*)st->u;
	st->v16 = (int16_t*)st->v;
	st->x16 = (int16_t*)st->x;
	st->y16 = (int16_t*)st->y;
	st->n_vec = (int)n_vec;
	return 1;
}

static inline void psw_sse_state_destroy(void *km, psw_sse_state_t *st)
{
	if (st && st->mem) kfree(km, st->mem);
}

static inline __m128i psw_set_low_i16(int16_t x)
{
	return _mm_cvtsi32_si128((int)(uint16_t)x);
}

static inline void psw_store_flags8(uint8_t *dst, __m128i d)
{
	int16_t tmp[8];
	int i;
	_mm_storeu_si128((__m128i*)tmp, d);
	for (i = 0; i < 8; ++i) dst[i] = (uint8_t)tmp[i];
}

static inline void psw_store_flags8_band(uint8_t *dst, int tbase, int st0, int en0, __m128i d)
{
	int16_t tmp[8];
	int i;
	_mm_storeu_si128((__m128i*)tmp, d);
	if (tbase >= st0 && tbase + 7 <= en0) {
		uint8_t *p = dst + (tbase - st0);
		for (i = 0; i < 8; ++i) p[i] = (uint8_t)tmp[i];
		return;
	}
	for (i = 0; i < 8; ++i) {
		int t = tbase + i;
		if (t >= st0 && t <= en0)
			dst[t - st0] = (uint8_t)tmp[i];
	}
}

static inline __m128i psw_load_band_go_t(const int16_t *go_t, int tbase, int st0, int en0, int tlen)
{
	if (tbase >= st0 && tbase + 7 <= en0 && tbase + 7 < tlen)
		return _mm_loadu_si128((const __m128i*)(go_t + tbase));
	else {
		int16_t tmp[8];
		int i;
		for (i = 0; i < 8; ++i) {
			int t = tbase + i;
			tmp[i] = (t >= st0 && t <= en0 && t < tlen) ? go_t[t] : 0;
		}
		return _mm_loadu_si128((const __m128i*)tmp);
	}
}

// m-specific helpers for SSE hot loop, avoid per-lane m branching
static inline void psw_fill_pp_score_block_m5(int16_t *dst, int tbase, int r, int st0, int en0,
                                              int qlen, int tlen, int8_t scale_shift,
                                              const int16_t *qp, const int16_t *tf,
                                              const int16_t *go_ge_q, const int16_t *go_ge_t)
{
    int lane;
    (void)qlen;

    if (tbase >= st0 && tbase + 7 <= en0 && tbase + 7 < tlen) {
        int t = tbase;
        int j = r - tbase;
        const int16_t *qpj = qp + (size_t)j * 5;
        const int16_t *tfi = tf + (size_t)t * 5;
        for (lane = 0; lane < 8; ++lane, ++t, --j, qpj -= 5, tfi += 5)
            dst[lane] = (int16_t)(psw_dot_scaled_sse_m5(qpj, tfi, scale_shift) + go_ge_t[t] + go_ge_q[j]);
        return;
    }

    for (lane = 0; lane < 8; ++lane) {
        int t = tbase + lane;
        if (t >= st0 && t <= en0 && t < tlen) {
            int j = r - t;
            const int16_t *qpj = qp + (size_t)j * 5;
            const int16_t *tfi = tf + (size_t)t * 5;
            dst[lane] = (int16_t)(psw_dot_scaled_sse_m5(qpj, tfi, scale_shift) + go_ge_t[t] + go_ge_q[j]);
        } else dst[lane] = 0;
    }
}

static inline void psw_fill_pp_score_block_m20(int16_t *dst, int tbase, int r, int st0, int en0,
                                               int qlen, int tlen, int8_t scale_shift,
                                               const int16_t *qp, const int16_t *tf,
                                               const int16_t *go_ge_q, const int16_t *go_ge_t)
{
    int lane;
    (void)qlen;

    if (tbase >= st0 && tbase + 7 <= en0 && tbase + 7 < tlen) {
        int t = tbase;
        int j = r - tbase;
        const int16_t *qpj = qp + (size_t)j * 20;
        const int16_t *tfi = tf + (size_t)t * 20;
        for (lane = 0; lane < 8; ++lane, ++t, --j, qpj -= 20, tfi += 20)
            dst[lane] = (int16_t)(psw_dot_scaled_sse_m20(qpj, tfi, scale_shift) + go_ge_t[t] + go_ge_q[j]);
        return;
    }

    for (lane = 0; lane < 8; ++lane) {
        int t = tbase + lane;
        if (t >= st0 && t <= en0 && t < tlen) {
            int j = r - t;
            const int16_t *qpj = qp + (size_t)j * 20;
            const int16_t *tfi = tf + (size_t)t * 20;
            dst[lane] = (int16_t)(psw_dot_scaled_sse_m20(qpj, tfi, scale_shift) + go_ge_t[t] + go_ge_q[j]);
        } else dst[lane] = 0;
    }
}

static inline void psw_fill_pp_score_block_m21(int16_t *dst, int tbase, int r, int st0, int en0,
                                               int qlen, int tlen, int8_t scale_shift,
                                               const int16_t *qp, const int16_t *tf,
                                               const int16_t *go_ge_q, const int16_t *go_ge_t)
{
    int lane;
    (void)qlen;

    if (tbase >= st0 && tbase + 7 <= en0 && tbase + 7 < tlen) {
        int t = tbase;
        int j = r - tbase;
        const int16_t *qpj = qp + (size_t)j * 21;
        const int16_t *tfi = tf + (size_t)t * 21;
        for (lane = 0; lane < 8; ++lane, ++t, --j, qpj -= 21, tfi += 21)
            dst[lane] = (int16_t)(psw_dot_scaled_sse_m21(qpj, tfi, scale_shift) + go_ge_t[t] + go_ge_q[j]);
        return;
    }

    for (lane = 0; lane < 8; ++lane) {
        int t = tbase + lane;
        if (t >= st0 && t <= en0 && t < tlen) {
            int j = r - t;
            const int16_t *qpj = qp + (size_t)j * 21;
            const int16_t *tfi = tf + (size_t)t * 21;
            dst[lane] = (int16_t)(psw_dot_scaled_sse_m21(qpj, tfi, scale_shift) + go_ge_t[t] + go_ge_q[j]);
        } else dst[lane] = 0;
    }
}

static inline void psw_fill_pp_score_block_generic(int16_t *dst, int tbase, int r, int st0, int en0,
                                                   int qlen, int tlen, int m, int8_t scale_shift,
                                                   const int16_t *qp, const int16_t *tf,
                                                   const int16_t *go_ge_q, const int16_t *go_ge_t)
{
    int lane;
    (void)qlen;

    if (tbase >= st0 && tbase + 7 <= en0 && tbase + 7 < tlen) {
        int t = tbase;
        int j = r - tbase;
        const int16_t *qpj = qp + (size_t)j * m;
        const int16_t *tfi = tf + (size_t)t * m;
        for (lane = 0; lane < 8; ++lane, ++t, --j, qpj -= m, tfi += m)
            dst[lane] = (int16_t)(psw_dot_scaled_sse(qpj, tfi, m, scale_shift) + go_ge_t[t] + go_ge_q[j]);
        return;
    }

    for (lane = 0; lane < 8; ++lane) {
        int t = tbase + lane;
        if (t >= st0 && t <= en0 && t < tlen) {
            int j = r - t;
            const int16_t *qpj = qp + (size_t)j * m;
            const int16_t *tfi = tf + (size_t)t * m;
            dst[lane] = (int16_t)(psw_dot_scaled_sse(qpj, tfi, m, scale_shift) + go_ge_t[t] + go_ge_q[j]);
        } else dst[lane] = 0;
    }
}

static inline void psw_fill_pp_score_block(int16_t *dst, int tbase, int r, int st0, int en0,
                                           int qlen, int tlen, int m, int8_t scale_shift,
                                           const int16_t *qp, const int16_t *tf,
                                           const int16_t *go_ge_q, const int16_t *go_ge_t)
{
    if (m == 20) {
        psw_fill_pp_score_block_m20(dst, tbase, r, st0, en0, qlen, tlen, scale_shift, qp, tf, go_ge_q, go_ge_t);
        return;
    }
    if (m == 21) {
        psw_fill_pp_score_block_m21(dst, tbase, r, st0, en0, qlen, tlen, scale_shift, qp, tf, go_ge_q, go_ge_t);
        return;
    }
    if (m == 5) {
        psw_fill_pp_score_block_m5(dst, tbase, r, st0, en0, qlen, tlen, scale_shift, qp, tf, go_ge_q, go_ge_t);
        return;
    }
    psw_fill_pp_score_block_generic(dst, tbase, r, st0, en0, qlen, tlen, m, scale_shift, qp, tf, go_ge_q, go_ge_t);
}

static inline void psw_fill_ps_score_block(int16_t *dst, int tbase, int r, int st0, int en0,
                                           int qlen, int tlen, const uint8_t *query,
                                           const int16_t *tp, int16_t go_q, int16_t ge_q,
                                           const int16_t *go_t, const int16_t *ge_t)
{
    int lane;
    const int16_t go_qe = (int16_t)(go_q + ge_q);
    (void)qlen;

    if (tbase >= st0 && tbase + 7 <= en0 && tbase + 7 < tlen) {
        int t = tbase;
        int j = r - tbase;
        for (lane = 0; lane < 8; ++lane, ++t, --j) {
            int aidx = (int)query[j];
            dst[lane] = (int16_t)(tp[(size_t)aidx * tlen + t] + go_t[t] + ge_t[t] + go_qe);
        }
        return;
    }

    for (lane = 0; lane < 8; ++lane) {
        int t = tbase + lane;
        if (t >= st0 && t <= en0 && t < tlen) {
            int j = r - t;
            int aidx = (int)query[j];
            dst[lane] = (int16_t)(tp[(size_t)aidx * tlen + t] + go_t[t] + ge_t[t] + go_qe);
        } else dst[lane] = 0;
    }
}

static inline void psw_fill_rev_qgap_block(int16_t *dst, int tbase, int r, int st0, int en0,
                                           int qlen, const int16_t *go_q)
{
    int lane;
    (void)qlen;
    if (tbase >= st0 && tbase + 7 <= en0) {
        const int16_t *g = go_q + (r - tbase);
        for (lane = 0; lane < 8; ++lane) dst[lane] = g[-lane];
        return;
    }
    for (lane = 0; lane < 8; ++lane) {
        int t = tbase + lane;
        if (t >= st0 && t <= en0) {
            int j = r - t;
            dst[lane] = go_q[j];
        } else dst[lane] = 0;
    }
}

static inline void psw_sse_core_pp(int r, int st0, int en0, int st, int en,
                                   psw_sse_state_t *sv, uint8_t *zr,
                                   const int16_t *go_t, const int16_t *go_q,
                                   int qlen, int tlen, int m, int8_t scale_shift,
                                   const int16_t *qp, const int16_t *tf,
                                   const int16_t *go_ge_q, const int16_t *go_ge_t,
                                   int16_t x1, int16_t v1)
{
	int t;
	const __m128i zero_ = _mm_setzero_si128();
	const __m128i flag1_ = _mm_set1_epi16(1);
	const __m128i flag2_ = _mm_set1_epi16(2);
	const __m128i flag8_ = _mm_set1_epi16(0x08);
	const __m128i flag16_ = _mm_set1_epi16(0x10);
	__m128i x1_ = psw_set_low_i16(x1);
	__m128i v1_ = psw_set_low_i16(v1);
	int16_t score0[8], topen[8];

	if (zr) {
		for (t = st; t <= en; t += 8) {
			__m128i z, a, b, xt1, vt1, ut, yt, tmp, qopen_, topen_, d;
			psw_fill_pp_score_block(score0, t, r, st0, en0, qlen, tlen, m, scale_shift, qp, tf, go_ge_q, go_ge_t);
			qopen_ = psw_load_band_go_t(go_t, t, st0, en0, tlen);
			psw_fill_rev_qgap_block(topen, t, r, st0, en0, qlen, go_q);
			z = _mm_loadu_si128((const __m128i*)score0);
			topen_ = _mm_loadu_si128((const __m128i*)topen);

			xt1 = _mm_load_si128(&sv->x[t >> 3]);
			tmp = _mm_srli_si128(xt1, 14);
			xt1 = _mm_or_si128(_mm_slli_si128(xt1, 2), x1_);
			x1_ = tmp;
			vt1 = _mm_load_si128(&sv->v[t >> 3]);
			tmp = _mm_srli_si128(vt1, 14);
			vt1 = _mm_or_si128(_mm_slli_si128(vt1, 2), v1_);
			v1_ = tmp;
			a = _mm_add_epi16(xt1, vt1);

			ut = _mm_load_si128(&sv->u[t >> 3]);
			yt = _mm_load_si128(&sv->y[t >> 3]);
			b = _mm_add_epi16(yt, ut);

			d = _mm_and_si128(_mm_cmpgt_epi16(a, z), flag1_);
			z = _mm_max_epi16(z, a);
			tmp = _mm_cmpgt_epi16(b, z);
			d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, flag2_));
			z = _mm_max_epi16(z, b);

			_mm_store_si128(&sv->u[t >> 3], _mm_sub_epi16(z, vt1));
			_mm_store_si128(&sv->v[t >> 3], _mm_sub_epi16(z, ut));

			tmp = _mm_add_epi16(_mm_sub_epi16(a, z), qopen_);
			d = _mm_or_si128(d, _mm_and_si128(flag8_, _mm_cmpgt_epi16(tmp, zero_)));
			_mm_store_si128(&sv->x[t >> 3], _mm_max_epi16(tmp, zero_));
			tmp = _mm_add_epi16(_mm_sub_epi16(b, z), topen_);
			d = _mm_or_si128(d, _mm_and_si128(flag16_, _mm_cmpgt_epi16(tmp, zero_)));
			_mm_store_si128(&sv->y[t >> 3], _mm_max_epi16(tmp, zero_));
			psw_store_flags8_band(zr, t, st0, en0, d);
		}
	} else {
		for (t = st; t <= en; t += 8) {
			__m128i z, a, b, xt1, vt1, ut, yt, tmp, qopen_, topen_;
			psw_fill_pp_score_block(score0, t, r, st0, en0, qlen, tlen, m, scale_shift, qp, tf, go_ge_q, go_ge_t);
			qopen_ = psw_load_band_go_t(go_t, t, st0, en0, tlen);
			psw_fill_rev_qgap_block(topen, t, r, st0, en0, qlen, go_q);
			z = _mm_loadu_si128((const __m128i*)score0);
			topen_ = _mm_loadu_si128((const __m128i*)topen);

			xt1 = _mm_load_si128(&sv->x[t >> 3]);
			tmp = _mm_srli_si128(xt1, 14);
			xt1 = _mm_or_si128(_mm_slli_si128(xt1, 2), x1_);
			x1_ = tmp;
			vt1 = _mm_load_si128(&sv->v[t >> 3]);
			tmp = _mm_srli_si128(vt1, 14);
			vt1 = _mm_or_si128(_mm_slli_si128(vt1, 2), v1_);
			v1_ = tmp;
			a = _mm_add_epi16(xt1, vt1);

			ut = _mm_load_si128(&sv->u[t >> 3]);
			yt = _mm_load_si128(&sv->y[t >> 3]);
			b = _mm_add_epi16(yt, ut);

			z = _mm_max_epi16(z, a);
			z = _mm_max_epi16(z, b);
			_mm_store_si128(&sv->u[t >> 3], _mm_sub_epi16(z, vt1));
			_mm_store_si128(&sv->v[t >> 3], _mm_sub_epi16(z, ut));
			tmp = _mm_add_epi16(_mm_sub_epi16(a, z), qopen_);
			_mm_store_si128(&sv->x[t >> 3], _mm_max_epi16(tmp, zero_));
			tmp = _mm_add_epi16(_mm_sub_epi16(b, z), topen_);
			_mm_store_si128(&sv->y[t >> 3], _mm_max_epi16(tmp, zero_));
		}
	}
}

static inline void psw_sse_core_ps(int r, int st0, int en0, int st, int en,
                                   psw_sse_state_t *sv, uint8_t *zr,
                                   int16_t go_q, int16_t ge_q,
                                   const int16_t *go_t, const int16_t *ge_t,
                                   int qlen, int tlen, const uint8_t *query,
                                   const int16_t *tp,
                                   int16_t x1, int16_t v1)
{
	int t;
	const __m128i zero_ = _mm_setzero_si128();
	const __m128i flag1_ = _mm_set1_epi16(1);
	const __m128i flag2_ = _mm_set1_epi16(2);
	const __m128i flag8_ = _mm_set1_epi16(0x08);
	const __m128i flag16_ = _mm_set1_epi16(0x10);
	const __m128i topen_const_ = _mm_set1_epi16(go_q);
	__m128i x1_ = psw_set_low_i16(x1);
	__m128i v1_ = psw_set_low_i16(v1);
	int16_t score0[8];

	if (zr) {
		for (t = st; t <= en; t += 8) {
			__m128i z, a, b, xt1, vt1, ut, yt, tmp, qopen_, d;
			psw_fill_ps_score_block(score0, t, r, st0, en0, qlen, tlen, query, tp, go_q, ge_q, go_t, ge_t);
			qopen_ = psw_load_band_go_t(go_t, t, st0, en0, tlen);
			z = _mm_loadu_si128((const __m128i*)score0);

			xt1 = _mm_load_si128(&sv->x[t >> 3]);
			tmp = _mm_srli_si128(xt1, 14);
			xt1 = _mm_or_si128(_mm_slli_si128(xt1, 2), x1_);
			x1_ = tmp;
			vt1 = _mm_load_si128(&sv->v[t >> 3]);
			tmp = _mm_srli_si128(vt1, 14);
			vt1 = _mm_or_si128(_mm_slli_si128(vt1, 2), v1_);
			v1_ = tmp;
			a = _mm_add_epi16(xt1, vt1);

			ut = _mm_load_si128(&sv->u[t >> 3]);
			yt = _mm_load_si128(&sv->y[t >> 3]);
			b = _mm_add_epi16(yt, ut);

			d = _mm_and_si128(_mm_cmpgt_epi16(a, z), flag1_);
			z = _mm_max_epi16(z, a);
			tmp = _mm_cmpgt_epi16(b, z);
			d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, flag2_));
			z = _mm_max_epi16(z, b);

			_mm_store_si128(&sv->u[t >> 3], _mm_sub_epi16(z, vt1));
			_mm_store_si128(&sv->v[t >> 3], _mm_sub_epi16(z, ut));

			tmp = _mm_add_epi16(_mm_sub_epi16(a, z), qopen_);
			d = _mm_or_si128(d, _mm_and_si128(flag8_, _mm_cmpgt_epi16(tmp, zero_)));
			_mm_store_si128(&sv->x[t >> 3], _mm_max_epi16(tmp, zero_));
			tmp = _mm_add_epi16(_mm_sub_epi16(b, z), topen_const_);
			d = _mm_or_si128(d, _mm_and_si128(flag16_, _mm_cmpgt_epi16(tmp, zero_)));
			_mm_store_si128(&sv->y[t >> 3], _mm_max_epi16(tmp, zero_));
			psw_store_flags8_band(zr, t, st0, en0, d);
		}
	} else {
		for (t = st; t <= en; t += 8) {
			__m128i z, a, b, xt1, vt1, ut, yt, tmp, qopen_, topen_;
			psw_fill_ps_score_block(score0, t, r, st0, en0, qlen, tlen, query, tp, go_q, ge_q, go_t, ge_t);
			qopen_ = psw_load_band_go_t(go_t, t, st0, en0, tlen);
			z = _mm_loadu_si128((const __m128i*)score0);

			xt1 = _mm_load_si128(&sv->x[t >> 3]);
			tmp = _mm_srli_si128(xt1, 14);
			xt1 = _mm_or_si128(_mm_slli_si128(xt1, 2), x1_);
			x1_ = tmp;
			vt1 = _mm_load_si128(&sv->v[t >> 3]);
			tmp = _mm_srli_si128(vt1, 14);
			vt1 = _mm_or_si128(_mm_slli_si128(vt1, 2), v1_);
			v1_ = tmp;
			a = _mm_add_epi16(xt1, vt1);

			ut = _mm_load_si128(&sv->u[t >> 3]);
			yt = _mm_load_si128(&sv->y[t >> 3]);
			b = _mm_add_epi16(yt, ut);

			z = _mm_max_epi16(z, a);
			z = _mm_max_epi16(z, b);
			_mm_store_si128(&sv->u[t >> 3], _mm_sub_epi16(z, vt1));
			_mm_store_si128(&sv->v[t >> 3], _mm_sub_epi16(z, ut));
			tmp = _mm_add_epi16(_mm_sub_epi16(a, z), qopen_);
			_mm_store_si128(&sv->x[t >> 3], _mm_max_epi16(tmp, zero_));
			tmp = _mm_add_epi16(_mm_sub_epi16(b, z), topen_const_);
			_mm_store_si128(&sv->y[t >> 3], _mm_max_epi16(tmp, zero_));
		}
	}
}

float psw_gg3_sse_pp(void *km, int qlen, const psw_prof_t *query,
                     int tlen, const psw_prof_t *target,
                     int8_t m, const int8_t *mat,
                     int8_t gapo, int8_t gape, int w,
                     int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	int16_t *qp = 0, *tf = 0, *qbf = 0, *tbf = 0;
	int16_t *go_q = 0, *ge_q = 0, *go_t = 0, *ge_t = 0, *go_ge_q = 0, *go_ge_t = 0;
	psw_prof_t query_n = {0, 0, 0, 0}, target_n = {0, 0, 0, 0};
	const psw_prof_t *q_use, *t_use;
	int32_t r, n_col, *off = 0;
	int32_t score32 = INT32_MIN;
	int32_t H0 = 0, last_H0_t = 0;
	uint8_t *z = 0;
	int16_t scale;
	int8_t scale_shift;
	psw_sse_state_t sv = {0};
	int ok = 0;
	int failed = 0;

	if (gapo < 0 || gape < 0) return PSW_NEG_INF_F;
	if (query == 0 || target == 0 || mat == 0) return PSW_NEG_INF_F;
	if (query->prof == 0 || target->prof == 0) return PSW_NEG_INF_F;
	if (qlen < 0 || tlen < 0 || m <= 0) return PSW_NEG_INF_F;
	if (query->len < qlen || target->len < tlen) return PSW_NEG_INF_F;
	if (query->dim < m || target->dim < m) return PSW_NEG_INF_F;
	scale = psw_pick_safe_scale_pow2(m, mat, gapo, gape);
	if (scale < 2) return PSW_NEG_INF_F;
	scale_shift = psw_scale_to_shift(scale);
	if (scale_shift <= 0) return PSW_NEG_INF_F;
	if (!psw_make_norm_prof(km, query, &query_n, m, scale)) return PSW_NEG_INF_F;
	if (!psw_make_norm_prof(km, target, &target_n, m, scale)) {
		psw_free_norm_prof(km, &query_n);
		return PSW_NEG_INF_F;
	}
	q_use = &query_n;
	t_use = &target_n;

	if (w < 0) w = tlen > qlen ? tlen : qlen;
	n_col = w + 1 < tlen ? w + 1 : tlen;

	do {
		qp = psw_gen_qp_i16(km, qlen, q_use, m, mat);
		if (qp == 0) { failed = 1; break; }
		tf = psw_gen_tf_i16(km, tlen, t_use, m);
		if (tf == 0) { failed = 1; break; }
		qbf = psw_gen_base_freq_i16(km, qlen, q_use, m, scale);
		if (qbf == 0) { failed = 1; break; }
		tbf = psw_gen_base_freq_i16(km, tlen, t_use, m, scale);
		if (tbf == 0) { failed = 1; break; }

		go_q = (int16_t*)kmalloc(km, (size_t)qlen * sizeof(int16_t));
		ge_q = (int16_t*)kmalloc(km, (size_t)qlen * sizeof(int16_t));
		go_t = (int16_t*)kmalloc(km, (size_t)tlen * sizeof(int16_t));
		ge_t = (int16_t*)kmalloc(km, (size_t)tlen * sizeof(int16_t));
		go_ge_q = (int16_t*)kmalloc(km, (size_t)qlen * sizeof(int16_t));
		go_ge_t = (int16_t*)kmalloc(km, (size_t)tlen * sizeof(int16_t));
		if (go_q == 0 || ge_q == 0 || go_t == 0 || ge_t == 0 || go_ge_q == 0 || go_ge_t == 0) { failed = 1; break; }
		for (r = 0; r < qlen; ++r) {
			go_q[r] = psw_sat16((int32_t)gapo * qbf[r]);
			ge_q[r] = psw_sat16((int32_t)gape * qbf[r]);
			go_ge_q[r] = (int16_t)(go_q[r] + ge_q[r]);
		}
		for (r = 0; r < tlen; ++r) {
			go_t[r] = psw_sat16((int32_t)gapo * tbf[r]);
			ge_t[r] = psw_sat16((int32_t)gape * tbf[r]);
			go_ge_t[r] = (int16_t)(go_t[r] + ge_t[r]);
		}
		if (!psw_sse_state_init(km, tlen, &sv)) { failed = 1; break; }
		ok = 1;

		if (m_cigar_ && n_cigar_ && cigar_) {
			*n_cigar_ = 0;
			z = (uint8_t*)kcalloc(km, (size_t)(qlen + tlen) * n_col, 1);
			off = (int32_t*)kmalloc(km, (size_t)(qlen + tlen) * sizeof(int32_t));
			if (z == 0 || off == 0) { failed = 1; break; }
		}

		if (qlen == 0 || tlen == 0) {
			int32_t s32 = qlen == 0 ? psw_gap_only_target_i16(go_t, ge_t, tlen)
			                       : psw_gap_only_query_i16(go_q, ge_q, qlen);
			score32 = s32;
			if (m_cigar_ && n_cigar_ && cigar_) {
				*n_cigar_ = 0;
				if (qlen == 0 && tlen > 0)
					*cigar_ = psw_push_cigar(km, n_cigar_, m_cigar_, *cigar_, PSW_CIGAR_DEL, tlen);
				else if (tlen == 0 && qlen > 0)
					*cigar_ = psw_push_cigar(km, n_cigar_, m_cigar_, *cigar_, PSW_CIGAR_INS, qlen);
			}
			break;
		}

		for (r = 0; r < qlen + tlen - 1; ++r) {
			int32_t st = 0, en = tlen - 1, st0, en0;
			int16_t x1, v1;
			if (st < r - qlen + 1) st = r - qlen + 1;
			if (en > r) en = r;
			if (st < (r - w + 1) >> 1) st = (r - w + 1) >> 1;
			if (en > (r + w) >> 1) en = (r + w) >> 1;
			if (st > en) continue;
			st0 = st; en0 = en;
			st = (st / 8) * 8;
			en = ((en + 8) / 8) * 8 - 1;

			if (st0 != 0) {
				if (r > st0 + st0 + w - 1) x1 = v1 = 0;
				else {
					x1 = sv.x16[st0 - 1];
					v1 = sv.v16[st0 - 1];
				}
			} else {
				x1 = 0;
				v1 = r ? go_q[r - 1] : 0;
			}
			if (en0 != r) {
				if (r < en0 + en0 - w - 1) sv.y16[en0] = sv.u16[en0] = 0;
			} else {
				sv.y16[r] = 0;
				sv.u16[r] = r ? go_t[r - 1] : 0;
			}
			if (z) off[r] = st0;
			psw_sse_core_pp((int)r, (int)st0, (int)en0, (int)st, (int)en, &sv, z ? z + (size_t)r * n_col : 0,
			                go_t, go_q, qlen, tlen, m, scale_shift, qp, tf, go_ge_q, go_ge_t, x1, v1);

			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0) {
					int32_t jh = r - last_H0_t;
					H0 += (int32_t)sv.v16[last_H0_t] - go_ge_q[jh];
				} else {
					++last_H0_t;
					H0 += (int32_t)sv.u16[last_H0_t] - go_ge_t[last_H0_t];
				}
			} else {
				H0 = (int32_t)sv.v16[0] - ((int32_t)go_ge_t[0] + go_ge_q[0]);
				last_H0_t = 0;
			}
		}
		score32 = H0;
		if (z && off)
			psw_backtrack(km, 1, 0, 0, z, off, 0, n_col, tlen - 1, qlen - 1, m_cigar_, n_cigar_, cigar_);
	} while (0);

	if (failed) score32 = INT32_MIN;
	if (z) kfree(km, z);
	if (off) kfree(km, off);
	if (ok) psw_sse_state_destroy(km, &sv);
	if (go_q) kfree(km, go_q);
	if (ge_q) kfree(km, ge_q);
	if (go_t) kfree(km, go_t);
	if (ge_t) kfree(km, ge_t);
	if (go_ge_q) kfree(km, go_ge_q);
	if (go_ge_t) kfree(km, go_ge_t);
	if (qp) kfree(km, qp);
	if (tf) kfree(km, tf);
	if (qbf) kfree(km, qbf);
	if (tbf) kfree(km, tbf);
	psw_free_norm_prof(km, &query_n);
	psw_free_norm_prof(km, &target_n);
	if (score32 == INT32_MIN) return PSW_NEG_INF_F;
	return (float)score32 / (float)scale;
}

float psw_gg3_sse_ps(void *km, int qlen, const uint8_t *query,
                     int tlen, const psw_prof_t *target,
                     int8_t m, const int8_t *mat,
                     int8_t gapo, int8_t gape, int w,
                     int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	int16_t *tp = 0, *tbf = 0, *go_t = 0, *ge_t = 0;
	psw_prof_t target_n = {0, 0, 0, 0};
	const psw_prof_t *t_use;
	int32_t r, n_col, *off = 0;
	int32_t score32 = INT32_MIN;
	int32_t H0 = 0, last_H0_t = 0;
	uint8_t *z = 0;
	int16_t scale, go_q, ge_q;
	int8_t scale_shift;
	psw_sse_state_t sv = {0};
	int ok = 0;
	int failed = 0;

	if (gapo < 0 || gape < 0) return PSW_NEG_INF_F;
	if (query == 0 || target == 0 || mat == 0) return PSW_NEG_INF_F;
	if (target->prof == 0) return PSW_NEG_INF_F;
	if (qlen < 0 || tlen < 0 || m <= 0) return PSW_NEG_INF_F;
	if (target->len < tlen || target->dim < m) return PSW_NEG_INF_F;
	for (r = 0; r < qlen; ++r)
		if ((int)query[r] >= m) return PSW_NEG_INF_F;

	scale = psw_pick_safe_scale_pow2(m, mat, gapo, gape);
	if (scale < 2) return PSW_NEG_INF_F;
	scale_shift = psw_scale_to_shift(scale);
	if (scale_shift <= 0) return PSW_NEG_INF_F;
	go_q = psw_sat16((int32_t)gapo * scale);
	ge_q = psw_sat16((int32_t)gape * scale);
	if (!psw_make_norm_prof(km, target, &target_n, m, scale)) return PSW_NEG_INF_F;
	t_use = &target_n;

	if (w < 0) w = tlen > qlen ? tlen : qlen;
	n_col = w + 1 < tlen ? w + 1 : tlen;

	do {
		tp = psw_gen_tp_i16(km, tlen, t_use, m, mat);
		if (tp == 0) { failed = 1; break; }
		tbf = psw_gen_base_freq_i16(km, tlen, t_use, m, scale);
		if (tbf == 0) { failed = 1; break; }
		go_t = (int16_t*)kmalloc(km, (size_t)tlen * sizeof(int16_t));
		ge_t = (int16_t*)kmalloc(km, (size_t)tlen * sizeof(int16_t));
		if (go_t == 0 || ge_t == 0) { failed = 1; break; }
		for (r = 0; r < tlen; ++r) {
			go_t[r] = psw_sat16((int32_t)gapo * tbf[r]);
			ge_t[r] = psw_sat16((int32_t)gape * tbf[r]);
		}
		if (!psw_sse_state_init(km, tlen, &sv)) { failed = 1; break; }
		ok = 1;

		if (m_cigar_ && n_cigar_ && cigar_) {
			*n_cigar_ = 0;
			z = (uint8_t*)kcalloc(km, (size_t)(qlen + tlen) * n_col, 1);
			off = (int32_t*)kmalloc(km, (size_t)(qlen + tlen) * sizeof(int32_t));
			if (z == 0 || off == 0) { failed = 1; break; }
		}

		if (qlen == 0 || tlen == 0) {
			int32_t s32 = qlen == 0 ? psw_gap_only_target_i16(go_t, ge_t, tlen)
			                       : psw_gap_only_query_scalar_i16(go_q, ge_q, qlen);
			score32 = s32;
			if (m_cigar_ && n_cigar_ && cigar_) {
				*n_cigar_ = 0;
				if (qlen == 0 && tlen > 0)
					*cigar_ = psw_push_cigar(km, n_cigar_, m_cigar_, *cigar_, PSW_CIGAR_DEL, tlen);
				else if (tlen == 0 && qlen > 0)
					*cigar_ = psw_push_cigar(km, n_cigar_, m_cigar_, *cigar_, PSW_CIGAR_INS, qlen);
			}
			break;
		}

		for (r = 0; r < qlen + tlen - 1; ++r) {
			int32_t st = 0, en = tlen - 1, st0, en0;
			int16_t x1, v1;
			if (st < r - qlen + 1) st = r - qlen + 1;
			if (en > r) en = r;
			if (st < (r - w + 1) >> 1) st = (r - w + 1) >> 1;
			if (en > (r + w) >> 1) en = (r + w) >> 1;
			if (st > en) continue;
			st0 = st; en0 = en;
			st = (st / 8) * 8;
			en = ((en + 8) / 8) * 8 - 1;

			if (st0 != 0) {
				if (r > st0 + st0 + w - 1) x1 = v1 = 0;
				else {
					x1 = sv.x16[st0 - 1];
					v1 = sv.v16[st0 - 1];
				}
			} else {
				x1 = 0;
				v1 = r ? go_q : 0;
			}
			if (en0 != r) {
				if (r < en0 + en0 - w - 1) sv.y16[en0] = sv.u16[en0] = 0;
			} else {
				sv.y16[r] = 0;
				sv.u16[r] = r ? go_t[r - 1] : 0;
			}
			if (z) off[r] = st0;
			psw_sse_core_ps((int)r, (int)st0, (int)en0, (int)st, (int)en, &sv, z ? z + (size_t)r * n_col : 0,
			                go_q, ge_q, go_t, ge_t, qlen, tlen, query, tp, x1, v1);

			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0)
					H0 += (int32_t)sv.v16[last_H0_t] - ((int32_t)go_q + ge_q);
				else {
					++last_H0_t;
					H0 += (int32_t)sv.u16[last_H0_t] - ((int32_t)go_t[last_H0_t] + ge_t[last_H0_t]);
				}
			} else {
				H0 = (int32_t)sv.v16[0] - ((int32_t)go_t[0] + ge_t[0] + go_q + ge_q);
				last_H0_t = 0;
			}
		}
		score32 = H0;
		if (z && off)
			psw_backtrack(km, 1, 0, 0, z, off, 0, n_col, tlen - 1, qlen - 1, m_cigar_, n_cigar_, cigar_);
	} while (0);

	if (failed) score32 = INT32_MIN;
	if (z) kfree(km, z);
	if (off) kfree(km, off);
	if (ok) psw_sse_state_destroy(km, &sv);
	if (go_t) kfree(km, go_t);
	if (ge_t) kfree(km, ge_t);
	if (tp) kfree(km, tp);
	if (tbf) kfree(km, tbf);
	psw_free_norm_prof(km, &target_n);
	if (score32 == INT32_MIN) return PSW_NEG_INF_F;
	return (float)score32 / (float)scale;
}

#else

float psw_gg3_sse_pp(void *km, int qlen, const psw_prof_t *query,
                     int tlen, const psw_prof_t *target,
                     int8_t m, const int8_t *mat,
                     int8_t gapo, int8_t gape, int w,
                     int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	return psw_gg3_pp(km, qlen, query, tlen, target, m, mat, gapo, gape, w, m_cigar_, n_cigar_, cigar_);
}

float psw_gg3_sse_ps(void *km, int qlen, const uint8_t *query,
                     int tlen, const psw_prof_t *target,
                     int8_t m, const int8_t *mat,
                     int8_t gapo, int8_t gape, int w,
                     int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	return psw_gg3_ps(km, qlen, query, tlen, target, m, mat, gapo, gape, w, m_cigar_, n_cigar_, cigar_);
}

#endif
