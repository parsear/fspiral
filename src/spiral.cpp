#include <stddef.h>
#include <vector>
#include <assert.h>
#include <math.h>

#include "encrypt.h"
#include "answer.h"
#include "utils.h"
#include "spiral.h"
#include "samples.h"

#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

// include __m512i
#if defined(__AVX512F__)
#include <immintrin.h>
#include <avx512fintrin.h>
#endif

#if defined(__AVX2__)
#define USE_AVX2
#include <immintrin.h>
#endif

// #define SPIRAL_TIME_COLLECT

#ifdef SPIRAL_TIME_COLLECT
#include <iostream>
#include <chrono>
#endif

/**
 * RGSWCiphertext class
*/

RNSRGSWCiphertext::RNSRGSWCiphertext(/* args */)
{
    this->a.resize(2 * ellnum);
    this->b.resize(2 * ellnum);
    this->arns.resize(2 * ellnum);
    this->brns.resize(2 * ellnum);
    uint64_t length = this->getLength();
    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        this->a[i].resize(length);
        this->b[i].resize(length);
        this->arns[i].resize(length);
        this->brns[i].resize(length);
    }
#ifdef INTEL_HEXL

    if (modulus == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, modulus);
        this->ntts = ntts;
    }

    intel::hexl::NTT ntts2(length, modulus2);
    this->ntts2 = ntts2;
#endif

}

RNSRGSWCiphertext::RNSRGSWCiphertext(int32_t len, uint64_t module) : length(len), modulus(module) 
{
    this->a.resize(2 * ellnum);
    this->b.resize(2 * ellnum);
    this->arns.resize(2 * ellnum);
    this->brns.resize(2 * ellnum);
    uint64_t length = this->getLength();
    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        this->a[i].resize(length);
        this->b[i].resize(length);
        this->arns[i].resize(length);
        this->brns[i].resize(length);
    }
    // sample_random(this->a, module);
    // TODO: should be removed
    // sample_random(this->b, module);


#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }

    intel::hexl::NTT ntts2(length, modulus2);
    this->ntts2 = ntts2;
#endif
}

RNSRGSWCiphertext::RNSRGSWCiphertext(int32_t len, uint64_t module, int32_t elln, uint64_t base, uint64_t bg):
    length(len), modulus(module), ellnum(elln), PP(base), BBg(bg)
{
    this->a.resize(2 * ellnum);
    this->b.resize(2 * ellnum);
    this->arns.resize(2 * ellnum);
    this->brns.resize(2 * ellnum);
    uint64_t length = this->getLength();
    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        this->a[i].resize(length);
        this->b[i].resize(length);
        this->arns[i].resize(length);
        this->brns[i].resize(length);
    }

#ifdef INTEL_HEXL
    if (module == crtMod)
    {
        intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
        this->ntts = ntts;
    } else {
        intel::hexl::NTT ntts(length, module);
        this->ntts = ntts;
    }

    intel::hexl::NTT ntts2(length, modulus2);
    this->ntts2 = ntts2;
#endif
} 

void RNSRGSWCiphertext::keyGen(Secret& secret)
{    
    int32_t length = secret.getLength();
    uint64_t modulus = secret.getModulus();
    uint64_t modulus2 = this->getSecModulus();

    std::vector<uint64_t> message(length);
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }
    // encrypt -s to convertKey
    for (size_t i = 0; i < length; i++)
    {
        // the value `modulus' can be reduced in the following steps?
        message[i] = modulus - secret.getData(i);
    }

    std::vector<uint64_t> message2(length);
    copy(message.begin(), message.end(), message2.begin());
    guass_to_modulus(message2.data(), modulus, modulus2);

    // uint64_t mult_factor;
    std::vector<uint64_t> temp1(length), temp2(length);

    // m(s, 1)
    intel::hexl::NTT ntts = this->getNTT();
    intel::hexl::NTT ntts2 = this->getNTT2();

    std::vector<uint64_t> secret_ntt1(length, 0);
    std::vector<uint64_t> secret_ntt2(length, 0);
    for (size_t i = 0; i < length; i++) secret_ntt1[i] = secret.getData(i);
    copy(secret_ntt1.begin(), secret_ntt1.end(), secret_ntt2.begin());
    guass_to_modulus(secret_ntt2.data(), modulus, modulus2);
    ntts.ComputeForward(secret_ntt1.data(), secret_ntt1.data(), 1, 1);
    ntts2.ComputeForward(secret_ntt2.data(), secret_ntt2.data(), 1, 1);

    std::vector<uint128_t> mult_factor = powerOfBg128(this->PP, this->BBg, this->ellnum);
    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        // faster by single barret reduce
        intel::hexl::EltwiseFMAMod(temp1.data(), message.data(),
                mult_factor[i] % this->modulus, nullptr, this->length, this->modulus, 1);
        intel::hexl::EltwiseFMAMod(temp2.data(), message2.data(),
                mult_factor[i] % this->modulus2, nullptr, this->length, this->modulus2, 1);

        // rlwe encryption, these rlwes is in ntt form.
        // encrypt(m)
        encrypt_rns(this->b[ellnum + i], this->a[ellnum + i],
                    this->brns[ellnum + i], this->arns[ellnum + i],
                    secret, temp1, temp2, 1, this->modulus2);

        // encrypt(-ms)
        negate(temp1.data(), this->length, this->modulus);
        negate(temp2.data(), this->length, this->modulus2);

        ntts.ComputeForward(temp1.data(), temp1.data(), 1, 1);
        intel::hexl::EltwiseMultMod(temp1.data(), temp1.data(), secret_ntt1.data(),
                this->length, this->modulus, 1);
        ntts.ComputeInverse(temp1.data(), temp1.data(), 1, 1);

        ntts2.ComputeForward(temp2.data(), temp2.data(), 1, 1);
        intel::hexl::EltwiseMultMod(temp2.data(), temp2.data(), secret_ntt2.data(),
                this->length, this->modulus2, 1);
        ntts2.ComputeInverse(temp2.data(), temp2.data(), 1, 1);

        // encrypt(this->b[i].data(), this->a[i].data(), secret, temp.data());
        encrypt_rns(this->b[i], this->a[i],
                    this->brns[i], this->arns[i],
                    secret, temp1, temp2, 1, this->modulus2);
#endif
        this->isntt = true;
    }
}

uint64_t RNSRGSWCiphertext::getModulus() const
{
    return this->modulus;
}

uint64_t RNSRGSWCiphertext::getSecModulus() const
{
    return this->modulus2;
}

int32_t RNSRGSWCiphertext::getLength() const
{
    return this->length;
}

int32_t RNSRGSWCiphertext::getEllnum() const
{
    return this->ellnum;
}

uint64_t RNSRGSWCiphertext::getBg() const
{
    return this->BBg;
}

uint64_t RNSRGSWCiphertext::getBase() const
{
    return this->PP;
}

bool RNSRGSWCiphertext::getIsNtt() const
{
    return this->isntt;
}

void RNSRGSWCiphertext::setIsNtt(bool is)
{
    this->isntt = is;
}

#ifdef INTEL_HEXL
intel::hexl::NTT RNSRGSWCiphertext::getNTT() const
{
    return this->ntts;
}

intel::hexl::NTT RNSRGSWCiphertext::getNTT2() const
{
    return this->ntts2;
}
#endif



// for Spiral
void get_spiral_parameters_one(spiral_parameters& para, int32_t rlwesNum, int32_t rgswsNum, uint64_t spiral_p)
{
    para.rlwesNum = rlwesNum;
    para.rgswsNum = rgswsNum;
    para.spiral_p = spiral_p;

    para.ellgsw = 4;
    para.bggsw = 0x01 << 14;
}

// for SpiralStream
void get_spiralstream_parameters_one(spiral_parameters& para, int32_t rlwesNum, int32_t rgswsNum, uint64_t spiral_p)
{
    para.rlwesNum = rlwesNum;
    para.rgswsNum = rgswsNum;
    para.spiral_p = spiral_p;

    para.ellgsw = 3;
    para.bggsw = 0x01 << 14;
    para.basegsw = 0x01 << 14;
}

void get_spiral_rns_parameters_one(spiral_parameters& para, int32_t rlwesNum, int32_t rgswsNum, uint64_t spiral_p, int32_t n)
{
    para.rlwesNum = rlwesNum;
    para.rgswsNum = rgswsNum;
    para.spiral_p = spiral_p;

    para.ell_rlwes_expansion = 6; // 8;
    para.bg_rlwes_expansion = 0x01 << 15; // 0x01 << 10;
    para.ell_rgsws_expansion = 6; // 8;
    para.bg_rgsws_expansion = 0x01 << 15; // 0x01 << 10;

    para.ell_convert = 7;
    para.bg_convert = 0x01 << 10;
    para.base_convert = 0x01 << 10;

    para.ellgsw = 4;
    para.bggsw = 0x01 << 14;

    para.n = n;
}

/****  the following algorithms are from Spiral   ****/
/**
constexpr __uint128_t b_inv_pa_i = 163640210UL * crtq2;
constexpr __uint128_t pa_inv_b_i = 97389680UL * crtq1;
constexpr uint64_t cr0_Q = 7906011006380390721UL;
constexpr uint64_t cr1_Q = 275UL; 

typedef struct __attribute__((packed, aligned(16))) {
    uint64_t x, y;
} ulonglong2_h;

inline ulonglong2_h umul64wide(uint64_t a, uint64_t b) {
    ulonglong2_h z;
    __uint128_t val = ((__uint128_t)a) * ((__uint128_t)b);
    z.x = (uint64_t)val;
    z.y = (uint64_t)(val >> 64);
    return z;
}

inline uint64_t cpu_add_u64(uint64_t operand1, uint64_t operand2,
                                uint64_t *result) {
    *result = operand1 + operand2;
    return (*result < operand1) ? 0x01 : 0x00; // detect overflow
}

inline uint64_t barrett_raw_u128(__uint128_t val, uint64_t const_ratio_0, uint64_t const_ratio_1, uint64_t modulus) {
    uint64_t zx = val & (((__uint128_t)1 << 64) - 1);
    uint64_t zy = val >> 64;

    uint64_t tmp1, tmp3, carry;
    ulonglong2_h prod = umul64wide(zx, const_ratio_0);
    carry = prod.y;
    ulonglong2_h tmp2 = umul64wide(zx, const_ratio_1);
    tmp3 = tmp2.y + cpu_add_u64(tmp2.x, carry, &tmp1);
    tmp2 = umul64wide(zy, const_ratio_0);
    carry = tmp2.y + cpu_add_u64(tmp1, tmp2.x, &tmp1);
    tmp1 = zy * const_ratio_1 + tmp3 + carry;
    tmp3 = zx - tmp1 * modulus;

    return tmp3;
}

inline uint64_t barrett_reduction_u128(__uint128_t val) {
    uint64_t reduced_val = barrett_raw_u128(val, cr0_Q, cr1_Q, crtMod);
    reduced_val -= (crtMod)*(static_cast<int>(reduced_val >= crtMod));
    return reduced_val;
}

uint64_t crt_compose(uint64_t x, uint64_t y) {
    // uint64_t val_ap = x * a_inv_p_i;
    // val_ap += y * p_inv_a_i;
    uint64_t val_ap_u64 = x;
    __uint128_t val = val_ap_u64 * b_inv_pa_i;
    val += y * pa_inv_b_i;
    // cout << "(" << x << ", " << y << ", " << z << ") -> " << (uint64_t)(val % Q_i) << endl;
    // return val % Q_i;
    return barrett_reduction_u128(val);
}
**/
inline uint64_t barrett_raw_u64(uint64_t input, uint64_t const_ratio_1, uint64_t modulus)
{
    unsigned long long tmp[2];
    tmp[1] = static_cast<unsigned long long>(((static_cast<__uint128_t>(input) * static_cast<__uint128_t>(const_ratio_1)) >> 64));

    // Barrett subtraction
    tmp[0] = input - tmp[1] * modulus;

    // One more subtraction is enough
    return (tmp[0] >= modulus ? (tmp[0] - modulus) : (tmp[0]));
}

inline uint64_t barrett_coeff(uint64_t val, size_t n)
{
    // return val % moduli[n];
    return (n == 0) ? barrett_raw_u64(val, cr1_p, crtq1) : barrett_raw_u64(val, cr1_b, crtq2);
    return val;
}

void fastMultiplyQueryByDatabaseDim1InvCRTCipher(
    // std::vector<std::vector<uint64_t> >& out,
    std::vector<RlweCiphertext>& out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    size_t dim0, size_t num_per)
{
    // v_firstdim: 
    //     poly_len, dim0, ct.rows, ct.cols
    // db: 
    //     poly_len, num_per, pt.cols, dim0, pt.rows
    size_t base_dim = 2;

    size_t ct_rows = base_dim;
    size_t ct_cols = 1;
    size_t pt_rows = 1;
    size_t pt_cols = 1;

    #if defined(__AVX512F__) && !defined(NO_CRT)
        assert(dim0 * ct_rows >= max_summed_pa_or_b_in_u64);

        for (size_t z = 0; z < poly_len; z++) {
            // auto start  = std::chrono::high_resolution_clock::now();
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            // if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                // for (size_t c = 0; c < pt_cols; c++) {
                    size_t uroll_max = max_summed_pa_or_b_in_u64;
                    size_t inner_limit = uroll_max;
                    size_t outer_limit = dim0 * ct_rows / inner_limit;

                    uint64_t sums_out_n0_u64_acc[4] = { 0, 0, 0, 0 };
                    uint64_t sums_out_n2_u64_acc[4] = { 0, 0, 0, 0 };

                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m512i sums_out_n0 = _mm512_setzero_si512();
                        __m512i sums_out_n2 = _mm512_setzero_si512();

                        // #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/8; i_jm++) {
                            size_t jm = o_jm * inner_limit + (8*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            uint64_t b_inp_3 = db[idx_b_base++];
                            uint64_t b_inp_4 = db[idx_b_base++];
                            __m512i b  = _mm512_set_epi64(b_inp_4, b_inp_4, b_inp_3, b_inp_3, b_inp_2, b_inp_2, b_inp_1, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m512i a = _mm512_load_si512(v_a);
                            __m512i a_lo = a;
                            __m512i a_hi_hi = _mm512_srli_epi64(a, packed_offset_2);
                            __m512i b_lo = b;
                            __m512i b_hi_hi = _mm512_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_mul_epu32(a_hi_hi, b_hi_hi));
                        }
                        
                        // reduce here, otherwise we will overflow
                        alignas(64) uint64_t sums_out_n0_u64[8];
                        alignas(64) uint64_t sums_out_n2_u64[8];
                        _mm512_store_si512(sums_out_n0_u64, sums_out_n0);
                        _mm512_store_si512(sums_out_n2_u64, sums_out_n2);

                        // TODO: check
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n0_u64[idx] + sums_out_n0_u64[4 + idx];
                            // sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val) % crtq1;
                            sums_out_n0_u64_acc[idx] += val;
                        }
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n2_u64[idx] + sums_out_n2_u64[4 + idx];
                            // sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val) % crtq2;
                            sums_out_n2_u64_acc[idx] += val;
                        }
                    }
                    // TODO: check
                    /**
                    for (size_t idx = 0; idx < 4; idx++) {
                        sums_out_n0_u64_acc[idx] %= crtq1;
                        sums_out_n2_u64_acc[idx] %= crtq2;
                    }
                    **/

                    /**
                    // output n0
                    size_t n = 0;
                    size_t idx_c =  n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);

                    idx_c += (crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);


                    // output n1
                    n = 1;
                    idx_c = n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);

                    idx_c += (crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);
                    **/

                    /**
                    uint64_t temp1 = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);
                    uint64_t temp2 = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);
                    uint64_t temp3 = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);
                    uint64_t temp4 = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);

                    out[i][z] = crt_compose(temp1, temp3, 0);
                    out[i][poly_len + z] = crt_compose(temp2, temp4, 0);;
                    **/
                   /**
                    out[i][z] = crt_compose(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2],
                            sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2]);
                    out[i][poly_len + z] = crt_compose(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3],
                            sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3]);
                    **/
                    out[i].a[z] = crt_compose(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2],
                            sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2]);
                    out[i].b[z] = crt_compose(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3],
                            sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3]);
                // }
            }
            // auto stop  = std::chrono::high_resolution_clock::now();

            // auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            // std::cout << " (in total) database multiplication total costs " << glapsed.count() << " us." << std::endl;
        }
        #elif defined(__AVX2__) && !defined(NO_CRT)
        assert(dim0 * ct_rows >= max_summed_pa_or_b_in_u64);

        for (size_t z = 0; z < poly_len; z++) {
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            // if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                for (size_t c = 0; c < pt_cols; c++) {
                    size_t uroll_max = max_summed_pa_or_b_in_u64;
                    size_t inner_limit = uroll_max;
                    size_t outer_limit = dim0 * ct_rows / inner_limit;

                    uint64_t sums_out_n0_u64_acc[4] = { 0, 0, 0, 0 };
                    uint64_t sums_out_n2_u64_acc[4] = { 0, 0, 0, 0 };

                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m256i sums_out_n0 = _mm256_setzero_si256();
                        __m256i sums_out_n2 = _mm256_setzero_si256();

                        #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/4; i_jm++) {
                            size_t jm = o_jm * inner_limit + (4*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            __m256i b  = _mm256_set_epi64x(b_inp_2, b_inp_2, b_inp_1, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m256i a = _mm256_load_si256((const __m256i *)v_a);
                            __m256i a_lo = a;
                            __m256i a_hi_hi = _mm256_srli_epi64(a, packed_offset_2);
                            __m256i b_lo = b;
                            __m256i b_hi_hi = _mm256_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm256_add_epi64(sums_out_n0, _mm256_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm256_add_epi64(sums_out_n2, _mm256_mul_epu32(a_hi_hi, b_hi_hi));
                        }
                        
                        // reduce here, otherwise we will overflow
                        alignas(64) uint64_t sums_out_n0_u64[4];
                        alignas(64) uint64_t sums_out_n2_u64[4];
                        _mm256_store_si256((__m256i *)sums_out_n0_u64, sums_out_n0);
                        _mm256_store_si256((__m256i *)sums_out_n2_u64, sums_out_n2);

                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n0_u64[idx];
                            // sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val) % crtq1;
                            sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val);
                        }
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n2_u64[idx];
                            // sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val) % crtq2;
                            sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val);
                        }
                    }
                /**
                    for (size_t idx = 0; idx < 4; idx++) {
                        sums_out_n0_u64_acc[idx] %= crtq1;
                        sums_out_n2_u64_acc[idx] %= crtq2;
                    }
                **/
                /**
                    // output n0
                    size_t n = 0;
                    size_t idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);

                    // output n1
                    n = 1;
                    idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i][idx_c] = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);
                **/
                /**
                    out[i][z] = crt_compose(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2],
                            sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2]);
                    out[i][poly_len + z] = crt_compose(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3],
                            sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3]);
                **/
                    out[i].a[z] = crt_compose(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2],
                            sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2]);
                    out[i].b[z] = crt_compose(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3],
                            sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3]);
                }
            }
        }
    #endif
}

void sample_spiral_database(std::vector<std::vector<uint64_t> >& database, uint64_t spiral_p)
{
    size_t total_size = database.size();
    for (size_t i = 0; i < total_size; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            int32_t rd = rand() % spiral_p;
            database[i][j] = rd > spiral_p/2 ? crtMod + rd - spiral_p : rd;
        }
    }

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
    for (size_t i = 0; i < total_size; i++)
    {
        ntts.ComputeForward(database[i].data(), database[i].data(), 1, 1);
    }
#endif
}

// encode to crt form
void spiral_sdatabase_tocrt(uint64_t* datacrt, std::vector<std::vector<uint64_t> >& data_ntt, int32_t N1, int32_t n)
{
    int32_t N2 = data_ntt.size() / n / n / N1;
    
    /** datacrt = [a, e, i, m, b, f, j, \cdots] 
    data_ntt = [a, b, c, d
                e, f, g, h
                i, j, k, l
                m, n, o, p]
    **/
   int64_t part_size = N * N2 * N1;
for (size_t part = 0; part < pow(n, 2); part++)
{
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N2; j++)
        {
            for (size_t k = 0; k < N1; k++)
            {
                uint64_t temp = data_ntt[j * N1 + k][i];
                // TODO: check
                datacrt[part_size * part + i * N1 * N2 + j * N1 + k] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
                // uint64_t temp1 = barrett_coeff(temp, 0);
                // uint64_t temp2 = barrett_coeff(temp, 1);
                // datacrt[i * N1 * N2 + j * N1 + k] = temp1 | temp2 << 32;
            }
        }
    }
}
}

/**
 * @brief the 
 * 
 * @param result 
 * @param constNum 
 */
void multConst(RlweCiphertext& result, uint64_t constNum)
{
    uint64_t modulus = result.getModulus();
#ifdef INTEL_HEXL
    // multiply constNum whatever it is ntt form or not
    intel::hexl::EltwiseFMAMod(result.a.data(), result.a.data(), constNum,
                nullptr, N, modulus, 1);
    intel::hexl::EltwiseFMAMod(result.b.data(), result.b.data(), constNum,
                nullptr, N, modulus, 1);
#endif
}


void queryGenerate(RlweCiphertext& cipher, Secret& queryKey, const int32_t row, const int32_t rlwesNum)
{
    // uint64_t length = queryKey.getLength();
    uint64_t length = N;
    // uint64_t modulus = queryKey.getModulus();
    std::vector<uint64_t> message(length);

    /**
    for (size_t i = 0; i < length; i++)
        message[i] = i * bsgsDelta;
    **/

    message[row * (length/rlwesNum)] = 1021937069741 * 100;
    // intel::hexl::NTT nttp(N, bsgsp);
    // nttp.ComputeInverse(message.data(), message.data(), 1, 1);

    // me \in (-p/2. p/2] * Delta
    // if me \in (p/2, p/2], (me - p) * Delta + bigMod
    /**/
    /**
    for (auto & me : message)
    {
        // VARIANT A
        me = round((long double)me * modulus / bsgsp); // TODO: it's same?
    }
    **/

    // store a in coefficient form
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = queryKey.getNTT();

    ntts.ComputeForward(message.data(), message.data(), 1, 1);
    // encrypt(cipher.b.data(), cipher.a.data(), queryKey, message.data());
    encrypt_special_rlwe(cipher.b.data(), cipher.a.data(), queryKey, message.data(), rlwesNum);

    // ntts.ComputeInverse(cipher.a.data(), cipher.a.data(), 1, 1);
    cipher.setIsNtt(true);
#endif
}


void encrypt_rgsw(RGSWCiphertext& result, Secret& secret, const std::vector<uint64_t>& message)
{    
    uint64_t length = secret.getLength();
    uint64_t modulus = secret.getModulus();

    // uint64_t mult_factor;
    std::vector<uint64_t> temp(length);

    // m(s, 1)
    if(!secret.isNttForm())
    {
        secret.toNttForm();
    }
    intel::hexl::NTT ntts = secret.getNTT();

    uint64_t PP = result.getBase();
    uint64_t BBg = result.getBg();
    int32_t ellnum = result.getEllnum();
    std::vector<uint64_t> mult_factor = powerOfBg(PP, BBg, ellnum);
    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        // faster by single barret reduce
        intel::hexl::EltwiseFMAMod(temp.data(), message.data(),
                mult_factor[i], nullptr, length, modulus, 1);

        ntts.ComputeForward(temp.data(), temp.data(), 1, 1);

        // rlwe encryption, these rlwes is in ntt form.
        // encrypt(m)
        encrypt(result.b[ellnum + i].data(), result.a[ellnum + i].data(), secret, temp.data());

        // encrypt(-ms)
        negate(temp.data(), length, modulus);

        intel::hexl::EltwiseMultMod(temp.data(), temp.data(), secret.getData().data(),
                            length, modulus, 1);
        encrypt(result.b[i].data(), result.a[i].data(), secret, temp.data());
#endif
    }
    result.setIsNtt(true);
}

/**
 * @brief 
 * 
 * @param cipher 
 * @param queryKey 
 * @param row 
 * @param col 
 * @param rlwesNum power of two
 * @param rgswsNum not always a power of two
 * @param ellgsw 
 * @param bg 
 */
void queryGenerate(RlweCiphertext& cipher, Secret& queryKey,
                   const int32_t row, const int32_t col,
                   const spiral_parameters& para)
                   // const int32_t rlwesNum, const int32_t rgswsNum,
                   // const int32_t ellgsw, const uint64_t bg, uint64_t spiral_p)
{
    // uint64_t length = queryKey.getLength();
    uint64_t length = N;
    uint64_t modulus = queryKey.getModulus();
    std::vector<uint64_t> message(length, 0);

    assert(para.rlwesNum + para.rgswsNum * para.ellgsw <= length); // TODO: if not, the client sends two query ciphertexts
    int32_t powof = 0x01 << (int32_t)ceil(log2(para.rgswsNum)); // e.g., 7 -> 8
    int32_t maxNums = para.rlwesNum >= para.rgswsNum * para.ellgsw ? para.rlwesNum : powof * para.ellgsw;

    // if (maxNums == rlwesNum)
    {
        // encode rlwe
        // uint64_t a = Delta;
        uint64_t a = modulus / para.spiral_p; // p = 2^8
        uint64_t b = intel::hexl::InverseMod(para.rlwesNum * 2, modulus);
        uint64_t encode;
        intel::hexl::EltwiseFMAMod(&encode, &a, b, nullptr, 1, modulus, 1);
        message[row * (length/para.rlwesNum)] = encode;

        // encode rgsw
        b = intel::hexl::InverseMod(powof * para.ellgsw * 2, modulus);
        std::vector<uint64_t> mult_factor = powerOfBg(0, para.bggsw, para.ellgsw);
        int32_t interval = length / (maxNums * 2);

        for (size_t i = 0; i < para.rgswsNum; i++)
        {
            int32_t encode_bit = (col >> i) & 0x01;
            for (size_t j = 0; j < para.ellgsw; j++)
            {
                int32_t index = para.ellgsw * i + j;
                // a = pow_mod(bg, j, modulus);
                a = mult_factor[j];
                intel::hexl::EltwiseFMAMod(&encode, &a, b, nullptr, 1, modulus, 1);
                // message[index * (length/(powof * ellgsw)) + length / (maxNums * 2)] = encode * encode_bit;
                message[index * (length/(powof * para.ellgsw)) + interval] = encode * encode_bit;
            }
        }
    }
    // for (size_t i = 1; i < length/2; i++) message[i] = 0;
    // for (size_t i = length/2 + 1; i < length; i++) message[i] = 0;

    // store a in coefficient form
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = queryKey.getNTT();

    ntts.ComputeForward(message.data(), message.data(), 1, 1);
    encrypt(cipher.b.data(), cipher.a.data(), queryKey, message.data());

    cipher.setIsNtt(true);
#endif
}

/**
 * @brief rns variant of queryGenerate
 * 
 */
void queryGenerate(std::vector<RlweCiphertext>& cipher, Secret& queryKey,
                   const int32_t row, const int32_t col,
                   const spiral_parameters& para)
{
    // uint64_t length = queryKey.getLength();
    uint64_t length = N;
    uint64_t modulus = queryKey.getModulus();
    std::vector<uint64_t> message(length, 0);

    assert(para.rlwesNum + para.rgswsNum * para.ellgsw <= length); // TODO: if not, the client sends two query ciphertexts
    int32_t powof = 0x01 << (int32_t)ceil(log2(para.rgswsNum)); // e.g., 7 -> 8
    int32_t maxNums = para.rlwesNum >= para.rgswsNum * para.ellgsw ? para.rlwesNum : powof * para.ellgsw;

    {
        // encode rlwe
        // uint64_t a = Delta;
        uint64_t a = modulus / para.spiral_p; // p = 2^8
        uint64_t b = intel::hexl::InverseMod(para.rlwesNum * 2, modulus);
        uint64_t encode;
        intel::hexl::EltwiseFMAMod(&encode, &a, b, nullptr, 1, modulus, 1);
        message[row * (length/para.rlwesNum)] = encode;

        // encode rgsw
        b = intel::hexl::InverseMod(powof * para.ellgsw * 2, modulus);
        std::vector<uint64_t> mult_factor = powerOfBg(0, para.bggsw, para.ellgsw);
        int32_t interval = length / (maxNums * 2);

        for (size_t i = 0; i < para.rgswsNum; i++)
        {
            int32_t encode_bit = (col >> i) & 0x01;
            for (size_t j = 0; j < para.ellgsw; j++)
            {
                int32_t index = para.ellgsw * i + j;
                // a = pow_mod(bg, j, modulus);
                a = mult_factor[j];
                intel::hexl::EltwiseFMAMod(&encode, &a, b, nullptr, 1, modulus, 1);
                message[index * (length/(powof * para.ellgsw)) + interval] = encode * encode_bit;
            }
        }
    }

    // store a in coefficient form
#ifdef INTEL_HEXL
    encrypt_rns_bsgs(cipher, queryKey, message);
#endif
}

/**
 * @brief 
 * 
 * @param cipher 
 * @param queryKey 
 * @param row 
 * @param col 
 * @param rlwesNum power of two
 * @param rgswsNum not always a power of two
 * @param ellgsw 
 * @param bg 
 */
void queryGenerate(std::vector<RlweCiphertext>& query_rlwe, std::vector<RGSWCiphertext>& query_rgsw,
                   Secret& queryKey,
                   const int32_t row, const int32_t col,
                   const spiral_parameters& para)
{
    uint64_t length = queryKey.getLength();
    uint64_t modulus = queryKey.getModulus();
    std::vector<uint64_t> message(length, 0);
    std::vector<uint64_t> temp(length, 0);

    uint64_t encode = modulus / para.spiral_p; 

    // encrypt rlwe ciphertexts
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = queryKey.getNTT();
    for (size_t i = 0; i < para.rlwesNum; i++)
    {
        message[0] = i == row ? encode : 0;

        ntts.ComputeForward(temp.data(), message.data(), 1, 1);
        encrypt(query_rlwe[i].b.data(), query_rlwe[i].a.data(), queryKey, temp.data());

        query_rlwe[i].setIsNtt(true);
    }
#endif

    // encrypt rgsw ciphertexts
    for (size_t i = 0; i < para.rgswsNum; i++)
    {
        // int32_t encode_bit = (col >> i) & 0x01;
        // message[0] = i == col ? 1 : 0;
        message[0] = (col >> i) & 0x01;
        encrypt_rgsw(query_rgsw[i], queryKey, message);
    }
}

/**
 * @brief the variant of packingLWEs algorithm. We packing RLWE ciphertexts instead of only LWE ciphertexts.
 *  Refer to Algorithm 2 in [https://eprint.iacr.org/2020/015].
 * @param result
 * @param cipher1  the input RLWE ciphertext 2
 * @param curl 
 */
void expandTwoRLWEs(RlweCiphertext& result, RlweCiphertext& cipher,
                    const uint64_t twol, const uint64_t rlweNums,
                    const AutoKey& autokey)
{   
    int32_t XNdiv2l = N / twol;
    // int32_t XNdiv2l = rlweNums / twol;

    cipher = result; // copy it
    cipher.multNegConst(XNdiv2l); // x^XNdiv2l

    // int32_t index = XNdiv2l + 1;
    int32_t index = twol + 1;

    RlweCiphertext temp = result;
    evalAuto(temp, index, autokey);
    result.addAndEqual(temp);

    temp = cipher;
    evalAuto(temp, index, autokey);
    cipher.addAndEqual(temp);
}

void expandTwoRLWEs(RlweCiphertext& result, RlweCiphertext& result_rns,
                    RlweCiphertext& cipher, RlweCiphertext& cipher_rns,
                    const uint64_t twol, const uint64_t rlweNums,
                    const AutoKeyRNS& autokey)
{   
    int32_t XNdiv2l = N / twol;
    // int32_t XNdiv2l = rlweNums / twol;

    cipher = result; // copy it
    cipher.multNegConst(XNdiv2l); // x^XNdiv2l
    
    cipher_rns = result_rns;
    cipher_rns.multNegConst(XNdiv2l);

    int32_t index = twol + 1;

    RlweCiphertext temp = result;
    RlweCiphertext temp2 = result_rns;
    evalAutoRNS(temp, temp2, index, autokey);
    addRnsCiphertext(result, result_rns, temp, temp2);

    temp = cipher;
    temp2 = cipher_rns;
    evalAutoRNS(temp, temp2, index, autokey);
    addRnsCiphertext(cipher, cipher_rns, temp, temp2);
}

/**
 * @brief the one RLWE ciphertext variant of expandTwoRLWEs
 *        Sometimes, the number of output rgsw ciphetexts is not power-of-two.
 * 
 * @param result 
 * @param twol 
 * @param rlweNums 
 * @param autokey 
 */
void expandOneRLWEs(RlweCiphertext& result,
                    const uint64_t twol, const uint64_t rlweNums,
                    const AutoKey& autokey)
{
    int32_t index = twol + 1;

    RlweCiphertext temp = result;
    evalAuto(temp, index, autokey);
    result.addAndEqual(temp);
}

void expandOneRLWEs(RlweCiphertext& result, RlweCiphertext& result_rns,
                    const uint64_t twol, const uint64_t rlweNums,
                    const AutoKeyRNS& autokey)
{
    int32_t index = twol + 1;

    RlweCiphertext temp = result;
    RlweCiphertext temp2 = result_rns;
    evalAutoRNS(temp, temp2, index, autokey);
    addRnsCiphertext(result, result_rns, temp, temp2);
}

/**
 * @brief 
 * 
 * @param result response of answer
 * @param rlwes the requried packing rlwe ciphertexts.
 *      The size should be power of two. lwes should be store in coefficient form.
 */
void queryExpand(std::vector<RlweCiphertext>& result, RlweCiphertext& input,
                 const AutoKey& autokey, const int32_t rlwesNum)
{    
    uint64_t output_size = result.size(); // power of two
    assert(output_size == rlwesNum);
    
    uint64_t modulus = input.getModulus();
    // assert(output_size == rlwesNum, "The out expansion ciphetexts size should be equal to rlwesNum");

    // the input lwes should be in coeffcient form
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = autokey.getNTT();
    if (input.getIsNtt())
    {
        ntts.ComputeInverse(input.a.data(), input.a.data(), 1, 1);
        ntts.ComputeInverse(input.b.data(), input.b.data(), 1, 1);
        input.setIsNtt(false);
#endif
    }

    result[0] = input;

    // for (size_t i = rlwesNum/2; i > 0; i >>= 1)
    for (size_t i = 1; i < rlwesNum; i <<= 1)
    {
        for (size_t j = 0; j < i; j++)
        {
            expandTwoRLWEs(result[j], result[i + j], rlwesNum / i, rlwesNum, autokey);
        }
    }

#ifdef INTEL_HEXL
    uint64_t b = intel::hexl::InverseMod(rlwesNum, modulus);
#endif

    // multiply the inverse of rlwesNum
    // remove term 2^l
    for (size_t i = 0; i < rlwesNum; i++)
    {
        multConst(result[i], b);
    }
    /*
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "evalTrNn costs " << glapsed.count() << " us." << std::endl;
    */
}


/**
 * @brief 
 * 
 * @param result response of answer
 * @param rlwes the requried packing rlwe ciphertexts.
 *      The size should be power of two. lwes should be store in coefficient form.
 */
void queryExpand(std::vector<RlweCiphertext>& expansion_rlwes, std::vector<RlweCiphertext>& expansion_rgsws,
                 RlweCiphertext& input, const AutoKey& autokey,
                 const int32_t rlwesNum, const int32_t rgswsNum, const int32_t ellgsw)
{    
    int32_t rlwes_size = expansion_rlwes.size(); // power of two
    assert(rlwes_size == rlwesNum); // The expansion rlwes ciphetexts size should be equal to rlwesNum
    int32_t rgsws_size = expansion_rgsws.size();
    assert(rgsws_size == rgswsNum * ellgsw); // The expansion rgsws ciphetexts size should be equal to rlwesNum
    
    // uint64_t modulus = input.getModulus();

    // the input lwes should be in coeffcient form
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = autokey.getNTT();
    if (input.getIsNtt())
    {
        ntts.ComputeInverse(input.a.data(), input.a.data(), 1, 1);
        ntts.ComputeInverse(input.b.data(), input.b.data(), 1, 1);
        input.setIsNtt(false);
#endif
    }
    expansion_rlwes[0] = input;

    int32_t powof = 0x01 << (int32_t)ceil(log2(rgswsNum)); // e.g., 7 -> 8
    int32_t maxNums = rlwesNum >= rgswsNum * ellgsw ? rlwesNum : powof * ellgsw;

    expandTwoRLWEs(expansion_rlwes[0], expansion_rgsws[0], maxNums * 2, maxNums * 2, autokey);

    // expand rlwe ciphertexts
    for (size_t i = 1; i < rlwesNum; i <<= 1)
    {
        for (size_t j = 0; j < i; j++)
        {
            expandTwoRLWEs(expansion_rlwes[j], expansion_rlwes[i + j], rlwesNum / i, rlwesNum, autokey);
        }
    }

    // expand rgsw ciphertexts
    // TODO: check: support rgswsNum which is not power of two
    for (size_t i = 1; i < rgswsNum * ellgsw; i <<= 1)
    {
        for (size_t j = 0; j < i; j++)
        {
            if (i + j >= rgswsNum * ellgsw)
            {
                // if exceed the size, just perform for the first ciphertext.
                expandOneRLWEs(expansion_rgsws[j], powof * ellgsw / i, powof * ellgsw, autokey);
            } else {
                expandTwoRLWEs(expansion_rgsws[j], expansion_rgsws[i + j], powof * ellgsw / i, powof * ellgsw, autokey);
            }
        }
    }
}

void modSwitch(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input,
               const AutoKeyRNS& autokey)
{
    // rns modswitch:
    //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
    int32_t length = autokey.getLength();
    uint64_t modulus1 = autokey.getModulus1();
    // uint64_t modulus2 = autokey.getModulus2();

    intel::hexl::NTT ntts1 = autokey.getNTT1();
    intel::hexl::NTT ntts2 = autokey.getNTT2();

    // std::vector<uint64_t> temp1(N);
    // std::vector<uint64_t> temp2(N);

    if (input[0].getIsNtt())
    {
        for (size_t i = 0; i < input.size(); i++)
        {
            ntts1.ComputeInverse(result[i].a.data(), result[i].a.data(), 1, 1);
            ntts1.ComputeInverse(result[i].b.data(), result[i].b.data(), 1, 1);

            ntts2.ComputeInverse(input[i].a.data(), input[i].a.data(), 1, 1);
            ntts2.ComputeInverse(input[i].b.data(), input[i].b.data(), 1, 1);
        }
    }

    for (size_t i = 0; i < input.size(); i++)
    {
#ifdef INTEL_HEXL
        //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
        intel::hexl::EltwiseSubMod(result[i].a.data(), result[i].a.data(), input[i].a.data(),
                            length, modulus1);
        intel::hexl::EltwiseFMAMod(result[i].a.data(), result[i].a.data(), bsModInv,
                            nullptr, length, modulus1, 1);
        ntts1.ComputeForward(result[i].a.data(), result[i].a.data(), 1, 1);

        // handle b
        intel::hexl::EltwiseSubMod(result[i].b.data(), result[i].b.data(), input[i].b.data(),
                            length, modulus1);
        intel::hexl::EltwiseFMAMod(result[i].b.data(), result[i].b.data(), bsModInv,
                            nullptr, length, modulus1, 1);
        ntts1.ComputeForward(result[i].b.data(), result[i].b.data(), 1, 1);

        result[i].setIsNtt(true);
#endif
    }
}

void modSwitch(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& input)
{
    // rns modswitch:
    //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
    assert(!input[0].getIsNtt());
    assert(!result[0].getIsNtt());

    int32_t length = result[0].getLength();
    uint64_t modulus1 = result[0].getModulus();
    // uint64_t modulus2 = input[0].getModulus();

    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

    for (size_t i = 0; i < input.size(); i++)
    {
#ifdef INTEL_HEXL
        //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
        intel::hexl::EltwiseSubMod(result[i].a.data(), result[i].a.data(), input[i].a.data(),
                            length, modulus1);
        intel::hexl::EltwiseFMAMod(result[i].a.data(), result[i].a.data(), bsModInv,
                        nullptr, length, modulus1, 1);
        // ntts1.ComputeForward(result[i].a.data(), result[i].a.data(), 1, 1);

        // handle b
        intel::hexl::EltwiseSubMod(result[i].b.data(), result[i].b.data(), input[i].b.data(),
                            length, modulus1);
        intel::hexl::EltwiseFMAMod(result[i].b.data(), result[i].b.data(), bsModInv,
                        nullptr, length, modulus1, 1);
        // ntts1.ComputeForward(result[i].b.data(), result[i].b.data(), 1, 1);

        // result[i].setIsNtt();
#endif
    }
}

void modSwitchToGSW(RGSWCiphertext& result,
                    std::vector<RlweCiphertext>::iterator input, std::vector<RlweCiphertext>::iterator input_rns,
                    const RNSRGSWCiphertext& convertKey, const int32_t ellgsw)
{
    // rns modswitch:
    //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)

    int32_t length = input->getLength();
    uint64_t modulus1 = input->getModulus();

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts1 = convertKey.getNTT();
    intel::hexl::NTT ntts2 = convertKey.getNTT2();

    if (input[0].getIsNtt())
    {
        for (size_t i = 0; i < ellgsw; i++)
        {
            ntts1.ComputeInverse(input[i].a.data(), input[i].a.data(), 1, 1);
            ntts1.ComputeInverse(input[i].b.data(), input[i].b.data(), 1, 1);
            ntts2.ComputeInverse(input_rns[i].a.data(), input_rns[i].a.data(), 1, 1);
            ntts2.ComputeInverse(input_rns[i].b.data(), input_rns[i].b.data(), 1, 1);
        }  
    }

    for (size_t i = 0; i < ellgsw; i++)
    {

        //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
        intel::hexl::EltwiseSubMod(result.a[ellgsw + i].data(), input[i].a.data(),
                        input_rns[i].a.data(), length, modulus1);
        intel::hexl::EltwiseFMAMod(result.a[ellgsw + i].data(), result.a[ellgsw + i].data(), bsModInv,
                        nullptr, length, modulus1, 1);
        ntts1.ComputeForward(result.a[ellgsw + i].data(), result.a[ellgsw + i].data(), 1, 1);

        // handle b
        intel::hexl::EltwiseSubMod(result.b[ellgsw + i].data(), input[i].b.data(),
                        input_rns[i].b.data(), length, modulus1);
        intel::hexl::EltwiseFMAMod(result.b[ellgsw + i].data(), result.b[ellgsw + i].data(), bsModInv,
                        nullptr, length, modulus1, 1);
        ntts1.ComputeForward(result.b[ellgsw + i].data(), result.b[ellgsw + i].data(), 1, 1);

        // result[i].setIsNtt();
    }
#endif
}

/**
 * @brief rns variant of queryExpand
 * 
 */
/**/
void queryExpand(std::vector<RlweCiphertext>& expansion_rlwes, std::vector<RlweCiphertext>& expansion_rgsws,
                 std::vector<RlweCiphertext>& input, const AutoKeyRNS& autokey,
                 const int32_t rlwesNum, const int32_t rgswsNum, const int32_t ellgsw)
{    
    int32_t rlwes_size = expansion_rlwes.size(); // power of two
    assert(rlwes_size == rlwesNum); // The expansion rlwes ciphetexts size should be equal to rlwesNum
    int32_t rgsws_size = expansion_rgsws.size();
    assert(rgsws_size == 2 * rgswsNum * ellgsw); // The expansion rgsws ciphetexts size should be equal to 2 * rlwesNum * ellgsw
    
    // uint64_t modulus = input.getModulus();
    int32_t length = autokey.getLength();
    uint64_t rnsModulus = autokey.getModulus2();
    std::vector<RlweCiphertext> expansion_rlwes_rns(rlwes_size, RlweCiphertext(length, rnsModulus));

    // the input lwes should be in coeffcient form
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts = autokey.getNTT1();
    if (input[0].getIsNtt())
    {
        ntts.ComputeInverse(input[0].a.data(), input[0].a.data(), 1, 1);
        ntts.ComputeInverse(input[0].b.data(), input[0].b.data(), 1, 1);
        input[0].setIsNtt(false);
    }
    intel::hexl::NTT ntts2 = autokey.getNTT2();
    if (input[1].getIsNtt())
    {
        ntts2.ComputeInverse(input[1].a.data(), input[1].a.data(), 1, 1);
        ntts2.ComputeInverse(input[1].b.data(), input[1].b.data(), 1, 1);
        input[1].setIsNtt(false);
    }
#endif
    expansion_rlwes[0] = input[0];
    expansion_rlwes_rns[0] = input[1];

    int32_t powof = 0x01 << (int32_t)ceil(log2(rgswsNum)); // e.g., 7 -> 8
    int32_t maxNums = rlwesNum >= rgswsNum * ellgsw ? rlwesNum : powof * ellgsw;

    int32_t rgsws = rgswsNum * ellgsw;

#ifdef SPIRAL_TIME_COLLECT
    auto start  = std::chrono::high_resolution_clock::now();
#endif
    expandTwoRLWEs(expansion_rlwes[0], expansion_rlwes_rns[0],
                   expansion_rgsws[0], expansion_rgsws[rgsws], maxNums * 2, maxNums * 2, autokey);

    // expand rlwe ciphertexts
    for (size_t i = 1; i < rlwesNum; i <<= 1)
    {
        for (size_t j = 0; j < i; j++)
        {
            expandTwoRLWEs(expansion_rlwes[j], expansion_rlwes_rns[j],
                           expansion_rlwes[i + j], expansion_rlwes_rns[i + j], rlwesNum / i, rlwesNum, autokey);
        }
    }

#ifdef SPIRAL_TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " rlwe expansion costs " << glapsed.count() << " us." << std::endl;
#endif

    // modswitch
    modSwitch(expansion_rlwes, expansion_rlwes_rns, autokey);

    // expand rgsw ciphertexts
    for (size_t i = 1; i < rgsws; i <<= 1)
    {
        for (size_t j = 0; j < i; j++)
        {
            if (i + j >= rgsws)
            {
                // if exceed the size, just perform for the first ciphertext.
                expandOneRLWEs(expansion_rgsws[j], expansion_rgsws[rgsws + j], powof * ellgsw / i, powof * ellgsw, autokey);
            } else {
                expandTwoRLWEs(expansion_rgsws[j], expansion_rgsws[rgsws + j],
                               expansion_rgsws[i + j], expansion_rgsws[rgsws + i + j],
                               powof * ellgsw / i, powof * ellgsw, autokey);
            }
        }
    }
}

void keyGen(RGSWCiphertext& result, Secret& secret)
{    
    uint64_t length = secret.getLength();
    uint64_t modulus = secret.getModulus();
    int32_t ellnum = result.getEllnum();

    std::vector<uint64_t> message(length);
    if (secret.isNttForm())
    {
        secret.toCoeffForm();
    }
    // why the following code doesn't work?
    // copy(secret.getData().begin(), secret.getData().end(), message.begin());
    // encrypt -s to convertKey
    for (size_t i = 0; i < length; i++)
    {
        // the value `modulus' can be reduced in the following steps?
        message[i] = modulus - secret.getData(i);
    }

    // uint64_t mult_factor;
    std::vector<uint64_t> temp(length);

    // m(s, 1)
    if(!secret.isNttForm())
    {
        secret.toNttForm();
    }
    intel::hexl::NTT ntts = secret.getNTT();

    std::vector<uint64_t> mult_factor = powerOfBg(result.getBase(), result.getBg(), result.getEllnum());
    for (size_t i = 0; i < ellnum; i++)
    {
#ifdef INTEL_HEXL
        // faster by single barret reduce
        intel::hexl::EltwiseFMAMod(temp.data(), message.data(),
                mult_factor[i], nullptr, length, modulus, 1);

        ntts.ComputeForward(temp.data(), temp.data(), 1, 1);

        // rlwe encryption, these rlwes is in ntt form.
        // encrypt(m)
        encrypt(result.b[ellnum + i].data(), result.a[ellnum + i].data(), secret, temp.data());

        // encrypt(-ms)
        // intel::hexl::EltwiseSubMod(temp.data(), this->b[ellnum - 1].data(), temp.data(),
        //        this->length, this->modulus);
        negate(temp.data(), length, modulus); // negating can be performed in its ntt form

        intel::hexl::EltwiseMultMod(temp.data(), temp.data(), secret.getData().data(),
                length,modulus, 1);
        encrypt(result.b[i].data(), result.a[i].data(), secret, temp.data());
#endif
    }
    result.setIsNtt(true);
}

/**
 * @brief 
 * 
 * @param result coefficient forms
 * @param firstDim  form
 * @param query ntt form
 */
void externalProduct(std::vector<uint64_t>& result_a, std::vector<uint64_t>& result_b,
                     RlweCiphertext& firstDim, const RGSWCiphertext& query)
{
    assert(query.getIsNtt());

    uint64_t ellnum = query.getEllnum();
    uint64_t length = query.getLength();
    uint64_t modulus = query.getModulus();

#ifdef INTEL_HEXL
    // keyswitch
    std::vector<std::vector<uint64_t> > decFirstDim(2 * ellnum, std::vector<uint64_t>(length, 0));
    std::vector<uint64_t> temp(length, 0);
    intel::hexl::NTT ntts = query.getNTT();

    if (firstDim.getIsNtt())
    {
        ntts.ComputeInverse(firstDim.a.data(), firstDim.a.data(), 1, 1);
        ntts.ComputeInverse(firstDim.b.data(), firstDim.b.data(), 1, 1);
        firstDim.setIsNtt(false);
    }

    // g^-1(a, b) * key
    decompose_rlwe(decFirstDim, firstDim.a, firstDim.b, ellnum, query.getBase(), query.getBg(), modulus);
    for (size_t i = 0; i <  2 * ellnum; i++)
    {
        ntts.ComputeForward(decFirstDim[i].data(), decFirstDim[i].data(), 1, 1);

        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    query.a[i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result_a.data(), result_a.data(), temp.data(), length, modulus);

    }

    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    query.b[i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result_b.data(), result_b.data(), temp.data(), length, modulus);
    }

    // return to ntt form
    // ntts.ComputeInverse(result_a.data(), result_a.data(), 1, 1);
    // ntts.ComputeInverse(result_b.data(), result_b.data(), 1, 1);
    // result.setIsNtt(true);
#endif
}

/**
 * @brief 
 * rns variant of externalProduct
 */
/**
void externalProduct(RlweCiphertext& result, RlweCiphertext& result_rns,
                     RlweCiphertext& firstDim, RlweCiphertext& firstDim_rns, const RNSRGSWCiphertext& query)
**/
/**
 * @brief external product and modulus switching
 */
void externalProductSwitch(std::vector<uint64_t>& result_a, std::vector<uint64_t>& result_b,
                           RlweCiphertext& firstDim, RlweCiphertext& firstDim_rns, const RNSRGSWCiphertext& query)
{
    assert(query.getIsNtt());
    int32_t ellnum = query.getEllnum();
    int32_t length = query.getLength();
    uint64_t modulus = query.getModulus();
    uint64_t modulus2 = query.getSecModulus();

    RlweCiphertext result_rns(length, modulus2);

#ifdef INTEL_HEXL
    // keyswitch
    std::vector<uint64_t> temp(length, 0);
    intel::hexl::NTT ntts = query.getNTT();
    intel::hexl::NTT ntts2 = query.getNTT2();

    if (firstDim.getIsNtt())
    {
        ntts.ComputeInverse(firstDim.a.data(), firstDim.a.data(), 1, 1);
        ntts.ComputeInverse(firstDim.b.data(), firstDim.b.data(), 1, 1);
        firstDim.setIsNtt(false);
    }

    if (firstDim_rns.getIsNtt())
    {
        ntts2.ComputeInverse(firstDim_rns.a.data(), firstDim_rns.a.data(), 1, 1);
        ntts2.ComputeInverse(firstDim_rns.b.data(), firstDim_rns.b.data(), 1, 1);
        firstDim_rns.setIsNtt(false);
    }

    // g^-1(a, b) * key
    std::vector<std::vector<uint64_t> > dec_a(ellnum, std::vector<uint64_t>(length, 0));
    std::vector<std::vector<uint64_t> > dec_a_rns(ellnum, std::vector<uint64_t>(length, 0));
    // decompose_crt(dec_a, dec_a_rns, firstDim.a, firstDim_rns.a, ellnum, query.getBase(), query.getBg());
    decompose_bsgs(dec_a, dec_a_rns, firstDim.a, firstDim_rns.a, ellnum, query.getBase(), query.getBg());

    for (size_t i = 0; i < ellnum; i++)
    {
        ntts.ComputeForward(dec_a[i].data(), dec_a[i].data(), 1, 1);

        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    query.a[i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result_a.data(), result_a.data(), temp.data(), length, modulus);
    }
    for (size_t i = 0; i < ellnum; i++)
    {
        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    query.b[i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result_b.data(), result_b.data(), temp.data(), length, modulus);
    }

    // another part of rns
    for (size_t i = 0; i < ellnum; i++)
    {
        ntts2.ComputeForward(dec_a_rns[i].data(), dec_a_rns[i].data(), 1, 1);

        intel::hexl::EltwiseMultMod(temp.data(), dec_a_rns[i].data(),
                    query.arns[i].data(), length, modulus2, 1);

        intel::hexl::EltwiseAddMod(result_rns.a.data(), result_rns.a.data(), temp.data(), length, modulus2);
    }
    for (size_t i = 0; i < ellnum; i++)
    {
        intel::hexl::EltwiseMultMod(temp.data(), dec_a_rns[i].data(),
                    query.brns[i].data(), length, modulus2, 1);

        intel::hexl::EltwiseAddMod(result_rns.b.data(), result_rns.b.data(), temp.data(), length, modulus2);
    }

    decompose_bsgs(dec_a, dec_a_rns, firstDim.b, firstDim_rns.b, ellnum, query.getBase(), query.getBg());
    for (size_t i = 0; i < ellnum; i++)
    {
        ntts.ComputeForward(dec_a[i].data(), dec_a[i].data(), 1, 1);

        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    query.a[ellnum + i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result_a.data(), result_a.data(), temp.data(), length, modulus);
    }
    for (size_t i = 0; i < ellnum; i++)
    {
        intel::hexl::EltwiseMultMod(temp.data(), dec_a[i].data(),
                    query.b[ellnum + i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result_b.data(), result_b.data(), temp.data(), length, modulus);
    }

    // another part of rns
    for (size_t i = 0; i < ellnum; i++)
    {
        ntts2.ComputeForward(dec_a_rns[i].data(), dec_a_rns[i].data(), 1, 1);

        intel::hexl::EltwiseMultMod(temp.data(), dec_a_rns[i].data(),
                    query.arns[ellnum + i].data(), length, modulus2, 1);

        intel::hexl::EltwiseAddMod(result_rns.a.data(), result_rns.a.data(), temp.data(), length, modulus2);
    }
    for (size_t i = 0; i < ellnum; i++)
    {
        intel::hexl::EltwiseMultMod(temp.data(), dec_a_rns[i].data(),
                    query.brns[ellnum + i].data(), length, modulus2, 1);

        intel::hexl::EltwiseAddMod(result_rns.b.data(), result_rns.b.data(), temp.data(), length, modulus2);
    }

    // return to coefficient form
    ntts.ComputeInverse(result_a.data(), result_a.data(), 1, 1);
    ntts.ComputeInverse(result_b.data(), result_b.data(), 1, 1);
    ntts2.ComputeInverse(result_rns.a.data(), result_rns.a.data(), 1, 1);
    ntts2.ComputeInverse(result_rns.b.data(), result_rns.b.data(), 1, 1);

    // modulus switching
    //  q_2^-1 * ((b1, a1) - (b2, a2)) (mod q_1)
    intel::hexl::EltwiseSubMod(result_a.data(), result_a.data(), result_rns.a.data(),
                        length, modulus);
    intel::hexl::EltwiseFMAMod(result_a.data(), result_a.data(), bsModInv,
                        nullptr, length, modulus, 1);
    ntts.ComputeForward(result_a.data(), result_a.data(), 1, 1);

    // handle b
    intel::hexl::EltwiseSubMod(result_b.data(), result_b.data(), result_rns.b.data(),
                        length, modulus);
    intel::hexl::EltwiseFMAMod(result_b.data(), result_b.data(), bsModInv,
                        nullptr, length, modulus, 1);
    ntts.ComputeForward(result_b.data(), result_b.data(), 1, 1);
    // result.setIsNtt(true);
#endif
}

void regevToGSW(RGSWCiphertext& result, std::vector<RlweCiphertext>& input,
                const RGSWCiphertext& convertKey, const int32_t ellgsw, const int32_t start_index)
{
    // int32_t ellgsw = input.size();

    if (input[0].getIsNtt()) // we assume that all getIsNtt() is equal to that of input[0]
    {
        // if input ciphertext is ntt form, copy first
        for (size_t i = 0; i < ellgsw; i++)
        {
            copy(input[i].a.begin(), input[i].a.end(), result.a[ellgsw + i].begin());
            copy(input[i].b.begin(), input[i].b.end(), result.b[ellgsw + i].begin());
            externalProduct(result.a[i], result.b[i], input[i], convertKey);
        }
    } else {
        // if input ciphertext is not ntt form, 
        intel::hexl::NTT ntts = convertKey.getNTT();
        for (size_t i = 0; i < ellgsw; i++)
        {
            ntts.ComputeForward(result.a[ellgsw + i].data(), input[i].a.data(), 1, 1);
            ntts.ComputeForward(result.a[ellgsw + i].data(), input[i].b.data(), 1, 1);
            externalProduct(result.a[i], result.b[i], input[i], convertKey);
        }
    }
    // the RGSW ciphertext is store in coefficient form
    result.setIsNtt(true);
}

void regevToGSW(RGSWCiphertext& result, std::vector<RlweCiphertext>::iterator input,
                std::vector<RlweCiphertext>::iterator input_end,
                const RGSWCiphertext& convertKey, const int32_t ellgsw)
{
    assert(input_end - input == ellgsw);
    
    if (input[0].getIsNtt()) // we assume that all getIsNtt() is equal to that of input[0]
    {
        // if input ciphertext is ntt form, copy first
        for (size_t i = 0; i < ellgsw; i++)
        {
            copy(input[i].a.begin(), input[i].a.end(), result.a[ellgsw + i].begin());
            copy(input[i].b.begin(), input[i].b.end(), result.b[ellgsw + i].begin());
            externalProduct(result.a[i], result.b[i], input[i], convertKey);
        }
    } else {
        // if input ciphertext is not ntt form, 
        intel::hexl::NTT ntts = convertKey.getNTT();
        for (size_t i = 0; i < ellgsw; i++)
        {
            ntts.ComputeForward(result.a[ellgsw + i].data(), input[i].a.data(), 1, 1);
            ntts.ComputeForward(result.b[ellgsw + i].data(), input[i].b.data(), 1, 1);
            externalProduct(result.a[i], result.b[i], input[i], convertKey);
        }
    }
    // the RGSW ciphertext is store in coefficient form
    result.setIsNtt(true);
}

// the output 
void regevToGSW(RGSWCiphertext& result, std::vector<RlweCiphertext>::iterator input,
                std::vector<RlweCiphertext>::iterator input_end,
                std::vector<RlweCiphertext>::iterator input_rns,
                std::vector<RlweCiphertext>::iterator input_rns_end,
                const RNSRGSWCiphertext& convertKey, const int32_t ellgsw)
{
    assert(input_end - input == ellgsw);

    modSwitchToGSW(result, input, input_rns, convertKey, ellgsw);
    for (size_t i = 0; i < ellgsw; i++)
    {
        externalProductSwitch(result.a[i], result.b[i], input[i], input_rns[i], convertKey);
    }

    // the RGSW ciphertext is store in ntt form
    result.setIsNtt(true);
}

/**
 * @brief 
 * 
 * @param result coefficient forms
 * @param firstDim  form
 * @param query ntt form
 */
void externalProduct(RlweCiphertext& result, const RGSWCiphertext& query)
{
    assert(query.getIsNtt());

    uint64_t ellnum = query.getEllnum();
    uint64_t length = query.getLength();
    uint64_t modulus = query.getModulus();

#ifdef INTEL_HEXL
    // keyswitch
    std::vector<std::vector<uint64_t> > decFirstDim(2 * ellnum, std::vector<uint64_t>(length, 0));
    std::vector<uint64_t> temp(length, 0);
    intel::hexl::NTT ntts = query.getNTT();

    if (result.getIsNtt())
    {
        ntts.ComputeInverse(result.a.data(), result.a.data(), 1, 1);
        ntts.ComputeInverse(result.b.data(), result.b.data(), 1, 1);
        result.setIsNtt(false);
    }

    // g^-1(a, b) * key
    decompose_rlwe(decFirstDim, result.a, result.b, ellnum, query.getBase(), query.getBg(), modulus);
    
    copy(temp.begin(), temp.end(), result.a.begin());
    copy(temp.begin(), temp.end(), result.b.begin());
    for (size_t i = 0; i <  2 * ellnum; i++)
    {
        ntts.ComputeForward(decFirstDim[i].data(), decFirstDim[i].data(), 1, 1);

        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    query.a[i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result.a.data(), result.a.data(), temp.data(), length, modulus);

    }

    for (size_t i = 0; i < 2 * ellnum; i++)
    {
        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    query.b[i].data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), temp.data(), length, modulus);
    }

    // return to coeffs form
    ntts.ComputeInverse(result.a.data(), result.a.data(), 1, 1);
    ntts.ComputeInverse(result.b.data(), result.b.data(), 1, 1);
    result.setIsNtt(false);
#endif
}

/**
 * @brief the ComplementingGSWEncodings
 * 
 * @param result coefficient forms
 * @param firstDim  form
 * @param query ntt form
 */
void externalProductComp(RlweCiphertext& result, const RGSWCiphertext& query)
{
    assert(query.getIsNtt());

    uint64_t ellnum = query.getEllnum();
    uint64_t length = query.getLength();
    uint64_t modulus = query.getModulus();

#ifdef INTEL_HEXL
    // keyswitch
    std::vector<std::vector<uint64_t> > decFirstDim(2 * ellnum, std::vector<uint64_t>(length, 0));
    std::vector<uint64_t> temp(length, 0);
    std::vector<uint64_t> zeros(length, 0);
    intel::hexl::NTT ntts = query.getNTT();

    if (result.getIsNtt())
    {
        ntts.ComputeInverse(result.a.data(), result.a.data(), 1, 1);
        ntts.ComputeInverse(result.b.data(), result.b.data(), 1, 1);
        result.setIsNtt(false);
    }

    // g^-1(a, b) * key
    decompose_rlwe(decFirstDim, result.a, result.b, ellnum, query.getBase(), query.getBg(), modulus);
    std::vector<uint64_t> mult_factor = powerOfBg(query.getBase(), query.getBg(), ellnum);

    copy(zeros.begin(), zeros.end(), result.a.begin());
    copy(zeros.begin(), zeros.end(), result.b.begin());
    // for (size_t i = 0; i <  2 * ellnum; i++)
    for (size_t i = 0; i < ellnum; i++)
    {
        ntts.ComputeForward(decFirstDim[i].data(), decFirstDim[i].data(), 1, 1);
        // G - query
        zeros[0] = mult_factor[i];
        ntts.ComputeForward(temp.data(), zeros.data(), 1, 1);
        intel::hexl::EltwiseSubMod(temp.data(), temp.data(), query.a[i].data(), length, modulus);
        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    temp.data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result.a.data(), result.a.data(), temp.data(), length, modulus);

    } for (size_t i = ellnum; i <  2 * ellnum; i++)
    {
        ntts.ComputeForward(decFirstDim[i].data(), decFirstDim[i].data(), 1, 1);

        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    query.a[i].data(), length, modulus, 1);

        // note: substrct instead of addition
        intel::hexl::EltwiseSubMod(result.a.data(), result.a.data(), temp.data(), length, modulus);
    }

    // for (size_t i = 0; i < 2 * ellnum; i++)
    for (size_t i = 0; i < ellnum; i++)
    {
        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    query.b[i].data(), length, modulus, 1);

        // note: substrct instead of addition
        intel::hexl::EltwiseSubMod(result.b.data(), result.b.data(), temp.data(), length, modulus);
    } for (size_t i = ellnum; i < 2 * ellnum; i++)
    {
        // G - query
        zeros[0] = mult_factor[i - ellnum];
        ntts.ComputeForward(temp.data(), zeros.data(), 1, 1);
        intel::hexl::EltwiseSubMod(temp.data(), temp.data(), query.b[i].data(), length, modulus);
        intel::hexl::EltwiseMultMod(temp.data(), decFirstDim[i].data(),
                    temp.data(), length, modulus, 1);

        intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), temp.data(), length, modulus);
    }

    // return to coeffs form
    ntts.ComputeInverse(result.a.data(), result.a.data(), 1, 1);
    ntts.ComputeInverse(result.b.data(), result.b.data(), 1, 1);
    result.setIsNtt(false);
#endif
}

void reorientCipher(uint64_t* cipherbuf, const std::vector<RlweCiphertext>& rotated_cipher, int32_t N1)
{
    assert(rotated_cipher.size() == N1);

    for (size_t i = 0; i < N; i++)
    {
        // cipherbuf[i * 2 * N1] = barrett_coeff(input.a[i], 0) | barrett_coeff(input.a[i], 1) << 32;
        // cipherbuf[i * 2 * N1 + 1] = barrett_coeff(input.b[i], 0) | barrett_coeff(input.b[i], 1) << 32;
        for (size_t j = 0; j < N1; j++)
        {
            uint64_t temp = rotated_cipher[j].a[i];
            cipherbuf[i * 2 * N1 + 2 * j] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
            temp = rotated_cipher[j].b[i];
            cipherbuf[i * 2 * N1 + 2 * j + 1] = barrett_coeff(temp, 0) | barrett_coeff(temp, 1) << 32;
        }
    }
}

/**
 * @brief 
 * 
 * @param result 
 * @param query1 
 * @param query2 
 * @param data store in ntt form
 */
void dimensionFoldingOld(RlweCiphertext& result, const std::vector<RlweCiphertext>& query1,
                      const std::vector<RGSWCiphertext>& query2, const std::vector<std::vector<uint64_t> >& data)
{
    int32_t firstDim_size = query1.size();
    int32_t secondDim_size = 0x01 << query2.size();
    int32_t data_size = data.size();
    assert(firstDim_size * secondDim_size == data_size);
    assert(query1[0].getIsNtt());

    uint64_t length = query1[0].getLength();
    uint64_t modulus = query1[0].getModulus();

    std::vector<RlweCiphertext> rlwes(secondDim_size, RlweCiphertext(length, modulus));
    // first dimension folding
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

    for (size_t i = 0; i < secondDim_size; i++)
    {
        for (size_t j = 0; j < firstDim_size; j++)
        {
            
            intel::hexl::EltwiseMultMod(temp1.data(), query1[j].b.data(), data[i * firstDim_size + j].data(),
                    N, modulus, 1);
            intel::hexl::EltwiseMultMod(temp2.data(), query1[j].a.data(), data[i * firstDim_size + j].data(),
                    N, modulus, 1);

            intel::hexl::EltwiseAddMod(rlwes[i].b.data(), rlwes[i].b.data(), temp1.data(),
                    N, modulus);
            intel::hexl::EltwiseAddMod(rlwes[i].a.data(), rlwes[i].a.data(), temp2.data(),
                    N, modulus);
        }
        rlwes[i].setIsNtt(true);
    }

    // second dimension folding
    for (size_t i = 0; i < query2.size(); i++)
    {
        size_t interval = 0x01 << i;
        for (size_t j = 0; j < secondDim_size; j += (2 * interval))
        {
            externalProductComp(rlwes[j], query2[i]);
            externalProduct(rlwes[j + interval], query2[i]); // TODO: inverse
            // rlwes[j] += rlwes[j + interval];
            intel::hexl::EltwiseAddMod(rlwes[j].a.data(), rlwes[j].a.data(), rlwes[j + interval].a.data(),
                                    length, modulus);
            intel::hexl::EltwiseAddMod(rlwes[j].b.data(), rlwes[j].b.data(), rlwes[j + interval].b.data(),
                                    length, modulus);
        }
    }

    result = rlwes[0];
}

void dimensionFolding(RlweCiphertext& result, const std::vector<RlweCiphertext>& query1,
                      const std::vector<RGSWCiphertext>& query2, uint64_t* data, const int32_t dataSize)
{
    int32_t firstDim_size = query1.size();
    int32_t secondDim_size = 0x01 << query2.size();
    assert(firstDim_size * secondDim_size == dataSize);
    assert(query1[0].getIsNtt());

    uint64_t length = query1[0].getLength();
    uint64_t modulus = query1[0].getModulus();

    std::vector<RlweCiphertext> rlwes(secondDim_size, RlweCiphertext(length, modulus));
    // first dimension folding
    std::vector<uint64_t> temp1(N);
    std::vector<uint64_t> temp2(N);

#ifdef SPIRAL_TIME_COLLECT
    auto start_reorient  = std::chrono::high_resolution_clock::now();
#endif
    uint64_t* cipherbuf = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * firstDim_size * N * 2);
    reorientCipher(cipherbuf, query1, firstDim_size);
#ifdef SPIRAL_TIME_COLLECT
    auto stop_reorient  = std::chrono::high_resolution_clock::now();
    auto glapsed_reorient = std::chrono::duration_cast<std::chrono::microseconds>(stop_reorient - start_reorient);
    std::cout << " reorienting ciphertexts costs " << glapsed_reorient.count() << " us." << std::endl;
#endif

#ifdef SPIRAL_TIME_COLLECT
    auto start  = std::chrono::high_resolution_clock::now();
#endif
    fastMultiplyQueryByDatabaseDim1InvCRTCipher(rlwes, data, cipherbuf, firstDim_size, secondDim_size);
#ifdef SPIRAL_TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " first dimension folding costs " << glapsed.count() << " us." << std::endl;
#endif

    for (size_t i = 0; i < secondDim_size; i++)
    {
        rlwes[i].setIsNtt(true);
    }

#ifdef SPIRAL_TIME_COLLECT
    auto start_seconddim = std::chrono::high_resolution_clock::now();
#endif
    // second dimension folding
    for (size_t i = 0; i < query2.size(); i++)
    {
        size_t interval = 0x01 << i;
        for (size_t j = 0; j < secondDim_size; j += (2 * interval))
        {
            externalProductComp(rlwes[j], query2[i]);
            externalProduct(rlwes[j + interval], query2[i]); // TODO: inverse
            // rlwes[j] += rlwes[j + interval];
            intel::hexl::EltwiseAddMod(rlwes[j].a.data(), rlwes[j].a.data(), rlwes[j + interval].a.data(),
                                    length, modulus);
            intel::hexl::EltwiseAddMod(rlwes[j].b.data(), rlwes[j].b.data(), rlwes[j + interval].b.data(),
                                    length, modulus);
        }
    }
#ifdef SPIRAL_TIME_COLLECT
    auto stop_seconddim  = std::chrono::high_resolution_clock::now();
    auto glapsed_seconddim = std::chrono::duration_cast<std::chrono::microseconds>(stop_seconddim - start_seconddim);
    std::cout << " second dimension folding costs " << glapsed_seconddim.count() << " us." << std::endl;
#endif
    result = rlwes[0];
}

void keygen(AutoKey& autokey, RGSWCiphertext& convertKey, Secret& answerKey)
{
    // autokey.keyGen(answerKey, N, true);
    std::vector<int32_t> indexLists;
    for (size_t i = 1; i < N; i <<= 1)
    {
        indexLists.push_back(N/i + 1);
    }
    autokey.keyGen(answerKey, indexLists);
    keyGen(convertKey, answerKey);
}

/**
 * @brief rns variant of keygen
 * 
 */
void keygen(AutoKeyRNS& autokey, RNSRGSWCiphertext& convertKey, Secret& answerKey)
{
    // autokey.keyGen(answerKey, N, true);
    std::vector<int32_t> indexLists;
    for (size_t i = 1; i < N; i <<= 1)
    {
        indexLists.push_back(N/i + 1);
    }
    // autokey.keyGen(answerKey, indexLists);
    autokey.keyGen(answerKey, N);
    convertKey.keyGen(answerKey);
}

void answer(RlweCiphertext& result, RlweCiphertext& query, AutoKey& autokey, RGSWCiphertext& convertKey,
            const std::vector<std::vector<uint64_t> >& data, const spiral_parameters& para)
{
    uint64_t length = autokey.getLength();
    uint64_t modulus = autokey.getModulus();
    int32_t ellgsw = para.ellgsw;

    std::vector<RlweCiphertext> expansion_rlwes(para.rlwesNum, RlweCiphertext(length, modulus));
    std::vector<RlweCiphertext> expansion_rgsws(para.rgswsNum * ellgsw, RlweCiphertext(length, modulus));

int ntimes = 1;
#ifdef SPIRAL_TIME_COLLECT
    auto start = std::chrono::high_resolution_clock::now();
#endif
for (size_t i = 0; i < ntimes; i++)
{
    queryExpand(expansion_rlwes, expansion_rgsws, query, autokey, para.rlwesNum, para.rgswsNum, ellgsw);
}
#ifdef SPIRAL_TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " query expansion costs " << glapsed.count() << " us." << std::endl;
#endif

    // evaluate external product
    // std::vector<RGSWCiphertext> converted_rgsw(rgswsNum, RGSWCiphertext(N, modulus, 8, 0, 0x01 << 7));
    std::vector<RGSWCiphertext> converted_rgsw(para.rgswsNum, RGSWCiphertext(N, modulus, ellgsw, 0, 0x01 << (int32_t)ceil(56./ ellgsw)));

    for (size_t i = 0; i < para.rgswsNum; i++)
        regevToGSW(converted_rgsw[i], expansion_rgsws.begin() + ellgsw * i, expansion_rgsws.begin() + ellgsw * (i+1), convertKey, ellgsw);


    // first and subsequent dimensions folding
    if(!expansion_rlwes[0].getIsNtt())
    {
#ifdef INTEL_HEXL   
        intel::hexl::NTT ntts = autokey.getNTT();
        for (size_t i = 0; i < para.rlwesNum; i++)
        {
            ntts.ComputeForward(expansion_rlwes[i].a.data(), expansion_rlwes[i].a.data(), 1, 1);
            ntts.ComputeForward(expansion_rlwes[i].b.data(), expansion_rlwes[i].b.data(), 1, 1);
            expansion_rlwes[i].setIsNtt(true);
        }
#endif
    }
    dimensionFoldingOld(result, expansion_rlwes, converted_rgsw, data);
}

void answer(RlweCiphertext& result, RlweCiphertext& query, AutoKey& autokey, RGSWCiphertext& convertKey,
            uint64_t* data, const spiral_parameters& para)
{
    uint64_t length = autokey.getLength();
    uint64_t modulus = autokey.getModulus();
    int32_t ellgsw = para.ellgsw;

    std::vector<RlweCiphertext> expansion_rlwes(para.rlwesNum, RlweCiphertext(length, modulus));
    std::vector<RlweCiphertext> expansion_rgsws(para.rgswsNum * ellgsw, RlweCiphertext(length, modulus));

int ntimes = 1;
#ifdef SPIRAL_TIME_COLLECT
    auto start = std::chrono::high_resolution_clock::now();
#endif
for (size_t i = 0; i < ntimes; i++)
{
    queryExpand(expansion_rlwes, expansion_rgsws, query, autokey, para.rlwesNum, para.rgswsNum, ellgsw);
}
#ifdef SPIRAL_TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " query expansion costs " << glapsed.count() << " us." << std::endl;
#endif

    // evaluate external product
    // std::vector<RGSWCiphertext> converted_rgsw(rgswsNum, RGSWCiphertext(N, modulus, 8, 0, 0x01 << 7));
    std::vector<RGSWCiphertext> converted_rgsw(para.rgswsNum, RGSWCiphertext(N, modulus, ellgsw, 0, 0x01 << (int32_t)ceil(56./ ellgsw)));

    for (size_t i = 0; i < para.rgswsNum; i++)
        regevToGSW(converted_rgsw[i], expansion_rgsws.begin() + ellgsw * i, expansion_rgsws.begin() + ellgsw * (i+1), convertKey, ellgsw);


    // first and subsequent dimensions folding
    if(!expansion_rlwes[0].getIsNtt())
    {
#ifdef INTEL_HEXL   
        intel::hexl::NTT ntts = autokey.getNTT();
        for (size_t i = 0; i < para.rlwesNum; i++)
        {
            ntts.ComputeForward(expansion_rlwes[i].a.data(), expansion_rlwes[i].a.data(), 1, 1);
            ntts.ComputeForward(expansion_rlwes[i].b.data(), expansion_rlwes[i].b.data(), 1, 1);
            expansion_rlwes[i].setIsNtt(true);
        }
#endif
    }
    
    int32_t dataSize = para.rlwesNum * (0x01 << para.rgswsNum);
    dimensionFolding(result, expansion_rlwes, converted_rgsw, data, dataSize);
}

void answer(RlweCiphertext& result, std::vector<RlweCiphertext>& query, AutoKeyRNS& autokey, RNSRGSWCiphertext& convertKey,
            uint64_t* data, const spiral_parameters& para)
{
    uint64_t length = autokey.getLength();
    uint64_t modulus = autokey.getModulus1();
    int32_t ellgsw = para.ellgsw;

    std::vector<RlweCiphertext> expansion_rlwes(para.rlwesNum, RlweCiphertext(length, modulus));
    std::vector<RlweCiphertext> expansion_rgsws(para.rgswsNum * ellgsw, RlweCiphertext(length, modulus));
    for (size_t i = 0; i < para.rgswsNum * ellgsw; i++) expansion_rgsws.push_back(RlweCiphertext(length, convertKey.getSecModulus()));

#ifdef SPIRAL_TIME_COLLECT
    auto start_expand  = std::chrono::high_resolution_clock::now();
#endif

    queryExpand(expansion_rlwes, expansion_rgsws, query, autokey, para.rlwesNum, para.rgswsNum, ellgsw);

#ifdef SPIRAL_TIME_COLLECT
    auto stop_expand  = std::chrono::high_resolution_clock::now();

    auto glapsed_expand = std::chrono::duration_cast<std::chrono::microseconds>(stop_expand - start_expand);
    std::cout << " query expansion costs " << glapsed_expand.count() << " us." << std::endl;
#endif

#ifdef SPIRAL_TIME_COLLECT
    auto start  = std::chrono::high_resolution_clock::now();
#endif

    // evaluate external product;
    std::vector<RGSWCiphertext> converted_rgsw(para.rgswsNum, RGSWCiphertext(N, modulus, ellgsw, 0, 0x01 << (int32_t)ceil(56./ ellgsw)));

    int rgsws = para.rgswsNum * ellgsw;
    for (size_t i = 0; i < para.rgswsNum; i++)
    {
        regevToGSW(converted_rgsw[i], expansion_rgsws.begin() + ellgsw * i, expansion_rgsws.begin() + ellgsw * (i+1),
                   expansion_rgsws.begin() + rgsws + ellgsw * i, expansion_rgsws.begin() + rgsws + ellgsw * (i+1),   
                   convertKey, ellgsw);
    }

#ifdef SPIRAL_TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " convert to gsw costs " << glapsed.count() << " us." << std::endl;
#endif

    // first and subsequent dimensions folding
    if(!expansion_rlwes[0].getIsNtt())
    {
#ifdef INTEL_HEXL   
        intel::hexl::NTT ntts = autokey.getNTT1();
        for (size_t i = 0; i < para.rlwesNum; i++)
        {
            ntts.ComputeForward(expansion_rlwes[i].a.data(), expansion_rlwes[i].a.data(), 1, 1);
            ntts.ComputeForward(expansion_rlwes[i].b.data(), expansion_rlwes[i].b.data(), 1, 1);
            expansion_rlwes[i].setIsNtt(true);
        }
#endif
    }
    
    int32_t dataSize = para.rlwesNum * (0x01 << para.rgswsNum);
    dimensionFolding(result, expansion_rlwes, converted_rgsw, data, dataSize);
}

void answer(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& query, AutoKeyRNS& autokey, RNSRGSWCiphertext& convertKey,
            uint64_t* data, const spiral_parameters& para)
{
    uint64_t length = autokey.getLength();
    uint64_t modulus = autokey.getModulus1();
    int32_t ellgsw = para.ellgsw;

    std::vector<RlweCiphertext> expansion_rlwes(para.rlwesNum, RlweCiphertext(length, modulus));
    std::vector<RlweCiphertext> expansion_rgsws(para.rgswsNum * ellgsw, RlweCiphertext(length, modulus));
    for (size_t i = 0; i < para.rgswsNum * ellgsw; i++) expansion_rgsws.push_back(RlweCiphertext(length, convertKey.getSecModulus()));

#ifdef SPIRAL_TIME_COLLECT
    auto start_expand  = std::chrono::high_resolution_clock::now();
#endif

    queryExpand(expansion_rlwes, expansion_rgsws, query, autokey, para.rlwesNum, para.rgswsNum, ellgsw);

#ifdef SPIRAL_TIME_COLLECT
    auto stop_expand  = std::chrono::high_resolution_clock::now();

    auto glapsed_expand = std::chrono::duration_cast<std::chrono::microseconds>(stop_expand - start_expand);
    std::cout << " query expansion costs " << glapsed_expand.count() << " us." << std::endl;
#endif

#ifdef SPIRAL_TIME_COLLECT
    auto start  = std::chrono::high_resolution_clock::now();
#endif

    // evaluate external product;
    std::vector<RGSWCiphertext> converted_rgsw(para.rgswsNum, RGSWCiphertext(N, modulus, ellgsw, 0, 0x01 << (int32_t)ceil(56./ ellgsw)));

    int rgsws = para.rgswsNum * ellgsw;
    for (size_t i = 0; i < para.rgswsNum; i++)
    {
        regevToGSW(converted_rgsw[i], expansion_rgsws.begin() + ellgsw * i, expansion_rgsws.begin() + ellgsw * (i+1),
                   expansion_rgsws.begin() + rgsws + ellgsw * i, expansion_rgsws.begin() + rgsws + ellgsw * (i+1),   
                   convertKey, ellgsw);
    }

#ifdef SPIRAL_TIME_COLLECT
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " convert to gsw costs " << glapsed.count() << " us." << std::endl;
#endif

    // first and subsequent dimensions folding
    if(!expansion_rlwes[0].getIsNtt())
    {
#ifdef INTEL_HEXL   
        intel::hexl::NTT ntts = autokey.getNTT1();
        for (size_t i = 0; i < para.rlwesNum; i++)
        {
            ntts.ComputeForward(expansion_rlwes[i].a.data(), expansion_rlwes[i].a.data(), 1, 1);
            ntts.ComputeForward(expansion_rlwes[i].b.data(), expansion_rlwes[i].b.data(), 1, 1);
            expansion_rlwes[i].setIsNtt(true);
        }
#endif
    }
    
    int32_t dataSize = para.rlwesNum * (0x01 << para.rgswsNum);
    for (size_t i = 0; i < pow(para.n, 2); i++)
        dimensionFolding(result[i], expansion_rlwes, converted_rgsw, data, dataSize);
}

void answer(RlweCiphertext& result, const std::vector<RlweCiphertext>& expansion_rlwes,
            const std::vector<RGSWCiphertext>& expansion_rgsws, const std::vector<std::vector<uint64_t> >& data)
{
    dimensionFoldingOld(result, expansion_rlwes, expansion_rgsws, data);
}

void answer(RlweCiphertext& result, const std::vector<RlweCiphertext>& expansion_rlwes,
            const std::vector<RGSWCiphertext>& expansion_rgsws, uint64_t* data)
{
    int32_t dataSize = expansion_rlwes.size() * (0x01 << expansion_rgsws.size());
    dimensionFolding(result, expansion_rlwes, expansion_rgsws, data, dataSize);
}

void recover(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& answerKey)
{
    decrypt(message, cipher, answerKey);

    uint64_t modulus = cipher.getModulus();
    uint64_t halfMolulus = modulus >> 1;
    uint64_t delta = modulus / SPIRAL_P;

    int64_t temp;
    for (auto iter = message.begin(); iter != message.end(); iter++)
    {
        // TODO: check it
        if (*iter > halfMolulus)
            temp = *iter - modulus;
        else
            temp = *iter;
        temp = (int64_t)roundl((long double)temp / delta);

        if (temp < 0)
            temp += SPIRAL_P;
        *iter = temp;
    }
}

void recover_noise(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& answerKey, uint64_t spiral_p)
{
    decrypt(message, cipher, answerKey);

    int32_t length = cipher.getLength();
    uint64_t modulus = cipher.getModulus();
    uint64_t halfMolulus = modulus >> 1;
    uint64_t delta = modulus / spiral_p;

    std::vector<uint64_t> err(length, 0);
    std::vector<uint64_t> deltaM(length, 0);
    copy(message.begin(), message.end(), err.begin());

    int64_t temp;
    for (auto iter = message.begin(); iter != message.end(); iter++)
    {
        // TODO: check it
        if (*iter > halfMolulus)
            temp = *iter - modulus;
        else
            temp = *iter;
        temp = (int64_t)roundl((long double)temp / delta);

        if (temp < 0)
            temp += spiral_p;
        *iter = temp;
    }

    // print errors (noise)
#ifdef INTEL_HEXL
    intel::hexl::EltwiseFMAMod(deltaM.data(), message.data(), delta,
                nullptr, length, modulus, 1);
    intel::hexl::EltwiseSubMod(err.data(), err.data(), deltaM.data(),
                    length, modulus);
    showLargeVector(err, "The err is ");
#endif
}

/**
 * @brief print uint64_t
 * 
 * @param database 
 * @param target_row 
 * @param target_col 
 * @param rows 
 */
void print_target_record(std::vector<std::vector<uint64_t> >& database,
                         int32_t target_row, int32_t target_col,
                         int32_t rows, uint64_t spiral_p)
{
    std::vector<uint64_t> message(N, 0);
    int64_t index = target_col * rows + target_row;

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
    ntts.ComputeInverse(message.data(), database[index].data(), 1, 1);
#endif

    for (size_t j = 0; j < N; j++)
    {
        uint64_t db = message[j];
        message[j] = db > crtMod/2 ? db + spiral_p - crtMod : db;
    }

    showLargeVector(message, "The correct   value is ");
}

bool check_target_record(std::vector<std::vector<uint64_t> >& database,
                         int32_t target_row, int32_t target_col,
                         int32_t rows, uint64_t spiral_p, std::vector<uint64_t>& decryptd_message)
{
    std::vector<uint64_t> message(N, 0);
    int64_t index = target_col * rows + target_row;

#ifdef INTEL_HEXL
    intel::hexl::NTT ntts(N, crtMod, root_of_unity_crt);
    ntts.ComputeInverse(message.data(), database[index].data(), 1, 1);
#endif

    bool istrue = true;
    for (size_t j = 0; j < N; j++)
    {
        uint64_t db = message[j];
        message[j] = db > crtMod/2 ? db + spiral_p - crtMod : db;
        if (decryptd_message[j] != message[j]) istrue = false;
    }

    showLargeVector(message, "The correct   value is ");

    return istrue;
}
