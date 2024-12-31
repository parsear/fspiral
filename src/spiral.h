#ifndef SPIRAL_H
#define SPIRAL_H

#include <stddef.h>
#include <vector>

#include "lwe.h"
#include "secret.h"

#define SPIRAL_P 256 // 8 bits

struct spiral_parameters {
    int32_t rlwesNum = 0x01 << 9; // the number of RLWEs
    int32_t rgswsNum = 6; // the number of RGSWs
    // plaintext modulus
    uint64_t spiral_p = 0x01 << 8;
    // expansion for rlwe ciphertexts
    int32_t ell_rlwes_expansion = 14;
    uint64_t base_rlwes_expansion = 0;
    uint64_t bg_rlwes_expansion = 0x01 << 4;
    // expansion for rgsw ciphertexts
    int32_t ell_rgsws_expansion = 14;
    uint64_t base_rgsws_expansion = 0;
    uint64_t bg_rgsws_expansion = 0x01 << 4;
    // convert rlwe to rgsw
    int32_t ell_convert = 7;
    uint64_t base_convert = 0x01 << 7;
    uint64_t bg_convert = 0x01 << 7;
    
    // rgsw parameters
    int32_t ellgsw = 8;
    uint64_t basegsw = 0;
    uint64_t bggsw = 0x01 << 7;

    int32_t n = 1;
};

class RNSRGSWCiphertext
{
private:
    /* data */
    int32_t length = N;
    uint64_t modulus = crtMod;
    uint64_t modulus2 = bsMod;

    int32_t ellnum = 4;
    uint64_t PP = 0; // 0x01 << 20; // 18, the base
    uint64_t BBg = 0x01 << 20;

    bool isntt = false;

// for 2nd folding
#ifdef INTEL_HEXL
    intel::hexl::NTT ntts;
    intel::hexl::NTT ntts2;
#endif

public:
    std::vector<std::vector<uint64_t> > a;
    std::vector<std::vector<uint64_t> > b;
    std::vector<std::vector<uint64_t> > arns;
    std::vector<std::vector<uint64_t> > brns;

    RNSRGSWCiphertext(/* args */);

    RNSRGSWCiphertext(int32_t len, uint64_t module);

    RNSRGSWCiphertext(int32_t len, uint64_t module, int32_t elln, uint64_t base, uint64_t bg); 

    void keyGen(Secret& secret);

    uint64_t getModulus() const;

    uint64_t getSecModulus() const;

    int32_t getLength() const;

    int32_t getEllnum() const;

    uint64_t getBg() const;

    uint64_t getBase() const;

    bool getIsNtt() const;

    void setIsNtt(bool is);

#ifdef INTEL_HEXL
    intel::hexl::NTT getNTT() const;
    intel::hexl::NTT getNTT2() const;
#endif

    // ~RGSWCiphertext();
};

// void get_spiral_parameters_one(spiral_parameters& para);
void get_spiral_parameters_one(spiral_parameters& para, int32_t rlwesNum = 0x01 << 9,
                               int32_t rgswsNum = 6, uint64_t spiral_p = 0x01 << 8);

void get_spiralstream_parameters_one(spiral_parameters& para, int32_t rlwesNum = 0x01 << 9,
                                     int32_t rgswsNum = 6, uint64_t spiral_p = 0x01 << 21);

void get_spiral_rns_parameters_one(spiral_parameters& para, int32_t rlwesNum = 0x01 << 9,
                               int32_t rgswsNum = 6, uint64_t spiral_p = 0x01 << 8, int32_t n = 1);

void sample_spiral_database(std::vector<std::vector<uint64_t> >& database, uint64_t spiral_p = SPIRAL_P);

void spiral_sdatabase_tocrt(uint64_t* datacrt, std::vector<std::vector<uint64_t> >& data_ntt, int32_t N1, int32_t n = 1);

void queryGenerate(RlweCiphertext& cipher, Secret& queryKey, const int32_t row, const int32_t rlwesNum);

void queryGenerate(RlweCiphertext& cipher, Secret& queryKey,
                   const int32_t row, const int32_t col,
                   const spiral_parameters& para);

/**
 * @brief rns variant of queryGenerate
 * 
 */
void queryGenerate(std::vector<RlweCiphertext>& cipher, Secret& queryKey,
                   const int32_t row, const int32_t col,
                   const spiral_parameters& para);

/**
 * @brief SpiralStream variant of queryGenerate
 * 
 * @param query_rlwe 
 * @param query_rgsw 
 * @param queryKey 
 * @param row 
 * @param col 
 * @param rlwesNum 
 * @param rgswsNum 
 */
void queryGenerate(std::vector<RlweCiphertext>& query_rlwe, std::vector<RGSWCiphertext>& query_rgsw,
                   Secret& queryKey,
                   const int32_t row, const int32_t col,
                   const spiral_parameters& para);

void queryExpand(std::vector<RlweCiphertext>& result, RlweCiphertext& input,
                 const AutoKey& autokey, const int32_t rlwesNum);

void queryExpand(std::vector<RlweCiphertext>& expansion_rlwes, std::vector<RlweCiphertext>& expansion_rgsws,
                 RlweCiphertext& input, const AutoKey& autokey,
                 const int32_t rlwesNum, const int32_t rgswsNum, const int32_t ellgsw);

void queryExpand(std::vector<RlweCiphertext>& expansion_rlwes, std::vector<RlweCiphertext>& expansion_rgsws,
                 std::vector<RlweCiphertext>& input, const AutoKeyRNS& autokey,
                 const int32_t rlwesNum, const int32_t rgswsNum, const int32_t ellgsw);

void keyGen(RGSWCiphertext& result, Secret& secret);

void regevToGSW(RGSWCiphertext& result, std::vector<RlweCiphertext>& input,
                const RGSWCiphertext& convertKey, const int32_t ellgsw, const int32_t start_index = 0);

void regevToGSW(RGSWCiphertext& result, std::vector<RlweCiphertext>::iterator input,
                std::vector<RlweCiphertext>::iterator input_end,
                const RGSWCiphertext& convertKey, const int32_t ellgsw);

void regevToGSW(RGSWCiphertext& result, std::vector<RlweCiphertext>::iterator input,
                std::vector<RlweCiphertext>::iterator input_end,
                std::vector<RlweCiphertext>::iterator input_rns,
                std::vector<RlweCiphertext>::iterator input_rns_end,
                const RNSRGSWCiphertext& convertKey, const int32_t ellgsw);

void keygen(AutoKey& autokey, RGSWCiphertext& convertKey, Secret& answerKey);

void keygen(AutoKeyRNS& autokey, RNSRGSWCiphertext& convertKey, Secret& answerKey);

// void keygen(AutoKeyRNS& autokey, std::vector<RGSWCiphertext>& convertKey, Secret& answerKey);

void answer(RlweCiphertext& result, RlweCiphertext& query, AutoKey& autokey, RGSWCiphertext& convertKey,
            const std::vector<std::vector<uint64_t> >& data, const spiral_parameters& para);

void answer(RlweCiphertext& result, RlweCiphertext& query, AutoKey& autokey, RGSWCiphertext& convertKey,
            uint64_t* data, const spiral_parameters& para);

void answer(RlweCiphertext& result, std::vector<RlweCiphertext>& query, AutoKeyRNS& autokey, RNSRGSWCiphertext& convertKey,
            uint64_t* data, const spiral_parameters& para);

void answer(std::vector<RlweCiphertext>& result, std::vector<RlweCiphertext>& query, AutoKeyRNS& autokey, RNSRGSWCiphertext& convertKey,
            uint64_t* data, const spiral_parameters& para);

/**
void answer(RlweCiphertext& result, const std::vector<RlweCiphertext>& expansion_rlwes,
            const std::vector<RGSWCiphertext>& expansion_rgsws, const std::vector<std::vector<uint64_t> >& data);
**/

void answer(RlweCiphertext& result, const std::vector<RlweCiphertext>& expansion_rlwes,
            const std::vector<RGSWCiphertext>& expansion_rgsws, uint64_t* data);

/**
 * @brief SpiralStream variant of answer
 * 
 * @param result 
 * @param expansion_rlwes 
 * @param expansion_rgsws 
 * @param data 
 */
void answer(RlweCiphertext& result, const std::vector<RlweCiphertext>& expansion_rlwes,
            const std::vector<RGSWCiphertext>& expansion_rgsws, const std::vector<std::vector<uint64_t> >& data);

void recover(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& answerKey);

void recover_noise(std::vector<uint64_t>& message, RlweCiphertext& cipher, Secret& answerKey, uint64_t spiral_p = SPIRAL_P);

void print_target_record(std::vector<std::vector<uint64_t> >& data,
                         int32_t target_row, int32_t target_col,
                         int32_t rows, uint64_t spiral_p = SPIRAL_P);

bool check_target_record(std::vector<std::vector<uint64_t> >& database,
                         int32_t target_row, int32_t target_col,
                         int32_t rows, uint64_t spiral_p, std::vector<uint64_t>& decryptd_message);

#endif
