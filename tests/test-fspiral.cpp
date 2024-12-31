#include <stddef.h>
#include <vector>
#include <iostream>
#include <time.h>
#include <chrono>

#include "pir.h"
#include "spiral.h"

#ifdef INTEL_HEXL
#include <hexl/hexl.hpp>
#endif

void print_config(const spiral_parameters& para, int32_t target_row, int32_t target_col)
{
    long double database_size = (long double)para.rlwesNum * (0x01 << para.rgswsNum) * N / pow(2, 23) * log2(para.spiral_p) * para.n * para.n;
    std::cout << "Database: " << para.rlwesNum << " * " << (0x01 << para.rgswsNum) << " * R_p";
    std::cout << ", (total " << database_size << " MB)" << std::endl;
    std::cout << std::endl;

    std::cout << "target (row, col) = (" << target_row << ", " << target_col << ")" << std::endl;
}

void test_query_expansion()
{    
    int32_t l = 4; // 4; // 6; // 5
    int32_t rlwesNum = 0x01 << l; // the number of RLWEs

    int32_t target_row = rand() % rlwesNum;

    std::cout << "rlweNum: " << rlwesNum << std::endl;
    std::cout << "target row: " << target_row << std::endl;

    uint64_t modulus = crtMod;
    Secret answerKey(modulus, false);
    uint64_t length = answerKey.getLength();

    // query
    RlweCiphertext input(N, modulus);
    queryGenerate(input, answerKey, target_row, rlwesNum);

    std::vector<RlweCiphertext> expansion_rlwes(rlwesNum, RlweCiphertext(N, modulus));
    std::vector<std::vector<uint64_t> > decryptd_message(rlwesNum, std::vector<uint64_t>(N));

    auto start_pk = std::chrono::high_resolution_clock::now();
    AutoKey autokey(length, modulus, 14, 0, 0x01 << 4);

    // autokey.keyGen(answerKey, rlwesNum, true);
    autokey.keyGen(answerKey, N, true);
    auto stop_pk = std::chrono::high_resolution_clock::now();

    auto glapsed_pk = std::chrono::duration_cast<std::chrono::microseconds>(stop_pk - start_pk);
    std::cout << "building packing keys costs " << glapsed_pk.count() << " us." << std::endl;

    // the results
int ntimes = 1;
    auto start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; i++)
{
    queryExpand(expansion_rlwes, input, autokey, rlwesNum);
}
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " query expansion costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    decrypt(decryptd_message[target_row], expansion_rlwes[target_row], answerKey);

    std::cout << std::endl;
    showLargeVector(decryptd_message[target_row], "The decrypted value is ");
}

void check_expansion_rgsws(std::vector<RlweCiphertext>& expansion_rgsws, Secret& answerKey,
                           const int32_t rgswsNum, const spiral_parameters& para, const int32_t ellgsw = 8)
{
    uint64_t length = answerKey.getLength();
    uint64_t modulus = answerKey.getModulus();
    std::vector<uint64_t> message(length, 0), decryptd_message(length, 0);
    
    // evaluate external product
    RGSWCiphertext convertKey(N, modulus, 7, 0x01 << 7, 0x01 << 7);
    keyGen(convertKey, answerKey);
    std::cout << "keyGen Done." << std::endl;

    // std::vector<RGSWCiphertext> result(rgswsNum, RGSWCiphertext(N, modulus, 8, 0, 0x01 << 7));
    std::vector<RGSWCiphertext> result(para.rgswsNum, RGSWCiphertext(N, modulus, ellgsw, 0, 0x01 << (int32_t)ceil(56./ ellgsw)));

    for (size_t i = 0; i < rgswsNum; i++)
        regevToGSW(result[i], expansion_rgsws.begin() + ellgsw * i, expansion_rgsws.begin() + ellgsw * (i+1), convertKey, ellgsw);
    std::cout << "regevToGSW Done." << std::endl;

    // test if the RGSW ciphertext is valid.
    sample_random(message, modulus);
    showLargeVector(message, "The message is ");

    RlweCiphertext query(N, modulus), answer(N, modulus);
    encrypt(query, answerKey, message);
    externalProduct(answer, query, result[0]);

    // decrypt message
    decrypt(decryptd_message, answer, answerKey);

    showLargeVector(decryptd_message, "The decrypted value is ");
}

void check_expansion_rgsws_rns(std::vector<RlweCiphertext>& expansion_rgsws, Secret& answerKey,
                           const int32_t rgswsNum, const spiral_parameters& para, const int32_t ellgsw = 8)
{
    uint64_t length = answerKey.getLength();
    uint64_t modulus = answerKey.getModulus();
    std::vector<uint64_t> message(length, 0), decryptd_message(length, 0);
    
    // evaluate external product
    RNSRGSWCiphertext convertKey(N, modulus, 8, 0, 0x01 << 10);
    // keyGen(convertKey, answerKey);
    convertKey.keyGen(answerKey);
    std::cout << "keyGen Done." << std::endl;

    std::vector<RGSWCiphertext> result(para.rgswsNum, RGSWCiphertext(N, modulus, ellgsw, 0, 0x01 << (int32_t)ceil(56./ ellgsw)));

    int rgsws = rgswsNum * ellgsw;
    for (size_t i = 0; i < rgswsNum; i++)
    {
        regevToGSW(result[i], expansion_rgsws.begin() + ellgsw * i, expansion_rgsws.begin() + ellgsw * (i+1),
                   expansion_rgsws.begin() + rgsws + ellgsw * i, expansion_rgsws.begin() + rgsws + ellgsw * (i+1),
                   convertKey, ellgsw);
    }
    std::cout << "regevToGSW Done." << std::endl;

    // test if the RGSW ciphertext is valid.
    sample_random(message, modulus);
    showLargeVector(message, "The message is ");

    RlweCiphertext query(N, modulus), answer(N, modulus);
    encrypt(query, answerKey, message);
    externalProduct(answer, query, result[0]);

    // decrypt message
    decrypt(decryptd_message, answer, answerKey);

    showLargeVector(decryptd_message, "The decrypted value is ");
}

void test_query_expansion_split()
{    
    // total = 2^6 * 2^4
    int32_t rlwesNum = 0x01 << 6; // the number of RLWEs
    int32_t rgswsNum = 4; // the number of RGSWs. note regevToGSW is required

    int32_t target_row = rand() % rlwesNum;
    int32_t target_col = rand() % (0x01 << rgswsNum);

    spiral_parameters para;
    get_spiral_parameters_one(para, rlwesNum, rgswsNum);
    print_config(para, target_row, target_col);

    uint64_t modulus = crtMod;
    Secret answerKey(modulus, false);
    uint64_t length = answerKey.getLength();

    // query
    RlweCiphertext input(N, modulus);
    queryGenerate(input, answerKey, target_row, target_col, para);

    std::vector<RlweCiphertext> expansion_rlwes(rlwesNum, RlweCiphertext(N, modulus));
    std::vector<RlweCiphertext> expansion_rgsws(rgswsNum * para.ellgsw, RlweCiphertext(N, modulus));

    std::vector<std::vector<uint64_t> > decryptd_message(rlwesNum, std::vector<uint64_t>(N));

    auto start_pk = std::chrono::high_resolution_clock::now();
    AutoKey autokey(length, modulus, 14, 0, 0x01 << 4);

    std::vector<int32_t> indexLists;
    for (size_t i = 1; i < N; i <<= 1) indexLists.push_back(N/i + 1);
    autokey.keyGen(answerKey, indexLists);

    auto stop_pk = std::chrono::high_resolution_clock::now();

    auto glapsed_pk = std::chrono::duration_cast<std::chrono::microseconds>(stop_pk - start_pk);
    std::cout << "building packing keys costs " << glapsed_pk.count() << " us." << std::endl;

    // the results
int ntimes = 1;
    auto start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; i++)
{
    queryExpand(expansion_rlwes, expansion_rgsws, input, autokey, rlwesNum, rgswsNum, para.ellgsw);
}
    auto stop  = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << ntimes << " query expansion costs " << glapsed.count() << " us." << std::endl;

    // decrypt message
    decrypt(decryptd_message[target_row], expansion_rlwes[target_row], answerKey);

    std::cout << std::endl;
    showLargeVector(decryptd_message[target_row], "The decrypted rlwe value is ");

    int32_t test_index = rand() % (rgswsNum * para.ellgsw);
    decrypt(decryptd_message[target_row], expansion_rgsws[test_index], answerKey);
    std::cout << std::endl;
    showLargeVector(decryptd_message[target_row], "The decrypted rgsw value is ");
    
    if (expansion_rgsws[test_index].getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = answerKey.getNTT();
        ntts.ComputeInverse(expansion_rgsws[test_index].a.data(), expansion_rgsws[test_index].a.data(), 1, 1);
        ntts.ComputeInverse(expansion_rgsws[test_index].b.data(), expansion_rgsws[test_index].b.data(), 1, 1);
        expansion_rgsws[test_index].setIsNtt(false);
#endif
    }

    // test expansion rgsws
    check_expansion_rgsws(expansion_rgsws, answerKey, rgswsNum, para, para.ellgsw);
    std::cout << "lowest bit of col: " << (target_col & 0x01) << std::endl;
}

void test_query_expansion_split_rns()
{    
    // total = 2^6 * 2^4
    int32_t rlwesNum = 0x01 << 6; // the number of RLWEs
    int32_t rgswsNum = 4; // the number of RGSWs. note regevToGSW is required

    int32_t target_row = rand() % rlwesNum;
    int32_t target_col = rand() % (0x01 << rgswsNum);

    spiral_parameters para;
    get_spiral_rns_parameters_one(para, rlwesNum, rgswsNum);
    print_config(para, target_row, target_col);

    uint64_t modulus = crtMod;
    uint64_t modulus2 = bsMod;
    Secret answerKey(modulus, false);
    uint64_t length = answerKey.getLength();

    // query
    std::vector<RlweCiphertext> input(1, RlweCiphertext(N, modulus));
    input.push_back(RlweCiphertext(N, modulus2));
    queryGenerate(input, answerKey, target_row, target_col, para);

    std::vector<RlweCiphertext> expansion_rlwes(rlwesNum, RlweCiphertext(N, modulus));
    std::vector<RlweCiphertext> expansion_rgsws(2 * rgswsNum * para.ellgsw, RlweCiphertext(N, modulus));

    std::vector<std::vector<uint64_t> > decryptd_message(rlwesNum, std::vector<uint64_t>(N));

    AutoKeyRNS autokey(length, modulus, modulus2, para.ell_rlwes_expansion, para.base_rlwes_expansion, para.bg_rlwes_expansion);
    autokey.keyGen(answerKey, N);
    std::cout << "building packing keys done. " << std::endl;

    // the results
    queryExpand(expansion_rlwes, expansion_rgsws, input, autokey, rlwesNum, rgswsNum, para.ellgsw);
    std::cout << " query expansion Done. " << std::endl;

    // decrypt messages
    decrypt(decryptd_message[target_row], expansion_rlwes[target_row], answerKey);

    std::cout << std::endl;
    showLargeVector(decryptd_message[target_row], "The decrypted rlwe value is ");

    int32_t test_index = rand() % (rgswsNum * para.ellgsw);
    decrypt(decryptd_message[target_row], expansion_rgsws[test_index], answerKey);
    std::cout << std::endl;
    showLargeVector(decryptd_message[target_row], "The decrypted rgsw value is ");
    
    if (expansion_rgsws[test_index].getIsNtt())
    {
#ifdef INTEL_HEXL
        intel::hexl::NTT ntts = answerKey.getNTT();
        ntts.ComputeInverse(expansion_rgsws[test_index].a.data(), expansion_rgsws[test_index].a.data(), 1, 1);
        ntts.ComputeInverse(expansion_rgsws[test_index].b.data(), expansion_rgsws[test_index].b.data(), 1, 1);
        expansion_rgsws[test_index].setIsNtt(false);
#endif
    }

    // test expansion rgsws
    check_expansion_rgsws_rns(expansion_rgsws, answerKey, rgswsNum, para, para.ellgsw);
    std::cout << "lowest bit of col: " << (target_col & 0x01) << std::endl;
}

void test_regev_product()
{
    int32_t ellgsw = 8; // the ell for output RGSW ciphertext
    uint64_t output_bg = 0x01 << 7;

    uint64_t modulus = crtMod;
    
    Secret answerKey(modulus); // default rlwe secret, store in ntt form
    std::vector<RlweCiphertext> cipher(ellgsw, RlweCiphertext(N, modulus));

    std::vector<uint64_t> message(N), decryptd_message(N);
    // encrypt message
    for (size_t i = 0; i < ellgsw; i++)
    {
        message[0] = pow(output_bg, i); // here we assume without approximate decompisition
        for (size_t j = 1; j < N; j++) message[j] = 0; 
        encrypt(cipher[i], answerKey, message);
    }

    // evaluate external product
    RGSWCiphertext convertKey(N, modulus, 7, 0x01 << 7, 0x01 << 7);
    auto start_qu = std::chrono::high_resolution_clock::now(); 
    // convertKey.keyGen(answerKey, 3, true);
    keyGen(convertKey, answerKey);
    auto stop_qu = std::chrono::high_resolution_clock::now();

    auto glapsed_qu = std::chrono::duration_cast<std::chrono::microseconds>(stop_qu - start_qu);
    std::cout << " query costs " << glapsed_qu.count() << " us." << std::endl;

    RGSWCiphertext result(N, modulus, 8, 0, 0x01 << 7);
int32_t ntimes = 1; // 16 * 32;

    auto start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < ntimes; i++)
{
    regevToGSW(result, cipher, convertKey, ellgsw);
}
    auto stop = std::chrono::high_resolution_clock::now();

    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "each external product costs " << glapsed.count() << " us." << std::endl;


    // test if the RGSW ciphertext is valid.
    sample_random(message, modulus);
    showLargeVector(message, "The message is ");

    RlweCiphertext query(N, modulus), answer(N, modulus);
    encrypt(query, answerKey, message);
    externalProduct(answer, query, result);

/**/
    // decrypt message
    auto start_de = std::chrono::high_resolution_clock::now();
    decrypt(decryptd_message, answer, answerKey);
    auto stop_de = std::chrono::high_resolution_clock::now();

    auto glapsed_de = std::chrono::duration_cast<std::chrono::microseconds>(stop_de - start_de);
    std::cout << " decrypt costs " << glapsed_de.count() << " us." << std::endl;
/**/
    showLargeVector(decryptd_message, "The decrypted value is ");
}

void test_spiral_old()
{
    // total = 2^9 * 2^6
    int32_t rlwesNum = 0x01 << 9; // the number of RLWEs
    int32_t rgswsNum = 6; // the number of RGSWs. note regevToGSW is required
    uint64_t spiral_p = 0x01 << 11; // 0x01 << 16;

    int32_t target_row = rand() % rlwesNum;
    int32_t target_col = rand() % (0x01 << rgswsNum);

    spiral_parameters para;
    get_spiral_parameters_one(para, rlwesNum, rgswsNum, spiral_p);
    print_config(para, target_row, target_col);

    uint64_t modulus = crtMod;
    Secret answerKey(modulus, false);
    uint64_t length = answerKey.getLength();

    // sample database
    std::vector<std::vector<uint64_t> > data(para.rlwesNum * (0x01 << para.rgswsNum), std::vector<uint64_t>(N, 0));
    sample_spiral_database(data, para.spiral_p);

    // build evaluation keys
    AutoKey autokey(length, modulus, para.ell_rlwes_expansion, para.base_rlwes_expansion, para.bg_rlwes_expansion);
    RGSWCiphertext convertKey(N, modulus, para.ell_convert, para.base_convert, para.bg_convert);
    keygen(autokey, convertKey, answerKey); // evualuation keys
    std::cout << "keyGen Done." << std::endl;

    // query
    RlweCiphertext query(N, modulus);
    queryGenerate(query, answerKey, target_row, target_col, para);

    // answer
    RlweCiphertext result(length, modulus);

    auto start  = std::chrono::high_resolution_clock::now();
    answer(result, query, autokey, convertKey, data, para);
    auto stop  = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " online server response costs " << glapsed.count() << " us." << std::endl;

    // recover message
    std::vector<uint64_t> decryptd_message(N, 0);
    // recover(decryptd_message, result, answerKey);
    recover_noise(decryptd_message, result, answerKey, para.spiral_p);

    std::cout << std::endl;
    showLargeVector(decryptd_message, "The decrypted value is ");

    // print_target_record(data, target_row, target_col, rlwesNum, spiral_p);
    bool isCorrect = check_target_record(data, target_row, target_col, para.rlwesNum, para.spiral_p, decryptd_message);
    if (!isCorrect) std::cout << "Error." << std::endl; 
}

void test_spiral_stream_old()
{
    // total = 2^4 * 2^4
    int32_t rlwesNum = 0x01 << 9; // the number of RLWEs
    int32_t rgswsNum = 6; // the number of RGSWs. note regevToGSW is required
    uint64_t spiral_p = 0x01 << 21; // 0x01 << 16;

    int32_t target_row = rand() % rlwesNum;
    int32_t target_col = rand() % (0x01 << rgswsNum);

    spiral_parameters para;
    get_spiralstream_parameters_one(para, rlwesNum, rgswsNum, spiral_p);
    print_config(para, target_row, target_col);

    uint64_t modulus = crtMod;
    Secret answerKey(modulus, false);
    uint64_t length = answerKey.getLength();

    // sample database
    cout << "sample database ..." << std::endl;
    std::vector<std::vector<uint64_t> > data(para.rlwesNum * (0x01 << para.rgswsNum), std::vector<uint64_t>(N, 0));
    sample_spiral_database(data, spiral_p);

    // query
    cout << "build query ..." << std::endl;
    std::vector<RlweCiphertext> query_rlwes(para.rlwesNum, RlweCiphertext(N, modulus));
    std::vector<RGSWCiphertext> query_rgsws(para.rgswsNum, RGSWCiphertext(N, modulus, para.ellgsw, para.basegsw, para.bggsw));
    queryGenerate(query_rlwes, query_rgsws, answerKey, target_row, target_col, para);

    // answer
    cout << "answer ..." << std::endl;
    RlweCiphertext result(length, modulus);
    answer(result, query_rlwes, query_rgsws, data);

    // recover message
    cout << "recover ..." << std::endl;
    std::vector<uint64_t> decryptd_message(N, 0);
    // recover(decryptd_message, result, answerKey);
    recover_noise(decryptd_message, result, answerKey, para.spiral_p);

    std::cout << std::endl;
    showLargeVector(decryptd_message, "The decrypted value is ");

    bool isCorrect = check_target_record(data, target_row, target_col, para.rlwesNum, para.spiral_p, decryptd_message);
    if (!isCorrect) std::cout << "Error." << std::endl; 
}

void test_spiral()
{
    // total = 2^9 * 2^6
    int32_t rlwesNum = 0x01 << 9; // the number of RLWEs
    int32_t rgswsNum = 6; // the number of RGSWs. note regevToGSW is required
    uint64_t spiral_p = 0x01 << 11; // 0x01 << 16;

    int32_t target_row = rand() % rlwesNum;
    int32_t target_col = rand() % (0x01 << rgswsNum);

    spiral_parameters para;
    get_spiral_parameters_one(para, rlwesNum, rgswsNum, spiral_p);
    print_config(para, target_row, target_col);

    uint64_t modulus = crtMod;
    Secret answerKey(modulus, false);
    uint64_t length = answerKey.getLength();

    // sample database
    std::vector<std::vector<uint64_t> > data(para.rlwesNum * (0x01 << para.rgswsNum), std::vector<uint64_t>(N, 0));
    sample_spiral_database(data, para.spiral_p);

    // handle database
    size_t num_words = para.rlwesNum * (0x01 << para.rgswsNum) * length;
    uint64_t* data_crt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words);
    spiral_sdatabase_tocrt(data_crt, data, para.rlwesNum);

    // build evaluation keys
    AutoKey autokey(length, modulus, para.ell_rlwes_expansion, para.base_rlwes_expansion, para.bg_rlwes_expansion);
    RGSWCiphertext convertKey(N, modulus, para.ell_convert, para.base_convert, para.bg_convert);
    keygen(autokey, convertKey, answerKey); // evualuation keys
    std::cout << "keyGen Done." << std::endl;

    // query
    RlweCiphertext query(N, modulus);
    queryGenerate(query, answerKey, target_row, target_col, para);

    // answer
    RlweCiphertext result(length, modulus);

    auto start  = std::chrono::high_resolution_clock::now();
    answer(result, query, autokey, convertKey, data_crt, para);
    // answer(result, query_rlwes, query_rgsws, data_crt);
    auto stop  = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " online server response costs " << glapsed.count() << " us." << std::endl;

    // recover message
    std::vector<uint64_t> decryptd_message(N, 0);
    // recover(decryptd_message, result, answerKey);
    recover_noise(decryptd_message, result, answerKey, para.spiral_p);

    std::cout << std::endl;
    showLargeVector(decryptd_message, "The decrypted value is ");

    // print_target_record(data, target_row, target_col, rlwesNum, spiral_p);
    bool isCorrect = check_target_record(data, target_row, target_col, para.rlwesNum, para.spiral_p, decryptd_message);
    if (!isCorrect) std::cout << "Error." << std::endl; 
}

void test_spiral_rns()
{
    // total = 2^9 * 2^6
    int32_t rlwesNum = 0x01 << 7; // the number of RLWEs
    int32_t rgswsNum = 8; // the number of RGSWs. note regevToGSW is required
    uint64_t spiral_p = 0x01 << 18; // 0x01 << 16;

    int32_t target_row = rand() % rlwesNum;
    int32_t target_col = rand() % (0x01 << rgswsNum);

    spiral_parameters para;
    get_spiral_rns_parameters_one(para, rlwesNum, rgswsNum, spiral_p);
    print_config(para, target_row, target_col);

    uint64_t modulus = crtMod;
    uint64_t modulus2 = bsMod;
    Secret answerKey(modulus, false);
    uint64_t length = answerKey.getLength();

    // sample database
    std::vector<std::vector<uint64_t> > data(para.rlwesNum * (0x01 << para.rgswsNum), std::vector<uint64_t>(N, 0));
    sample_spiral_database(data, para.spiral_p);

    // handle database
    size_t num_words = para.rlwesNum * (0x01 << para.rgswsNum) * length;
    uint64_t* data_crt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words);
    auto start_pre  = std::chrono::high_resolution_clock::now();
    spiral_sdatabase_tocrt(data_crt, data, para.rlwesNum);
    auto stop_pre  = std::chrono::high_resolution_clock::now();
    auto glapsed_pre = std::chrono::duration_cast<std::chrono::microseconds>(stop_pre - start_pre);
    std::cout << " preprocess costs " << glapsed_pre.count() << " us." << std::endl;

    // build evaluation keys
    // AutoKeyRNS autokey(length, modulus, para.ell_rlwes_expansion, para.base_rlwes_expansion, para.bg_rlwes_expansion);
    AutoKeyRNS autokey(length, modulus, modulus2, para.ell_rlwes_expansion, para.base_rlwes_expansion, para.bg_rgsws_expansion);
    RNSRGSWCiphertext convertKey(N, modulus, para.ell_convert, para.base_convert, para.bg_convert);
    keygen(autokey, convertKey, answerKey); // evualuation keys
    std::cout << "keyGen Done." << std::endl;

    // query
    std::vector<RlweCiphertext> query(1, RlweCiphertext(N, modulus));
    query.push_back(RlweCiphertext(N, modulus2));
    auto start_qu  = std::chrono::high_resolution_clock::now();
    queryGenerate(query, answerKey, target_row, target_col, para);
    auto stop_qu  = std::chrono::high_resolution_clock::now();
    auto glapsed_qu = std::chrono::duration_cast<std::chrono::microseconds>(stop_qu - start_qu);
    std::cout << " query generate costs " << glapsed_qu.count() << " us." << std::endl;

    // answer
    RlweCiphertext result(length, modulus);

    auto start  = std::chrono::high_resolution_clock::now();
    answer(result, query, autokey, convertKey, data_crt, para);
    auto stop  = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " online server response costs " << glapsed.count() << " us." << std::endl;

    // recover message
    std::vector<uint64_t> decryptd_message(N, 0);
    auto start_re  = std::chrono::high_resolution_clock::now();
    // recover(decryptd_message, result, answerKey);
    recover_noise(decryptd_message, result, answerKey, para.spiral_p);
    auto stop_re  = std::chrono::high_resolution_clock::now();
    auto glapsed_re = std::chrono::duration_cast<std::chrono::microseconds>(stop_re - start_re);
    std::cout << " recover costs " << glapsed_re.count() << " us." << std::endl;

    std::cout << std::endl;
    showLargeVector(decryptd_message, "The decrypted value is ");

    // print_target_record(data, target_row, target_col, rlwesNum, spiral_p);
    bool isCorrect = check_target_record(data, target_row, target_col, para.rlwesNum, para.spiral_p, decryptd_message);
    if (!isCorrect) std::cout << "Error." << std::endl; 
}

void test_spiral_pack_rns()
{
    // total = 2^9 * 2^6
    int32_t rlwesNum = 0x01 << 8; // the number of RLWEs
    int32_t rgswsNum = 7; // the number of RGSWs. note regevToGSW is required
    uint64_t spiral_p = 0x01 << 18; // 0x01 << 16;
    int32_t n = 2;

    int32_t target_row = rand() % rlwesNum;
    int32_t target_col = rand() % (0x01 << rgswsNum);

    spiral_parameters para;
    get_spiral_rns_parameters_one(para, rlwesNum, rgswsNum, spiral_p, n);
    print_config(para, target_row, target_col);

    uint64_t modulus = crtMod;
    uint64_t modulus2 = bsMod;
    Secret answerKey(modulus, false);
    uint64_t length = answerKey.getLength();

    // sample database
    std::vector<std::vector<uint64_t> > data(para.rlwesNum * (0x01 << para.rgswsNum) * pow(para.n, 2), std::vector<uint64_t>(N, 0));
    sample_spiral_database(data, para.spiral_p);

    // handle database
    size_t num_words = para.rlwesNum * (0x01 << para.rgswsNum) * length * pow(para.n, 2);
    uint64_t* data_crt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words);
    spiral_sdatabase_tocrt(data_crt, data, para.rlwesNum, para.n);

    // build evaluation keys
    // AutoKeyRNS autokey(length, modulus, para.ell_rlwes_expansion, para.base_rlwes_expansion, para.bg_rlwes_expansion);
    AutoKeyRNS autokey(length, modulus, modulus2, para.ell_rlwes_expansion, para.base_rlwes_expansion, para.bg_rgsws_expansion);
    RNSRGSWCiphertext convertKey(N, modulus, para.ell_convert, para.base_convert, para.bg_convert);
    keygen(autokey, convertKey, answerKey); // evualuation keys
    std::cout << "keyGen Done." << std::endl;

    // query
    std::vector<RlweCiphertext> query(1, RlweCiphertext(N, modulus));
    query.push_back(RlweCiphertext(N, modulus2));
    queryGenerate(query, answerKey, target_row, target_col, para);

    // answer
    std::vector<RlweCiphertext> result(para.n * para.n, RlweCiphertext(length, modulus));

    auto start  = std::chrono::high_resolution_clock::now();
    answer(result, query, autokey, convertKey, data_crt, para);
    auto stop  = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " online server response costs " << glapsed.count() << " us." << std::endl;

    // recover message
    std::vector<uint64_t> decryptd_message(N, 0);
    // recover(decryptd_message, result, answerKey);
    recover_noise(decryptd_message, result[0], answerKey, para.spiral_p);

    std::cout << std::endl;
    showLargeVector(decryptd_message, "The decrypted value is ");

    // print_target_record(data, target_row, target_col, rlwesNum, spiral_p);
    bool isCorrect = check_target_record(data, target_row, target_col, para.rlwesNum, para.spiral_p, decryptd_message);
    if (!isCorrect) std::cout << "Error." << std::endl; 
}

void test_spiral_stream()
{
    // total = 2^4 * 2^4
    int32_t rlwesNum = 0x01 << 9; // 0x01 << 9; // the number of RLWEs
    int32_t rgswsNum = 9; // 6; // the number of RGSWs. note regevToGSW is required
    uint64_t spiral_p = 0x01 << 20; // 0x01 << 16;

    int32_t target_row = rand() % rlwesNum;
    int32_t target_col = rand() % (0x01 << rgswsNum);

    spiral_parameters para;
    get_spiralstream_parameters_one(para, rlwesNum, rgswsNum, spiral_p);
    print_config(para, target_row, target_col);

    uint64_t modulus = crtMod;
    Secret answerKey(modulus, false);
    int32_t length = answerKey.getLength();

    // sample database
    cout << "sample database ..." << std::endl;
    std::vector<std::vector<uint64_t> > data(para.rlwesNum * (0x01 << para.rgswsNum), std::vector<uint64_t>(N, 0));
    sample_spiral_database(data, spiral_p);

    // handle database
    size_t num_words = para.rlwesNum * (0x01 << para.rgswsNum) * length;
    uint64_t* data_crt = (uint64_t *)aligned_alloc(64, sizeof(uint64_t) * num_words);

    auto start_pre   = std::chrono::high_resolution_clock::now();    
    spiral_sdatabase_tocrt(data_crt, data, para.rlwesNum);
    auto stop_pre   = std::chrono::high_resolution_clock::now();
    auto glapsed_pre = std::chrono::duration_cast<std::chrono::microseconds>(stop_pre - start_pre);
    std::cout << " processing costs " << glapsed_pre.count() << " us." << std::endl;

    // query
    cout << "build query ..." << std::endl;
    std::vector<RlweCiphertext> query_rlwes(para.rlwesNum, RlweCiphertext(N, modulus));
    std::vector<RGSWCiphertext> query_rgsws(para.rgswsNum, RGSWCiphertext(N, modulus, para.ellgsw, para.basegsw, para.bggsw));
    auto start_qu  = std::chrono::high_resolution_clock::now();
    queryGenerate(query_rlwes, query_rgsws, answerKey, target_row, target_col, para);
    auto stop_qu   = std::chrono::high_resolution_clock::now();
    auto glapsed_qu = std::chrono::duration_cast<std::chrono::microseconds>(stop_qu - start_qu);
    std::cout << " generate query costs " << glapsed_qu.count() << " us." << std::endl;

    // answer
    cout << "answer ..." << std::endl;
    RlweCiphertext result(length, modulus);
    
    auto start = std::chrono::high_resolution_clock::now();
    answer(result, query_rlwes, query_rgsws, data_crt);
    auto stop  = std::chrono::high_resolution_clock::now();
    auto glapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << " online server response costs " << glapsed.count() << " us." << std::endl;

    // recover message
    cout << "recover ..." << std::endl;
    std::vector<uint64_t> decryptd_message(N, 0);
    auto start_re   = std::chrono::high_resolution_clock::now();
    // recover(decryptd_message, result, answerKey);
    recover_noise(decryptd_message, result, answerKey, para.spiral_p);
    auto stop_re   = std::chrono::high_resolution_clock::now();
    auto glapsed_re = std::chrono::duration_cast<std::chrono::microseconds>(stop_re - start_re);
    std::cout << " recover costs " << glapsed_re.count() << " us." << std::endl;


    std::cout << std::endl;
    showLargeVector(decryptd_message, "The decrypted value is ");

    bool isCorrect = check_target_record(data, target_row, target_col, para.rlwesNum, para.spiral_p, decryptd_message);
    if (!isCorrect) std::cout << "Error." << std::endl; 
}

int main(int argc, char** argv)
{
    srand(time(NULL));

    // test_query_expansion();
    // test_query_expansion_split();
    // test_regev_product();
    // test_spiral_old();
    // test_spiral_stream_old();
    // test_spiral();
    // test_spiral_stream();
    // test_query_expansion_split_rns();
    test_spiral_rns();
    // test_spiral_pack_rns();

    return 0;
}
