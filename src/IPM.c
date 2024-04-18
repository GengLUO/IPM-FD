// Implement Extended Inner Product Masking with faults detection : See PROOF submission (http://www.proofs-workshop.org/2019/)

// In all the files of this project
// n is the number of shares
// k is duplicating parameter for faults detection ( n > k)
// N = n - k + 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sodium.h> //for random numbers

#include "aes.h"
#include "IPM.h"

//precomputed squaring table in GF(256) GF(2^8)
byte sq[256] = {
        0x00, 0x01, 0x04, 0x05, 0x10, 0x11, 0x14, 0x15,
        0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55,
        0x1b, 0x1a, 0x1f, 0x1e, 0x0b, 0x0a, 0x0f, 0x0e,
        0x5b, 0x5a, 0x5f, 0x5e, 0x4b, 0x4a, 0x4f, 0x4e,
        0x6c, 0x6d, 0x68, 0x69, 0x7c, 0x7d, 0x78, 0x79,
        0x2c, 0x2d, 0x28, 0x29, 0x3c, 0x3d, 0x38, 0x39,
        0x77, 0x76, 0x73, 0x72, 0x67, 0x66, 0x63, 0x62,
        0x37, 0x36, 0x33, 0x32, 0x27, 0x26, 0x23, 0x22,
        0xab, 0xaa, 0xaf, 0xae, 0xbb, 0xba, 0xbf, 0xbe,
        0xeb, 0xea, 0xef, 0xee, 0xfb, 0xfa, 0xff, 0xfe,
        0xb0, 0xb1, 0xb4, 0xb5, 0xa0, 0xa1, 0xa4, 0xa5,
        0xf0, 0xf1, 0xf4, 0xf5, 0xe0, 0xe1, 0xe4, 0xe5,
        0xc7, 0xc6, 0xc3, 0xc2, 0xd7, 0xd6, 0xd3, 0xd2,
        0x87, 0x86, 0x83, 0x82, 0x97, 0x96, 0x93, 0x92,
        0xdc, 0xdd, 0xd8, 0xd9, 0xcc, 0xcd, 0xc8, 0xc9,
        0x9c, 0x9d, 0x98, 0x99, 0x8c, 0x8d, 0x88, 0x89,
        0x9a, 0x9b, 0x9e, 0x9f, 0x8a, 0x8b, 0x8e, 0x8f,
        0xda, 0xdb, 0xde, 0xdf, 0xca, 0xcb, 0xce, 0xcf,
        0x81, 0x80, 0x85, 0x84, 0x91, 0x90, 0x95, 0x94,
        0xc1, 0xc0, 0xc5, 0xc4, 0xd1, 0xd0, 0xd5, 0xd4,
        0xf6, 0xf7, 0xf2, 0xf3, 0xe6, 0xe7, 0xe2, 0xe3,
        0xb6, 0xb7, 0xb2, 0xb3, 0xa6, 0xa7, 0xa2, 0xa3,
        0xed, 0xec, 0xe9, 0xe8, 0xfd, 0xfc, 0xf9, 0xf8,
        0xad, 0xac, 0xa9, 0xa8, 0xbd, 0xbc, 0xb9, 0xb8,
        0x31, 0x30, 0x35, 0x34, 0x21, 0x20, 0x25, 0x24,
        0x71, 0x70, 0x75, 0x74, 0x61, 0x60, 0x65, 0x64,
        0x2a, 0x2b, 0x2e, 0x2f, 0x3a, 0x3b, 0x3e, 0x3f,
        0x6a, 0x6b, 0x6e, 0x6f, 0x7a, 0x7b, 0x7e, 0x7f,
        0x5d, 0x5c, 0x59, 0x58, 0x4d, 0x4c, 0x49, 0x48,
        0x1d, 0x1c, 0x19, 0x18, 0x0d, 0x0c, 0x09, 0x08,
        0x46, 0x47, 0x42, 0x43, 0x56, 0x57, 0x52, 0x53,
        0x06, 0x07, 0x02, 0x03, 0x16, 0x17, 0x12, 0x13};

byte hardrandom[4][4] = {
        { 43, 65, 63, 97 },
        { 123, 1, 239, 54 },
        { 78, 76, 127, 179 },
        { 222, 48, 74, 59 }
};
//IPM settings : constant public valuesTODO:important
byte **L;       // The orthogonal of n^2 matrix H
byte **L_prime; // orthogonal of H with the (k - 1) 0s removed (n - k +1)^2
byte **L_prime_inv; //Inverse of L_prime, used for IPM_MULT
byte L_prime_test_homo [2][4]= {
        {1, 19, 205, 142},
        {1,83,65,68}
};

//number of random numbers generated
static unsigned int randcount = 0;
static unsigned int multcount = 0;

//generate random non zero byte using sodium library
byte random_byte()  ///randNonZero()
{
    randcount++;
    byte b = randombytes_uniform(256);
    while (b == 0)
    {
        b = randombytes_uniform(256);
    }
    return b;
}

//return the current number of used random numbers
unsigned int get_randcount()
{
    return randcount;
}

/**
 * @brief Setup the IPM-FD masking scheme: TODO:important
 *        1. L
 *        2. L_prime
 *        3. L_hat
 *
 * @param n : the number of shares
 * n: Length of each codeword (number of symbols per codeword).
 * k: Dimension of the code (number of information symbols each codeword can represent).
 * d: Minimum distance (the error-correcting capability of the code).
 * @patam k duplicating parameter for faults detection ( n > k)
 */
byte **IPM_FD_Setup(int n, int k)
{
    if (n <= 0 || k <= 0 || n < k)
        error("n > 0 should be greater than k > 0\n");
    int i = 0;
    randcount = 0;
    while (sodium_init() < 0) //init sodium library
    {
        i++;
        if (i >= 5)
            error("panic! the sodium library couldn't be initialized, it is not safe to use\n");
    }
    L = allocate2D(k, n);
    int j;
    //use best known linear code for L (aka dual of H)

    if (n == 1 && k == 1)
    {
        L[0][0] = 1;
    }
    else if (n == 2 && k == 1)
    { //best known L for n=2 and k=1
        // 27 = X^8 = X ^ 4 + X ^ 3 + X + 1
        L[0][0] = 1;
        L[0][1] = 27;
    }
    else if (n == 3 && k == 1)
    { // 250 = X^26 = X^7 + X^6 + X^5 + X^4 + X^3 + X
        L[0][0] = 1;
        L[0][1] = 27;
        L[0][2] = 250;
    }
    else if (n == 4 && k == 1)
    { //188 = X^17 = X^7 + X^5 + X^4 + X^3 + X^2
        L[0][0] = 1;
        L[0][1] = 27;
        L[0][2] = 250;
        L[0][3] = 188;
    }
    else if (n == 2 && k == 2)
    {
        L[0][0] = 1;
        L[0][1] = 0;

        L[1][0] = 0;
        L[1][1] = 1;
    }
    else if (n == 3 && k == 2)
    { //best known L for n=3 and k=2
        L[0][0] = 1;
        L[0][1] = 0;
        L[0][2] = 27;

        L[1][0] = 0;
        L[1][1] = 1;
        L[1][2] = 188;
    }
    else if (n == 4 && k == 2)
    { //151 = x^20 = X^7 + X^4 + X^2 + X + 1
        //239 = x^27 = X^7 + X^6 + X^5 + X^3 + X^2 + X + 1
        // 128 = X^7
        L[0][0] = 1;
        L[0][1] = 0;
        L[0][2] = 27;
        L[0][3] = 151;

        L[1][0] = 0;
        L[1][1] = 1;
        L[1][2] = 239;
        L[1][3] = 128;
    }
    else if (n==5 && k == 2) //TODO: DELETE IT
    {
        L[0][0] = 1;
        L[0][1] = 0;
        L[0][2] = 43;
        L[0][3] = 122;
        L[0][4] = 199;

        L[1][0] = 0;
        L[1][1] = 1;
        L[1][2] = 27;
        L[1][3] = 250;
        L[1][4] = 188;
    }
    else
    {
        for (i = 0; i < k; i++)
        {
            for (j = 0; j < k; j++)
                L[i][j] = (i == j) ? 1 : 0;

            //choose the best values here
            for (j = k; j < n; j++)
                L[i][j] = random_byte(); //make sure they all non-zero
        }
    }
    int N = n - k + 1;
    L_prime = allocate2D(k, N);
    L_prime_inv = allocate2D(k, N);

    for (j = 0; j < k; j++)
    {
        L_prime[j][0] = 1;
        L_prime_inv[j][0] = 1;
        for (i = 1; i < N; i++)
        {
            L_prime[j][i] = L[j][i + k - 1];
            L_prime_inv[j][i] = GF256_Inverse(L[j][i + k - 1]);
        }
    }

    printf("L:\n");
    for (i = 0; i<k; i++){
        for (j = 0; j < n; j++){
            printf("%d,",L[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("L':\n");
    for (i = 0; i<k; i++){
        for (j = 0; j < N; j++){
            printf("%d,",L_prime[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("L'inv:\n");
    for (i = 0; i<k; i++){
        for (j = 0; j < N; j++){
            printf("%d,",L_prime_inv[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    return L;
}

//generate the squaring table for lookup
void gen_square_table()
{
    int i;
    byte x = 0;
    printf("byte sq[256]={");
    for (i = 0; i < 256; i++)
    {
        if ((i % 8) == 0)
            printf("\n");
        printf("0x%02x", GF256_Mult(x, x));
        x++;
        if (i < 255)
            printf(",");
    }
    printf("};\n");
}

byte GF256_Square(byte x)
{
    //return GF256_Mult(x,x);
    return sq[x];
}

//Iterative constant time multiplication over GF(2^8) with unconditional reduction
//see paper "On the Performance and Security of Multiplication in GF(2^N)"
byte GF256_Mult(byte a, byte b)
{
    int x = a, y = b, m, res = 0;
    unsigned char i;
    for (i = 0; i < 8; i++)
    {
        // Determine if the least significant bit of y is 1 (using bitwise AND and negation)
        m = -(y & 1); // m is either 0xffff (if y's LSB is 1 -> -(y & 1) = -1 = 0xffff in 2's compliment) or 0x0000 (if y's LSB is 0)
        // Conditional XOR based on the LSB of y
        // If y's LSB is 1 (m = 0xffff), res is XORed with x. If y's LSB is 0 (m = 0x0000), res remains unchanged.
        // On each iteration of the loop, res accumulates the partial product of the multiplication. This accumulation is done through the XOR operation.
        // This approach mimics the addition step in traditional multiplication but is done in a bitwise manner suitable for field arithmetic in GF(2^8).
        res ^= (x & m);
        // Right shift y, preparing for the next iteration (next bit)
        y >>= 1;
        // Left shift x, equivalent to multiplying by 2 in GF(2^8)
        x <<= 1;
        // Reduction step: if the result of the shift left operation on x
        // results in a byte with more than 8 bits (i.e., if the MSB is 1),
        // reduce it using the irreducible polynomial (0x1b).
        m = -((x >> 8) & 1); // Check MSB after the shift (if overflow or not)
        // If x after the shift has a MSB of 1 (i.e., x > 255), it's XORed with 0x1b for reduction.
        // If the MSB is 0 (i.e., x <= 255), no change is made to x.
        x ^= (m & 0x1b);
    }
    multcount++;
    return (byte)res;
}

void IPM_Square(byte *x2, const byte *x, int N, int position)
{
    int i;
    for (i = 1; i < N; i++) {
        x2[i] = GF256_Mult(GF256_Square(x[i]), L_prime[position][i]);
//        printf("%2x,",GF256_Square(x[i]));
    }
    x2[0] = GF256_Square(x[0]); //L_prime[position][0] = 1

    printf("Square position = %d\n", position);
    print(x, "x: ", N);
    print(x2, "x^2: ", N);
    printf("\n");
}

//squaring of IPM_FD share, more efficient than Mult(x,x)
void IPM_FD_Square(byte *Z2, const byte *Z, int n, int k)
{
    int i, j, N = n - k + 1;
    byte P[k][N];
    byte Q[k - 1];

    byte Z__[k][N];

    for (i = 0; i < k; i++)
    {
        Z__[i][0] = Z[i];
        memcpy(Z__[i] + 1, Z + k, n - k); //copy first i-1 coord of Z
    }
    for (j = 0; j < k; j++)
        IPM_Square(P[j], Z__[j], N, j);

    //No need to refresh because IPM_Mult already does the job
    //for (j = 1; j < k; j++) //refresh only k-1 of them
    //    IPRefresh(P[j], N, j);

    for (j = 0; j < k - 1; j++)
        Q[j] = IPM_Homogenize(P[0], P[j + 1], N, k, j + 1);

    //return
    Z2[0] = P[0][0];
    memcpy(Z2 + k, P[0] + 1, n - k);
    for (j = 0; j < k - 1; j++)
    {
        Z2[j + 1] = Q[j];
    }
}

//Algo 2 : secure addition in IPM-FD √
void IPM_FD_Add(byte *res, const byte *op1, const byte *op2, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        res[i] = op1[i] ^ op2[i];
    }
}

void IPM_FD_Mult(byte *R, const byte *Z, const byte *Z_prime, int n, int k)
{
    int i, j, N = n - k + 1;
    byte P[k][N];

    byte Z__[k][N], Z_prime__[k][N]; //Z__[i] => i th in [0...k-1] coord of Z is kept, the other in [0...k-1] are dropped;

    for (i = 0; i < k; i++)
    {
        Z__[i][0] = Z[i];
        memcpy(Z__[i] + 1, Z + k, n - k); //copy first i-1 coord of Z

        Z_prime__[i][0] = Z_prime[i];
        memcpy(Z_prime__[i] + 1, Z_prime + k, n - k);
    }
    multcount = 0;
    for (j = 0; j < k; j++)
        IPM_Mult(P[j], Z__[j], Z_prime__[j], N, k, j);
//    printf("multcount = %d\n",multcount);

    for (j = 1; j < k; j++)
        P[j][0] = IPM_Homogenize(P[0], P[j], N, k, j);

    //return
    for (j = 0; j < k; j++)
        R[j] = P[j][0];
    memcpy(R + k, P[0] + 1, n - k);

}

//Add constant value to masked vector
void IPConstAdd(byte *res, const byte *x, byte c, int n, int k)
{
    memcpy(res + k, x + k, n - k);
    int i;
    for (i = 0; i < k; i++)
        res[i] = x[i] ^ c;
}

//Multiply constant value to masked vector
void IPConstMult(byte *res, const byte *x, byte c, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        res[i] = GF256_Mult(x[i], c);
    }
}

void IPM_Mult(byte *P, const byte *Z, const byte *Z_prime, int N, int k, int position)
{
//    byte T[N][N], U[N][N], V[N][N], U_prime[N][N];

    byte T, U, U_prime;
    byte next;
    byte R[N];
    memset(R,0,N);

    int i, j;
    int index_i, index_j;

    // Computation of the matrix T
    for (i = 0; i < N; i++) {
        for (j = i; j < N; (next)?j++:j) {

            //1. Compute T
            if (i == j){
                U_prime = 0;
                next = 1;

                index_i = i;
                index_j = j;
            }
            else {
                if (next) {
                    U_prime = hardrandom[i][j];
//                    U_prime = random_byte();
                    next = 0;

                    index_i = i;
                    index_j = j;
                } else
                {
                    next = 1;

                    index_i = j;
                    index_j = i;
                }
            }
            T = GF256_Mult(GF256_Mult(Z[index_i], Z_prime[index_j]) , L_prime[position][index_j]);
            U = GF256_Mult(U_prime, L_prime_inv[position][index_i]);
            R[index_i] ^= T ^ U;
        }

    }
    memcpy(P,R,N);

    printf("Mult position = %d\n", position);
    print(Z, "Z: ", N);
    print(Z_prime, "Z_prime: ", N);
    print(P, "P: ",N);
    printf("\n");
//    printf("P':\n");
//    for (int i = 0; i < N; i++) {
//        printf("%2x",P[i]);
//    }
//    printf("\n");
}

//void IPM_Mult(byte *P, const byte *Z, const byte *Z_prime, int N, int k, int position)
//{
//
//    //we have here n-k+1 compare to the original n in IPMsetup()
//    int i, j;
//
//    //step 1 -> get the orthogonal of L : A_hat
//    byte A_hat[N][N];
//    for (i = 0; i < N - 1; i++)
//        randombytes(A_hat[i], N);
//    randombytes(A_hat[N - 1], N - 1);
//
//    randcount += N * N - 1;
//    byte delta = 0;
//    for (j = 0; j < N - 1; j++)
//    {
//        byte sum = 0;
//        for (i = 0; i < N; i++)
//        {
//            sum ^= GF256_Mult(A_hat[i][j], L_hat[position][i][j]);
//        }
//        delta ^= sum;
//    }
////start with n*(n-1) multiplication
//    byte sum = 0;
//    for (i = 0; i < N - 1; i++) //row-wise
//    {
//        sum ^= GF256_Mult(A_hat[i][N - 1], L_hat[position][i][N - 1]);
//    }
////then (n-1) multiplication
//    A_hat[N - 1][N - 1] = GF256_Mult(delta ^ sum, GF256_Inverse(L_hat[position][N - 1][N - 1]));
////then 1 multiplication
//
////need to put together an orthogonal matrix -> sequential, a lot of MAC
//
//    //step 2 -> tensor product of R and Q : R_hat
//    byte R_hat[N][N];
//    //save values for i,j in {k, N}
//
//    for (i = 0; i < N; i++)
//    {
//        for (j = 0; j < N; j++)
//            R_hat[i][j] = GF256_Mult(Z[i], Z_prime[j]);
//    }
//
//    //step 3 -> addition of R_hat and A_hat
//    byte B_hat[N][N];
//    for (i = 0; i < N; i++)
//    {
//        for (j = 0; j < N; j++)
//        {
//            B_hat[i][j] = R_hat[i][j] ^ A_hat[i][j];
//        }
//    }
//    //step 4 -> get the b :
//    byte b = 0;
//    for (i = 1; i < N; i++)
//    {
//        byte tmp = 0;
//        for (j = 0; j < N; j++)
//        {
//            tmp ^= GF256_Mult(L_hat[position][i][j], B_hat[i][j]);
//        }
//        b ^= tmp;
//    }
//    //return
//    memcpy(P, B_hat[0], N);
//    P[0] ^= b;
//}

byte IPM_Homogenize(const byte *Z, const byte *Z_prime, int N, int k, int position)
{
    int i;
    byte res = Z_prime[0];
    // byte epsilon;
    for (i = 1; i < N; i++)
    {
        res ^= GF256_Mult(L_prime[position][i], Z[i] ^ Z_prime[i]);
    }
    printf("Homo position = %d\n", position);
    print(Z, "Z: ", N);
    print(Z_prime, "Z_prime: ", N);
    printf("Homo_result: %02x",res);
    printf("\n\n");
    return res;
}

/**
 * @brief inner product between the mask M and param L
 *
 * @param Z  = sum (L<i+k> * Mi)
 * @param M (the mask) has size = n-k
 */
void innerProduct(byte *Z, const byte *M, int n, int k)
{
    int i, j;
    for (j = 0; j < k; j++) // Left part of H: the random-only part of L
    {
        Z[j] = 0;
        for (i = 0; i < n - k; i++)
            Z[j] ^= GF256_Mult(L[j][i + k], M[i]); //Field multiplication
    }
    for (i = k; i < n; i++) // right part of H: just I_n-k
        Z[i] = M[i - k];
}

void mask(byte *Z, byte X, int n, int k)
{
    int i;
    byte M[n - k]; //n - k random masks
//    randombytes(M, n - k); //Get random masks, as in Equation 7 : M_k+1, ... , M_n
//    randcount += n - k;
    M[0] = 65;
    M[1] = 63;
    M[2] = 97;
    innerProduct(Z, M, n, k);
    for (i = 0; i < k; i++)
        Z[i] ^= X; //duplicate the data for k times
    //The first k个 Z[i] = X + last (n-k)个(L * M)
    //End up with k IPM sharings
    //Have k ways to demask
    //<L_1, Z> = X
    //<L_2, Z> = X
    //which means the Z can be used with first k L to get the X
}

// /**
//  * @brief mask a sensitive variable into a sharing w.r.t IPM-FD scheme
//  *
//  * @param Z result = the sharing
//  * @param X sensitive value to be masked
//  * @param n nb of shares
//  * @param k duplicating parameter
//  */
// void mask(byte *Z, byte X, int n, int k)
// {
//     int i;
//     byte M[n - k]; //n - k random masks
//     randombytes(M, n - k);
//     randcount += n - k; //number of random bytes used

//     for (int j = 0; j < k; j++)
//     {
//         Z[j] = 0;
//         for (i = 0; i < n - k; i++)
//             Z[j] ^= GF256_Mult(L[j][i + k], M[i]); //Field multiplication
//         Z[j] ^= X;                                 //Duplicate the data
//     }
//     for (i = k; i < n; i++)
//         Z[i] = M[i - k]; //Add the mask
// }

//return 1 if detected faults and 0 overwise
//TODO: return where the fault is
int detect_faults(const byte *X, int k)
{
    int i, n = 0;
    for (i = 1; i < k; i++)
        if (X[0] != X[i])
        {
            // printf("***Fault deteted : %d, %d\n", X[0], X[1]);
            n = 1;
            break;
        }
    return n;
}

//conditional swap
//return a if bool==1, b overwise
int cswap(int a, int b, int boo)
{
    return boo * a + (1 - boo) * b;
}

int unmask(const byte Z[], int n, int k)
{
    int i, j;
    byte res[k];
    memset(res, 0, k);
    for (j = 0; j < k; j++)
    {
        for (i = 0; i < n; i++)
            res[j] ^= GF256_Mult(L[j][i], Z[i]);
    }
    //check for faults
    i = detect_faults(res, k);
    //return NULL_BYTE if i=1 (fault detected) and correct byte otherwise
    return cswap(NULL_BYTE, res[0], i);
}

// int unmask(const byte *Z, int n, int k)
// {
//     int i, j;
//     byte res[k];
//     memset(res, 0, k);
//     for (j = 0; j < k; j++)
//     {
//         for (i = 0; i < n; i++)
//             res[j] ^= GF256_Mult(L[j][i], Z[i]); //Field multiplication
//     }

//     i = 0;
//     for (j = 1; j < k; j++) //check for faults
//         if (res[0] != res[j])
//         {
//             i = 1; //Fault deteted
//             break;
//         }

//     //return NULL_BYTE = -1 if i = 1 (fault detected) and the correct byte otherwise
//     return cswap(NULL_BYTE, res[0], i);
// }

void IPRefresh(byte *Z, int N, int position)
{
    byte R[N];
    randombytes(R + 1, N - 1); //n-k random bytes
    randcount += N - 1;
    int i;
    R[0] = 0;
    for (i = 1; i < N; i++)
        R[0] ^= GF256_Mult(R[i], L_prime[position][i]);

    IPM_FD_Add(Z, Z, R, N);
}
//Algo 3 √
void IPM_FD_Refresh(byte *Z, int n, int k)
{
    int i, j;
    byte epsilon[n - k], Z_prime[n];
    randombytes(epsilon, n - k);
    randcount += n - k;

    for (i = k; i < n; i++)
        Z_prime[i] = Z[i] ^ epsilon[i - k];

    for (i = 0; i < k; i++)
    {
        Z_prime[i] = Z[i];
        for (j = 0; j < n - k; j++)
            Z_prime[i] ^= GF256_Mult(L[i][j + k], epsilon[j]); //Field multiplication
    }
    memcpy(Z, Z_prime, n);
}

void error(const char *msg)
{
    perror(msg);
    exit(-1);
}

//allocate a 3D array
byte ***allocate3D(int k, int m, int n)
{
    byte ***arr3D;
    int i, j;

    arr3D = (byte ***)malloc(k * sizeof(byte **));
    if (!arr3D)
        error("malloc failed at allocate3D\n");
    for (i = 0; i < k; i++)
    {
        arr3D[i] = (byte **)malloc(m * sizeof(byte *));
        if (!arr3D[i])
            error("malloc failed at allocate3D\n");
        for (j = 0; j < m; j++)
        {
            arr3D[i][j] = (byte *)malloc(n * sizeof(byte));
            if (!arr3D[i][j])
                error("malloc failed at allocate3D\n");
        }
    }

    return arr3D;
}

//deallocate a 3D array
void deallocate3D(byte ***arr3D, int k, int m, int p)
{
    int i, j;

    for (i = 0; i < k; i++)
    {
        for (j = 0; j < m; j++)
        {
            memset(arr3D[i][j], 0, p);
            free(arr3D[i][j]);
        }
        memset(arr3D[i], 0, m);
        free(arr3D[i]);
    }
    memset(arr3D, 0, k);
    free(arr3D);
}

//allocate a 2D array
byte **allocate2D(int rows, int cols)
{
    byte **arr2D;
    int i;

    arr2D = (byte **)malloc(rows * sizeof(byte *));
    for (i = 0; i < rows; i++)
    {
        arr2D[i] = (byte *)malloc(cols * sizeof(byte));
    }
    return arr2D;
}

//deallocate a 2D array
void deallocate2D(byte **arr2D, int rows, int cols)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        memset(arr2D[i], 0, cols);
        free(arr2D[i]);
    }
    memset(arr2D, 0, rows);
    free(arr2D);
}

void freeMemory(int n, int k)
{
    deallocate2D(L, k, n);
    deallocate2D(L_prime, k, n);
    deallocate2D(L_prime_inv, k, n);
}
void print(const byte *a, const char *msg, int n)
{
    int j;
    printf("%s[", msg);
    for (j = 0; j < n - 1; j++)
    {
        printf("%02x,", a[j]);
    }
    printf("%02x]\n", a[n - 1]);
}