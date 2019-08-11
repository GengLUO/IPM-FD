#ifndef __IMP_h__
#define __IMP_h__

typedef unsigned char byte;
#define NULL_BYTE -1

extern byte **L;
extern byte **L_prime;

/**
 * @brief Use the sodium lib RNG to get a random byte
 * 
 * @return byte 
 */
byte random_byte();
unsigned int get_randcount();

/**
 * @brief Setup the masking scheme : choose the best code
 * for n<=4 and k<=2 we chosse the best known code H, overwise the randomly choose the code
 * 
 * @param n number of shares
 * @param k duplicating parameter
 * @return byte** the dual code of H
 */
byte **IPM_FD_Setup(int n, int k);

/**
 * @brief mask a sensitive variable into a sharing w.r.t extended IPM scheme
 * 
 * @param Z result = the sharing
 * @param X sensitive value to be masked
 * @param n nb of shares
 * @param k duplicating parameter
 */
void mask(byte *Z, byte X, int n, int k);

void IPRefresh(byte *Z, int N, int position);
/**
 * @brief refresh the mask
 * 
 * @param Z add to Z a random mask R such that <L,R>=<L,Z>
 * @param n 
 * @param k 
 */
void IPM_FD_Refresh(byte *Z, int n, int k);

/**
 * @brief unmask a sharing
 * We can also check for data integrity (faults)
 * @param Z the sharing
 * @param n nb of shares
 * @param k we have k ways to demask
 * @return the sensitive data or NULL_BYTE if fault detected
 */
int unmask(const byte Z[], int n, int k);

byte GF256_Square(byte x);

byte GF256_Mult(byte a, byte b);
/**
 * @brief Secure multiplication with fault detection
 * 
 * @param R = Z x Z' masked under the same L
 * @param Z 1st operand
 * @param Z_prime 2nd operand
 */
void IPM_FD_Mult(byte *R, const byte *Z, const byte *Z_prime, int n, int k);
void IPM_FD_Add(byte *res, const byte *op1, const byte *op2, int n);

void IPM_FD_Square(byte *x2, const byte *x, int n, int k);

/**
 * @brief compute squaring more efficiently to use in S-box
 * Use original IPSquare 
 * 
 * @param x2 
 * @param x 
 * @param N 
 * @param position indicate witch L to use
 */
void IPM_Square(byte *x2, const byte *x, int N, int position);

void innerProduct(byte *Z, const byte *M, int n, int k);
/**
 * @brief Homogenization of two sharings
 * 
 * @param position indicate whitch L to use w.r.t the the matrix L in extended IPM
 * @return an homegenized value of T and T'
 * @param N = n-k+1
 * @note : We improved this algo by only computing the needed value instead of the whole table R
 */
byte IPM_Homogenize(const byte *Z, const byte *Z_prime, int N, int k, int position);

/**
 * @brief standard IPM multiplication
 * multiply masked values under same vector L
 * 
 * @param P result of mutl
 * @param Z 1st operand
 * @param Z_prime 2nd operand
 * @param N size of operands, also equals to n-k+1, (n is the original number of shares)
 * @param k current duplicating parameter
 * @param position indicate with L to use w.r.t the the matrix L in extended IPM
 */
void IPM_Mult(byte *P, const byte *Z, const byte *Z_prime, int N, int k, int position);
void IPConstMult(byte *res, const byte *x, byte c, int n);
void IPConstAdd(byte *res, const byte *x, byte c, int n, int k);

void error(const char *msg);

byte ***allocate3D(int k, int m, int n);
byte **allocate2D(int rows, int cols);
void deallocate2D(byte **arr2D, int rows, int cols);
void deallocate3D(byte ***arr3D, int k, int m, int p);

void print(const byte *a, const char *msg, int n);

void freeMemory(int n, int k);
#endif
