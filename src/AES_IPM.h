#ifndef __AES_IPM_h__
#define __AES_IPM_h__

//#include "IPM.h"
typedef unsigned char byte;

void sbox_share(byte *res, const byte *x, int n, int k);

void printMes(const char *s, const byte *m);
void unmaskPrintMes(byte **stateshare, const char *s, int n, int k);

void x127_share(byte *res, const byte *x, int n, int k);
void AES_IPM_FD_Key_Expansion(const byte *key, byte *wshare[176], int n, int k);

int AES_IPM_FD(const byte in[16], byte out[16], const byte key[16], int n, int k);
int run_aes_share(const byte in[16], byte out[16], const byte key[16], int n, int k, int *status, int dt);
void inv_aes_share(const byte in[16], byte out[16], const byte key[16], int n, int k);
#endif
