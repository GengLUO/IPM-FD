#ifndef __common_h__
#define __common_h__

typedef unsigned char byte;

void report_time(int dt, int nt, int base, unsigned int nrand);
void check_ciphertext(byte *out, byte *outex, int nbyte);
int runalgo(void (*algo)(const byte *, byte *, const byte *, int, int), byte *in, byte *out, byte *key, int n, int l,
            byte *outex, int nbyte, int nt, int base);
#endif
