//secure AES implemention with IPM faults detection

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "aes.h"
#include "IPM.h"
#include "AES_IPM.h"

#define CIPHER_VERBOSE 0
#define DECIPHER_VERBOSE 0

//raise sharing to the power 127
void x127_share(byte *res, const byte *x, int n, int k)
{
	// use 6 squarings and 4 multiplications, can we do better?
	byte z[n];
	IPM_FD_Square(z, x, n, k); // z=x^2

	byte y[n];
	IPM_FD_Mult(y, x, z, n, k); // y=z*x=x^3

	byte w[n];
	IPM_FD_Square(w, y, n, k); //6

	IPM_FD_Mult(z, w, x, n, k); // z =x^7
	IPM_FD_Square(w, w, n, k);  // w=x^12

	byte y2[n];
	IPM_FD_Mult(y2, w, y, n, k); // y2=x^15

	IPM_FD_Square(y2, y2, n, k); //30
	IPM_FD_Square(y2, y2, n, k); //60
	IPM_FD_Square(y2, y2, n, k); //120

	IPM_FD_Mult(res, y2, z, n, k); // y=x^127
}

//perform the sbox transformation on masked vector
//using sbox polynomial expression
//63 + 8fx^127 + b5x^191 + 01x^223 + f4x^239 + 25x^247 + f9x^251 + 09x^253 + 05x^254
void sbox_share(byte *res, const byte *x, int n, int k)
{
	byte x254[n], x253[n], x251[n], x247[n], x239[n], x223[n], x191[n], x127[n];

	x127_share(x127, x, n, k);
	IPM_FD_Square(x254, x127, n, k);

	IPM_FD_Square(x253, x254, n, k);
	IPM_FD_Square(x251, x253, n, k); //251 = (253* 2)-255
	IPM_FD_Square(x247, x251, n, k);

	IPM_FD_Square(x239, x247, n, k);
	IPM_FD_Square(x223, x239, n, k); //x223 = (x239)^2
	IPM_FD_Square(x191, x223, n, k); //x^191 = (x^223)^2

	IPConstMult(x254, x254, 0x05, n);
	IPConstMult(x253, x253, 0x09, n);
	IPConstMult(x251, x251, 0xf9, n);

	IPConstMult(x247, x247, 0x25, n);
	IPConstMult(x239, x239, 0xf4, n);
	//IPConstMult(x223, x223,0x01, n);
	IPConstMult(x191, x191, 0xb5, n);
	IPConstMult(x127, x127, 0x8f, n);

	IPConstAdd(res, x127, 0x63, n, k);
	for (int i = 0; i < n; i++)
		res[i] ^= x191[i] ^ x223[i] ^ x239[i] ^ x247[i] ^ x251[i] ^ x253[i] ^ x254[i];
}

//subbyte operation on masked state
void subbyte_share(byte *stateshare[16], int n, int k)
{
	int j;
	for (j = 0; j < 16; j++)
		sbox_share(stateshare[j], stateshare[j], n, k);
}

//TODO complete
void inv_subbyte_share(byte *stateshare[16], int n)
{
	for (int j = 0; j < 16; j++)
	{
		//find the inv S-box polynomial form
	}
}

void shiftrows_share(byte *stateshare[16], int n)
{
	byte m;
	int i;
	for (i = 0; i < n; i++)
	{
		m = stateshare[1][i];
		stateshare[1][i] = stateshare[5][i];
		stateshare[5][i] = stateshare[9][i];
		stateshare[9][i] = stateshare[13][i];
		stateshare[13][i] = m;

		m = stateshare[2][i];
		stateshare[2][i] = stateshare[10][i];
		stateshare[10][i] = m;
		m = stateshare[6][i];
		stateshare[6][i] = stateshare[14][i];
		stateshare[14][i] = m;

		m = stateshare[3][i];
		stateshare[3][i] = stateshare[15][i];
		stateshare[15][i] = stateshare[11][i];
		stateshare[11][i] = stateshare[7][i];
		stateshare[7][i] = m;
	}
}

void inv_shiftrows_share(byte *stateshare[16], int n)
{
	byte m;
	int i;
	for (i = 0; i < n; i++)
	{
		m = stateshare[13][i];

		stateshare[13][i] = stateshare[9][i];
		stateshare[9][i] = stateshare[5][i];
		stateshare[5][i] = stateshare[1][i];
		stateshare[1][i] = m;

		m = stateshare[2][i];
		stateshare[2][i] = stateshare[10][i];
		stateshare[10][i] = m;
		m = stateshare[6][i];
		stateshare[6][i] = stateshare[14][i];
		stateshare[14][i] = m;

		m = stateshare[15][i];
		stateshare[15][i] = stateshare[3][i];
		stateshare[3][i] = stateshare[7][i];
		stateshare[7][i] = stateshare[11][i];
		stateshare[11][i] = m;
	}
}

void mixcolumns_share(byte *stateshare[16], int n)
{
	byte ns[16];
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 4; j++)
		{
			ns[j * 4] = multx(stateshare[j * 4][i]) ^ multx(stateshare[j * 4 + 1][i]) ^ stateshare[j * 4 + 1][i] ^ stateshare[j * 4 + 2][i] ^ stateshare[j * 4 + 3][i];
			ns[j * 4 + 1] = stateshare[j * 4][i] ^ multx(stateshare[j * 4 + 1][i]) ^ multx(stateshare[j * 4 + 2][i]) ^ stateshare[j * 4 + 2][i] ^ stateshare[j * 4 + 3][i];
			ns[j * 4 + 2] = stateshare[j * 4][i] ^ stateshare[j * 4 + 1][i] ^ multx(stateshare[j * 4 + 2][i]) ^ multx(stateshare[j * 4 + 3][i]) ^ stateshare[j * 4 + 3][i];
			ns[j * 4 + 3] = multx(stateshare[j * 4][i]) ^ stateshare[j * 4][i] ^ stateshare[j * 4 + 1][i] ^ stateshare[j * 4 + 2][i] ^ multx(stateshare[j * 4 + 3][i]);
		}
		for (j = 0; j < 16; j++)
			stateshare[j][i] = ns[j];
	}
}

void inv_mixcolumns_share(byte *stateshare[16], int n)
{
	byte ns[16];
	int i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 4; j++)
		{
			ns[j * 4 + 0] = GF256_Mult(stateshare[j * 4][i], 0x0e) ^ GF256_Mult(stateshare[j * 4 + 1][i], 0x0b) ^ GF256_Mult(stateshare[j * 4 + 2][i], 0x0d) ^ GF256_Mult(stateshare[j * 4 + 3][i], 0x09);
			ns[j * 4 + 1] = GF256_Mult(stateshare[j * 4][i], 0x09) ^ GF256_Mult(stateshare[j * 4 + 1][i], 0x0e) ^ GF256_Mult(stateshare[j * 4 + 2][i], 0x0b) ^ GF256_Mult(stateshare[j * 4 + 3][i], 0x0d);
			ns[j * 4 + 2] = GF256_Mult(stateshare[j * 4][i], 0x0d) ^ GF256_Mult(stateshare[j * 4 + 1][i], 0x09) ^ GF256_Mult(stateshare[j * 4 + 2][i], 0x0e) ^ GF256_Mult(stateshare[j * 4 + 3][i], 0x0b);
			ns[j * 4 + 3] = GF256_Mult(stateshare[j * 4][i], 0x0b) ^ GF256_Mult(stateshare[j * 4 + 1][i], 0x0d) ^ GF256_Mult(stateshare[j * 4 + 2][i], 0x09) ^ GF256_Mult(stateshare[j * 4 + 3][i], 0x0e);
		}
		for (j = 0; j < 16; j++)
			stateshare[j][i] = ns[j];
	}
}

void addroundkey_share(byte *stateshare[16], byte *wshare[176], int round, int n)
{
	int i, j;
	for (i = 0; i < 16; i++)
		for (j = 0; j < n; j++)
			stateshare[i][j] ^= wshare[16 * round + i][j];
}

int inv_aes_share_subkeys(const byte in[16], byte out[16], byte *wshare[176], int n, int k)
{
	int i;
	int round = 10;
	byte *stateshare[16];

	for (i = 0; i < 16; i++)
	{
		stateshare[i] = (byte *)malloc(n * sizeof(byte));
		mask(stateshare[i], in[i], n, k);
	}

	if (DECIPHER_VERBOSE)
	{
		printf("\n**round %d\n", round);
		printMes("input", in);
	}
	addroundkey_share(stateshare, wshare, 10, n);

	if (DECIPHER_VERBOSE)
		unmaskPrintMes(stateshare, "add round key", n, k);

	for (round = 9; round >= 1; round--)
	{
		if (DECIPHER_VERBOSE)
			printf("\n**round %d\n", round);

		inv_shiftrows_share(stateshare, n);
		if (DECIPHER_VERBOSE)
			unmaskPrintMes(stateshare, "inv shift rows", n, k);

		inv_subbyte_share(stateshare, n);
		if (DECIPHER_VERBOSE)
			unmaskPrintMes(stateshare, "inv subbyte", n, k);

		addroundkey_share(stateshare, wshare, round, n);
		if (DECIPHER_VERBOSE)
			unmaskPrintMes(stateshare, "add round key", n, k);

		inv_mixcolumns_share(stateshare, n);
		if (DECIPHER_VERBOSE)
			unmaskPrintMes(stateshare, "inv mix columns", n, k);
	}

	inv_shiftrows_share(stateshare, n);
	inv_subbyte_share(stateshare, n);
	addroundkey_share(stateshare, wshare, 0, n);

	int status = 1; //determine if computaion is correct (not faulted)

	for (i = 0; i < 16; i++)
	{ //check for faults before proceeding

		int unmasked = unmask(stateshare[i], n, k);
		if (unmasked != NULL_BYTE) //TODO : leakage !!
			out[i] = (byte)unmasked;
		else
		{ //raise alarm we can choose to restart the computations
			out[i] = 0x00;
			status = 0;
		}
		memset(stateshare[i], 0, n);
		free(stateshare[i]);
	}
	return status; //everything went well
}

void inv_aes_share(const byte in[16], byte out[16], const byte key[16], int n, int k)
{

	int i;
	byte w[176];
	byte *wshare[176];

	keyexpansion((byte *)key, w);

	for (i = 0; i < 176; i++)
	{
		wshare[i] = (byte *)malloc(n * sizeof(byte));
		mask(wshare[i], w[i], n, k);
	}

	inv_aes_share_subkeys(in, out, wshare, n, k);

	for (i = 0; i < 176; i++)
	{
		memset(wshare[i], 0, n);
		free(wshare[i]);
	}
}

/**
 * @brief AES with shares and in 10 rounds. 
 * Chipher function in AES standard
 * 
 * @param in plaitext
 * @param out result
 * @param wshare shared sub-keys
 * @param n nb of shares
 * @param k duplicating parameter
 * @return the status 1 if OK and 0 if fault detected
 */
int aes_share_subkeys(const byte in[16], byte out[16], byte *wshare[176], int n, int k)
{
	int i;
	int round = 0;
	byte *stateshare[16];

	for (i = 0; i < 16; i++)
	{
		stateshare[i] = (byte *)malloc(n * sizeof(byte));
		mask(stateshare[i], in[i], n, k);
	}
	if (CIPHER_VERBOSE)
	{
		printf("\n**round %d\n", round);
		printMes("input", in);
	}
	addroundkey_share(stateshare, wshare, round, n);

	if (CIPHER_VERBOSE)
		unmaskPrintMes(stateshare, "add round key", n, k);
	for (round = 1; round < 10; round++)
	{
		if (CIPHER_VERBOSE)
			printf("\n**round %d\n", round);
		subbyte_share(stateshare, n, k);
		if (CIPHER_VERBOSE)
			unmaskPrintMes(stateshare, "subbyte", n, k);
		shiftrows_share(stateshare, n);
		if (CIPHER_VERBOSE)
			unmaskPrintMes(stateshare, "shift rows", n, k);
		mixcolumns_share(stateshare, n);
		if (CIPHER_VERBOSE)
			unmaskPrintMes(stateshare, "mix columns", n, k);
		addroundkey_share(stateshare, wshare, round, n);
		if (CIPHER_VERBOSE)
			unmaskPrintMes(stateshare, "add round key", n, k);
	}

	subbyte_share(stateshare, n, k);
	shiftrows_share(stateshare, n);
	addroundkey_share(stateshare, wshare, 10, n);

	int status = 1; //determine if computaion is correct (not faulted)

	for (i = 0; i < 16; i++)
	{ //check for faults before proceeding

		int unmasked = unmask(stateshare[i], n, k);
		if (unmasked != NULL_BYTE) //TODO : leakage !!
			out[i] = (byte)unmasked;
		else
		{ //raise alarm we can choose to restart the computations
			out[i] = 0x00;
			status = 0;
		}
		memset(stateshare[i], 0, n);
		free(stateshare[i]);
	}
	return status; //everything went well
}

//key expansion and then delegating to aes_share_subkeys()
//put in status 0 if faults detected and 1 ovewise
//return the time without the key expansion
int run_aes_share(const byte in[16], byte out[16], const byte key[16], int n, int k, int *status, int nt)
{
	int i;
	byte w[176]; //Nb*(Nr + 1)
	byte *wshare[176];
	clock_t start, end;

	keyexpansion((byte *)key, w); //TODO key schedule not protected

	for (i = 0; i < 176; i++)
	{
		wshare[i] = (byte *)malloc(n * sizeof(byte)); //shares of the sub-keys
		mask(wshare[i], w[i], n, k);
	}

	start = clock();
	for (i = 0; i < nt; i++)
	{
		*status = aes_share_subkeys(in, out, wshare, n, k);
		// if (*status == 0)
		// {
		// 	printf("Faults detected => need to restart over\n");
		// }
		//else
		//	printf("Everything went well\n");
	}
	end = clock();

	for (i = 0; i < 176; i++)
	{
		memset(wshare[i], 0, n);
		free(wshare[i]);
	}
	//return status;
	return (int)(end - start);
}

//IPM-FD protected key expansion and then 10 rounds of AES
//return status= 0 if faults detected and 1 ovewise
int AES_IPM_FD(const byte in[16], byte out[16], const byte key[16], int n, int k)
{
	int i, status;
	byte *wshare[176];

	for (i = 0; i < 176; i++)
	{
		wshare[i] = (byte *)malloc(n * sizeof(byte)); //shares of the sub-keys
		if (!wshare[i])
			error("malloc failed at wshare\n");
	}
	AES_IPM_FD_Key_Expansion(key, wshare, n, k);

	status = aes_share_subkeys(in, out, wshare, n, k);
	for (i = 0; i < 176; i++)
	{
		memset(wshare[i], 0, n);
		free(wshare[i]);
	}
	return status;
}

void printMes(const char *s, const byte *m)
{
	printf("%s = ", s);
	int i;
	for (i = 0; i < 16; i++)
		printf("%02x", m[i]);
	printf("\n");
}

void unmaskPrintMes(byte **stateshare, const char *s, int n, int k)
{
	int i;
	byte out[16];
	for (i = 0; i < 16; i++)
	{
		out[i] = unmask(stateshare[i], n, k);
	}
	printMes(s, out);
}

/**
 * @brief protected AES-128 key expansion
 * @param key the 128-bit key
 * @param wshare the subkeys = result of this procedure
 * @param n number of shares
 * @pama k duplicating parameter
 */
void AES_IPM_FD_Key_Expansion(const byte *key, byte *wshare[176], int n, int k)
{
	int i, j;
	byte temp[4][n];
	byte rcon[10] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36};
	//byte rcon[10];
	//setrcon(rcon);

	for (i = 0; i < 16; i++)
		mask(wshare[i], key[i], n, k);

	byte w_temp[n];
	for (i = 16; i < 176; i += 4)
	{
		for (j = 0; j < 4; j++)
			memcpy(temp[j], wshare[i - 4 + j], n);

		if ((i % 16) == 0)
		{
			sbox_share(w_temp, wshare[i - 3], n, k);
			IPConstAdd(temp[0], w_temp, rcon[i / 16 - 1], n, k);

			sbox_share(temp[1], wshare[i - 2], n, k);
			sbox_share(temp[2], wshare[i - 1], n, k);
			sbox_share(temp[3], wshare[i - 4], n, k);
		}

		for (j = 0; j < 4; j++)
			IPM_FD_Add(wshare[i + j], wshare[i + j - 16], temp[j], n);
	}
}
