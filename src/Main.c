#include <stdio.h>
#include <time.h>
#include <string.h>

#include "AES_IPM.h"
#include "IPM.h"
#include "aes.h"
#include "common.h"

#define TEST 1

void test()
{
    int n = 5,
            k = 2; //Nb of shares; duplicating parameter
    IPM_FD_Setup(n, k);
    //test mask | unmask
    byte X = 33;//random_byte();

    printf("X = %02x\n", X);

    byte Z[n];
    mask(Z, X, n, k);
    print(Z, "masked data ", n);

    X = unmask(Z, n, k);
    printf("Unmasked = %02x\n\n", X);

//    printf("Expected Squared data: b68c9dc2\n");

    printf("---------------- Test IPM_FD_Mult --------------\n");
    byte Y = 55; //random_byte();
    printf("Y = %02x\n", Y);
    byte Q[n];
    mask(Q, Y, n, k);
    print(Q, "(Sharing) masked data Q ", n);
    int u = unmask(Q, n, k);
    printf("Unmasked = %02x\n\n", u);

    /////
    printf("X*Y = %d\n", GF256_Mult(X, Y));
    byte T[n];

    IPM_FD_Mult(T, Z, Q, n, k);

    print(T, "Z*Q ", n);
    u = unmask(T, n, k);
    printf("Unmasked = %d\n\n", u);

    printf("---------------- Test square---------------- \n");
    mask(Z, X, n, k);
    IPM_FD_Square(Z,Z,n,k);
    print(Z, "Squared data ", n);

}
int main()
{
    if (TEST)
        test();
    else
    {
        printf("-------AES CIPHER-------\n");
        // byte key[16] = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
        // byte in[16] = {0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34};
        // byte outex[16] = {0x39, 0x25, 0x84, 0x1d, 0x02, 0xdc, 0x09, 0xfb, 0xdc, 0x11, 0x85, 0x97, 0x19, 0x6a, 0x0b, 0x32};

        byte key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f};
        byte in[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff};
        byte outex[16] = {0x69, 0xc4, 0xe0, 0xd8, 0x6a, 0x7b, 0x04, 0x30, 0xd8, 0xcd, 0xb7, 0x80, 0x70, 0xb4, 0xc5, 0x5a};

        printMes("in ", in);
        printMes("key ", key);
        byte out[16];
        memset(out, 0, 16);

        int n = 3, k = 2, nbyte = 16; //Nb of shares; duplicating parameter

        // printf("\n---------IPM-FD with protected Key Expansion----------\n");
        // memset(out, 0, n);
        // IPM_FD_Setup(n, k);
        // AES_IPM_FD(in, out, key, n, k);
        // check_ciphertext(out, outex, nbyte);
        // printMes("\nout", out);
        // freeMemory(n, k);
        // //----------------------------------

        //----------------------------------
        int dt, status = -1, nTimes = 1000, base = 0;
        printf("Without countermeasure:\n");
        dt = run_aes(in, out, key, nTimes);
        base = dt;
        check_ciphertext(out, outex, nbyte);
        report_time(dt, nTimes, base, 0); //randcount=0
        printf("-------------------------\n");

        for (k = 1; k <= 2; k++)
        {
            for (n = 2; n <= 4; n++)
            {
                if (n <= k)
                    continue;
                memset(out, 0, nbyte);
                IPM_FD_Setup(n, k);
                printf("\nWith extended IPM-FD n = %d, k = %d :\n", n, k);

                dt = run_aes_share(in, out, key, n, k, &status, nTimes);
                if (status == 0)
                {
                    printf("Faults detected => need to restart over\n");
                }
                //else
                //	printf("Everything went well\n");

                report_time(dt, nTimes, base, get_randcount());
                check_ciphertext(out, outex, nbyte);
                freeMemory(n, k);
                printf("-------------------------\n");
            }
        }
        printMes("\nout", out);
    }
    return 0;
}
