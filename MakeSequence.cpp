#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#define length 1000000
#define samplesize 1000
#define BLOCK_SIZE 500
#define BUF_SIZE (BLOCK_SIZE / 8 + 1)
typedef unsigned char byte;
char s[length * samplesize] = {0};
byte S[length * samplesize / 8] = {0};

//This Module is to adjust a block of sequence with the prescribed Linear complexity Profile

class State {
public:
    State() {}
    ~State() {}
    void Init()
    {
        N = 0;
        m = -1;
        L = 0;
        memset(B, 0, sizeof(B));
        memset(C, 0, sizeof(C));
        B[0] = C[0] = 1;
    }
    void OneRound()
    {
        char d = q[N];
        for (int i = 1; i <= L; i++)d ^= C[i] & q[N - i];
        if (d) {
            if (N + 1 - L > L) {
                char Temp[BLOCK_SIZE + 1];
                for (int i = 1; i <= L; i++)Temp[i] = C[i];
                for (int i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
                for (int i = 1; i <= L; i++)B[i] = Temp[i];
                L = N + 1 - L;
                m = N;
            } else {
                for (int i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
            }
        }
        N++;
    }
    /*pull down the point that is above the line y=0.5x*/
    void PullDownOneRound()
    {
        char d = q[N];
        for (int i = 1; i <= L; i++)d ^= C[i] & q[N - i];

        q[N] ^= d;

        OneRound();
    }
    /*push up the point that is below the line y=0.5x*/
    void PushUpOneRound()
    {
        char d = q[N];
        for (int i = 1; i <= L; i++)d ^= C[i] & q[N - i];

        q[N] ^= ~d;

        OneRound();
    }
    /*whether the point is on the line y=0.5x*/
    bool IsInSlash()
    {
        return 2 * L == N;
    }

    void ShowCurrentLCP()
    {
        Init();
        for (int i = 0; i < BLOCK_SIZE; i++) {
            OneRound();
            printf("(%3d, %3d)%c", N, L, (i + 1) % 10 == 0 ? '\n' : ' ');
        }
    }
public:
    int N;
    int m;
    int L;
    char B[BLOCK_SIZE + 1];
    char C[BLOCK_SIZE + 1];
    char q[BLOCK_SIZE];
};


/*
 * [n1] is the number of PCT T1s, [n2] is the number of PCT T2s, 
 * mm is the horizontal length of the big PCT,  [gap] is the tail 
 * length.
 * fp: file descriptor of the source sequence
 */
 // Generate a sequence with each block having the prescribed linear complexity profile 
void MakeSequence(int n1, int n2, int mm, int gap, FILE *fp) {
    State curState;
    State backState;
    int frontlen = BLOCK_SIZE - n1 * 2 - n2 * 4 - 2 * mm - gap;
    byte buf[BUF_SIZE];
    int tmpLen = 0;
    printf("=========================================================================\n");
    printf("n1 = %d\n", n1);
    printf("n2 = %d\n", n2);
    printf("mm = %d\n", mm);
    printf("frontlen = %d\n", frontlen);
    printf("=========================================================================\n");

    for (int loop = 0; loop < length * samplesize / BLOCK_SIZE; loop++) {
        fread(buf, 1, BUF_SIZE, fp);
        for (int i = 0; i < BLOCK_SIZE / 8; i++) {
            byte tmp = buf[i];
            for (int j = 0; j < 8 && 8 * i + j < BLOCK_SIZE; j++) {
                curState.q[8 * i + j] = (tmp & 0x80) == 0 ? 0 : 1;
                tmp <<= 1;
            }
        }
        
        /*initialize the Berlekamp-Massey algorithm*/
        curState.Init();
        backState = curState;
        for (int i = 0; i < frontlen; i++) {
            curState.OneRound();
            if (curState.IsInSlash()) {
                backState = curState;
            }
        }

        curState = backState;


        tmpLen = (frontlen - curState.N) / 2;
        if (tmpLen > 0) {
            for (int i = 0; i < tmpLen - 1; i++) {
                curState.OneRound();
            }
            curState.PushUpOneRound();
            for (int i = 0; i < tmpLen; i++) {
                curState.OneRound();
            }
        }

        for (int i = 0; i < n2; i++) {
            curState.PullDownOneRound();
            curState.PushUpOneRound();
            curState.OneRound();
            curState.OneRound();
        }
        
        for (int i = 0; i < n1; i++) {
            curState.PushUpOneRound();
            curState.OneRound();
        }
        
        for (int i = 0; i < mm - 1; i++) {
            curState.PullDownOneRound();
        }
        curState.PushUpOneRound();

        // store the result block to sequence
        for (int i = 0; i < BLOCK_SIZE; i++) {
            s[BLOCK_SIZE * loop + i] = curState.q[i];
        }
    }
}

int main(int argc, char **argv) {
    FILE *fp_in = fopen("rdseed.data", "rb");
    if (fp_in == NULL) {
        printf("fopen file rdseed.data failed\n");
        return 1;
    }
    MakeSequence(111, 46, 22, 10, fp_in);  // Generate a block of sequence with the prescribed linear complexity profile, where the parameters are the number of PCTs like T1, T2 ...
    fclose(fp_in);

    FILE *fp_out  = fopen("sample.data", "wb");
    if (fp_out == NULL) {
        printf("fopen file sample.data for writing failed\n");
        return 2;
    }

    // transfer the bit sequence to byte sequence
    for (int i = 0; i < length * samplesize / 8; i++) {
        for (int j = 0; j < 8; j++) {
            S[i] <<= 1;
            S[i] ^= s[8 * i + j]; 
        }
    }

    fwrite(S, 1, length * samplesize / 8, fp_out);
    fclose(fp_out);

    printf("success\n");
    return 0;
}
