#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;
const int MaxN = 300;


/*
所有二进制使用大头表示法
即S_0表示S的最高位，S_n表示S的最低位
*/


//明文向量
int mr[MaxN * MaxN];
int mg[MaxN * MaxN];
int mb[MaxN * MaxN];
//密文向量
int cr[MaxN * MaxN];
int cg[MaxN * MaxN];
int cb[MaxN * MaxN];
//攻击得到密文
int ar[MaxN * MaxN];
int ag[MaxN * MaxN];
int ab[MaxN * MaxN];
//控制参数
double a = 2.0;
double b = 1.99;
//像素和
int m, n;
double SumR = 0, SumG = 0, SumB = 0;

double X1[2 * MaxN * MaxN], Y1[2 * MaxN * MaxN], Z1[2 * MaxN * MaxN], G1[2 * MaxN * MaxN];
double X2[2 * MaxN * MaxN], Y2[2 * MaxN * MaxN], Z2[2 * MaxN * MaxN], G2[2 * MaxN * MaxN];

//置换向量
struct sortseq
{
    int id;
    double val;
}
s23[MaxN * MaxN], s22[MaxN * MaxN], s21[MaxN * MaxN], s20[MaxN * MaxN],
s43[MaxN * MaxN], s42[MaxN * MaxN], s41[MaxN * MaxN], s40[MaxN * MaxN],
s13[MaxN * MaxN], s12[MaxN * MaxN], s11[MaxN * MaxN], s10[MaxN * MaxN],
s33[MaxN * MaxN], s32[MaxN * MaxN], s31[MaxN * MaxN], s30[MaxN * MaxN];
//逆置换向量
int invs43[MaxN * MaxN], invs42[MaxN * MaxN], invs41[MaxN * MaxN], invs40[MaxN * MaxN],
invs23[MaxN * MaxN], invs22[MaxN * MaxN], invs21[MaxN * MaxN], invs20[MaxN * MaxN],
invs13[MaxN * MaxN], invs12[MaxN * MaxN], invs11[MaxN * MaxN], invs10[MaxN * MaxN];
//异或向量
int u1[MaxN * MaxN], v1[MaxN * MaxN], w1[MaxN * MaxN];
int u2[MaxN * MaxN], v2[MaxN * MaxN], w2[MaxN * MaxN];
int u1bit[MaxN * MaxN][4];
int u2bit[MaxN * MaxN][4];

bool cmp(sortseq a0, sortseq b0)
{
    if (a0.val == b0.val) return a0.id < b0.id;
    return a0.val < b0.val;
}

//2D-LCLM
void LCLM(double x0, double y0, double z0, double* X, double* Y, double* Z, double* V)
{
    int k = 0;
    for (int i = 1; i <= 250; i++)
    {
        double x = b * x0 * (1 - z0);
        double y = b * y0 * (1 - z0);
        double z = a * x0 * x0 + y0 * y0;
        x0 = x; y0 = y; z0 = z;
    }
    for (int i = 0; i < 2 * m * n; i++)
    {
        double x = b * x0 * (1 - z0);
        double y = b * y0 * (1 - z0);
        double z = a * x0 * x0 + y0 * y0;
        x0 = x; y0 = y; z0 = z;
        X[i] = x0;
        Y[i] = y0;
        Z[i] = z0;
        V[i] = (x0 + y0 + z0) / 3.0;
    }
}

void init()
{

    LCLM(0.2 + SumR * 1e-9, 0.4 + SumG * 1e-9, 0.1 + SumB * 1e-9, X1, Y1, Z1, G1);
    LCLM(0.3 + SumR * 1e-9, 0.5 + SumG * 1e-9, 0.2 + SumB * 1e-9, X2, Y2, Z2, G2);
    for (int i = 0; i < n * m; i++)
    {
        u1[i] = (int)floor(fabs(X1[i]) * 1e15) % 16;
        v1[i] = (int)floor(((double)(Y1[i] * 1e3 - floor(Y1[i] * 1e3))) * 1e3) % 256;
        w1[i] = (int)floor(((double)(Z1[i] * 1e3 - floor(Z1[i] * 1e3))) * 1e3) % 256;

        u2[i] = (int)floor(fabs(X2[i]) * 1e15) % 16;
        v2[i] = (int)floor(((double)(Y2[i] * 1e3 - floor(Y2[i] * 1e3))) * 1e3) % 256;
        w2[i] = (int)floor(((double)(Z2[i] * 1e3 - floor(Z2[i] * 1e3))) * 1e3) % 256;

        for (int j = 3; j >= 0; j--)
        {
            u1bit[i][j] = u1[i] % 2;
            u2bit[i][j] = u2[i] % 2;
            u1[i] /= 2;
            u2[i] /= 2;
        }
    }
    for (int i = 0; i < n * m; i++)
    {
        s23[i].id = s13[i].id = i; s23[i].val = G1[i]; s13[i].val = G1[i + n * m];
        s22[i].id = s12[i].id = i; s22[i].val = Z1[i]; s12[i].val = Z1[i + n * m];
        s21[i].id = s11[i].id = i; s21[i].val = Y1[i]; s11[i].val = Y1[i + n * m];
        s20[i].id = s10[i].id = i; s20[i].val = X1[i]; s10[i].val = X1[i + n * m];

        s43[i].id = s33[i].id = i; s43[i].val = G2[i]; s33[i].val = G2[i + n * m];
        s42[i].id = s32[i].id = i; s42[i].val = Z2[i]; s32[i].val = Z2[i + n * m];
        s41[i].id = s31[i].id = i; s41[i].val = Y2[i]; s31[i].val = Y2[i + n * m];
        s40[i].id = s30[i].id = i; s40[i].val = X2[i]; s30[i].val = X2[i + n * m];
    }
    sort(s23, s23 + m * n, cmp); sort(s13, s13 + m * n, cmp);
    sort(s22, s22 + m * n, cmp); sort(s12, s12 + m * n, cmp);
    sort(s21, s21 + m * n, cmp); sort(s11, s11 + m * n, cmp);
    sort(s20, s20 + m * n, cmp); sort(s10, s10 + m * n, cmp);
    sort(s43, s43 + m * n, cmp); sort(s33, s33 + m * n, cmp);
    sort(s42, s42 + m * n, cmp); sort(s32, s32 + m * n, cmp);
    sort(s41, s41 + m * n, cmp); sort(s31, s31 + m * n, cmp);
    sort(s40, s40 + m * n, cmp); sort(s30, s30 + m * n, cmp);
}

//加密过程
int Pm[MaxN * MaxN];
int Cm[MaxN * MaxN];
int Pmbit[MaxN * MaxN][8];
int Cmbit[MaxN * MaxN][8];
int Tmbit[MaxN * MaxN][8];
void encryption(int* P, int* C)
{
    /*
    cout << "P*" << endl;
    for (int i = 0; i < n * m; i++)
    {
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << (P[i] + v1[i]) % 256 << endl;
        }
    }
    */
    for (int i = 0; i < n * m; i++)
    {
        Pm[i] = w1[i] ^ ((P[i] + v1[i]) % 256);
        for (int j = 7; j >= 0; j--)
        {
            Pmbit[i][j] = Pm[i] % 2;
            Pm[i] /= 2;
        }
    }
    //S1置换
    for (int i = 0; i < n * m; i++)
    {
        Cmbit[i][4] = Pmbit[s13[i].id][4];
        Cmbit[i][5] = Pmbit[s12[i].id][5];
        Cmbit[i][6] = Pmbit[s11[i].id][6];
        Cmbit[i][7] = Pmbit[s10[i].id][7];
    }
    //S2置换
    for (int i = 0; i < n * m; i++)
    {
        Cmbit[i][0] = Pmbit[s23[i].id][0];
        Cmbit[i][1] = Pmbit[s22[i].id][1];
        Cmbit[i][2] = Pmbit[s21[i].id][2];
        Cmbit[i][3] = Pmbit[s20[i].id][3];
    }

    for (int i = 0; i < n * m; i++)
    {
        Pmbit[i][4] = u1bit[i][0] ^ Cmbit[i][0] ^ Cmbit[i][4];
        Pmbit[i][5] = u1bit[i][1] ^ Cmbit[i][1] ^ Cmbit[i][5];
        Pmbit[i][6] = u1bit[i][2] ^ Cmbit[i][2] ^ Cmbit[i][6];
        Pmbit[i][7] = u1bit[i][3] ^ Cmbit[i][3] ^ Cmbit[i][7];
    }

    //S3置换
    for (int i = 0; i < n * m; i++)
    {
        Cmbit[i][4] = Pmbit[s33[i].id][4];
        Cmbit[i][5] = Pmbit[s32[i].id][5];
        Cmbit[i][6] = Pmbit[s31[i].id][6];
        Cmbit[i][7] = Pmbit[s30[i].id][7];
    }
    //S4置换
    for (int i = 0; i < n * m; i++)
    {
        Cmbit[i][0] = Pmbit[s43[i].id][0];
        Cmbit[i][1] = Pmbit[s42[i].id][1];
        Cmbit[i][2] = Pmbit[s41[i].id][2];
        Cmbit[i][3] = Pmbit[s40[i].id][3];
    }
    for (int i = 0; i < n * m; i++)
    {
        Pmbit[i][0] = u2bit[i][0] ^ Cmbit[i][0] ^ Cmbit[i][4];
        Pmbit[i][1] = u2bit[i][1] ^ Cmbit[i][1] ^ Cmbit[i][5];
        Pmbit[i][2] = u2bit[i][2] ^ Cmbit[i][2] ^ Cmbit[i][6];
        Pmbit[i][3] = u2bit[i][3] ^ Cmbit[i][3] ^ Cmbit[i][7];
    }

    for (int i = 0; i < n * m; i++)
    {
        Pm[i] = Pmbit[i][0] * 128 + Pmbit[i][1] * 64 + Pmbit[i][2] * 32 + Pmbit[i][3] * 16
            + Pmbit[i][4] * 8 + Pmbit[i][5] * 4 + Pmbit[i][6] * 2 + Pmbit[i][7] * 1;
        C[i] = w2[i] ^ ((Pm[i] + v2[i]) % 256);
    }
    /*
    cout << "Ps" << endl;
    for (int i = 0; i < n * m; i++)
    {
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << Pm[i] << endl;
        }
    }

    cout << "C" << endl;
    for (int i = 0; i < n * m; i++)
    {
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << C[i] << endl;
        }
    }
    */
}


//解密过程
void decryption(int* C, int* P)
{
    for (int i = 0; i < n * m; i++)
    {
        Cm[i] = ((C[i] ^ w2[i]) + 256 - v2[i]) % 256;
        //Cm[i] = C[i];
        for (int j = 7; j >= 0; j--)
        {
            Cmbit[i][j] = Cm[i] % 2;
            Cm[i] /= 2;
        }
    }

    //计算逆置换序列
    for (int i = 0; i < n * m; i++)
    {
        invs43[s43[i].id] = i;
        invs42[s42[i].id] = i;
        invs41[s41[i].id] = i;
        invs40[s40[i].id] = i;
        invs13[s13[i].id] = i;
        invs12[s12[i].id] = i;
        invs11[s11[i].id] = i;
        invs10[s10[i].id] = i;
    }

    //S3置换
    for (int i = 0; i < n * m; i++)
    {
        Pmbit[i][4] = Cmbit[s33[i].id][4];
        Pmbit[i][5] = Cmbit[s32[i].id][5];
        Pmbit[i][6] = Cmbit[s31[i].id][6];
        Pmbit[i][7] = Cmbit[s30[i].id][7];
    }

    for (int i = 0; i < n * m; i++)
    {
        Tmbit[i][0] = (u2bit[i][0] ^ Pmbit[i][4] ^ Cmbit[i][0]);
        Tmbit[i][1] = (u2bit[i][1] ^ Pmbit[i][5] ^ Cmbit[i][1]);
        Tmbit[i][2] = (u2bit[i][2] ^ Pmbit[i][6] ^ Cmbit[i][2]);
        Tmbit[i][3] = (u2bit[i][3] ^ Pmbit[i][7] ^ Cmbit[i][3]);
    }

    //S4逆置换
    for (int i = 0; i < n * m; i++)
    {
        Pmbit[i][0] = Tmbit[invs43[i]][0];
        Pmbit[i][1] = Tmbit[invs42[i]][1];
        Pmbit[i][2] = Tmbit[invs41[i]][2];
        Pmbit[i][3] = Tmbit[invs40[i]][3];
    }

    //S2置换
    for (int i = 0; i < n * m; i++)
    {
        Cmbit[i][0] = Pmbit[s23[i].id][0];
        Cmbit[i][1] = Pmbit[s22[i].id][1];
        Cmbit[i][2] = Pmbit[s21[i].id][2];
        Cmbit[i][3] = Pmbit[s20[i].id][3];
    }

    for (int i = 0; i < n * m; i++)
    {
        Tmbit[i][4] = (u1bit[i][0] ^ Cmbit[i][0] ^ Cmbit[i][4]);
        Tmbit[i][5] = (u1bit[i][1] ^ Cmbit[i][1] ^ Cmbit[i][5]);
        Tmbit[i][6] = (u1bit[i][2] ^ Cmbit[i][2] ^ Cmbit[i][6]);
        Tmbit[i][7] = (u1bit[i][3] ^ Cmbit[i][3] ^ Cmbit[i][7]);
    }

    //S1逆置换
    for (int i = 0; i < n * m; i++)
    {
        Pmbit[i][4] = Tmbit[invs13[i]][4];
        Pmbit[i][5] = Tmbit[invs12[i]][5];
        Pmbit[i][6] = Tmbit[invs11[i]][6];
        Pmbit[i][7] = Tmbit[invs10[i]][7];
    }


    for (int i = 0; i < n * m; i++)
    {
        Pm[i] = Pmbit[i][0] * 128 + Pmbit[i][1] * 64 + Pmbit[i][2] * 32 + Pmbit[i][3] * 16
            + Pmbit[i][4] * 8 + Pmbit[i][5] * 4 + Pmbit[i][6] * 2 + Pmbit[i][7] * 1;
        P[i] = ((Pm[i] ^ w1[i]) + 256 - v1[i]) % 256;
        //P[i] = Pm[i];
    }
}


//攻击(256*256为例)
int PA0[MaxN * MaxN];
int CA0[MaxN * MaxN];
int PA[16][MaxN * MaxN];
int CA[16][MaxN * MaxN];
int recordp[MaxN * MaxN];

int as23[MaxN * MaxN], as22[MaxN * MaxN], as21[MaxN * MaxN], as20[MaxN * MaxN],
as43[MaxN * MaxN], as42[MaxN * MaxN], as41[MaxN * MaxN], as40[MaxN * MaxN],
as13[MaxN * MaxN], as12[MaxN * MaxN], as11[MaxN * MaxN], as10[MaxN * MaxN],
as33[MaxN * MaxN], as32[MaxN * MaxN], as31[MaxN * MaxN], as30[MaxN * MaxN];

void CPAinit()
{
    for (int i = 0; i < n * m; i++) PA0[i] = 0;
    encryption(PA0, CA0);
}

void atts2()
{
    for (int k = 16; k <= 128; k *= 2)
    {
        if (k == 32) continue;
        for (int i = 0; i < n * m; i++)
        {
            int ts = i;
            for (int j = 15; j >= 0; j--)
            {
                PA[j][i] = (ts % 2) * k;
                ts /= 2;
            }
        }
        for (int j = 0; j < 16; j++)
        {
            encryption(PA[j], CA[j]);
            for (int i = 0; i < n * m; i++)
            {
                CA[j][i] = (CA[j][i] ^ CA0[i]);
                CA[j][i] /= (k / 16);
                CA[j][i] %= 2;
            }
        }
        for (int i = 0; i < n * m; i++)
        {
            int ts = 1;
            recordp[i] = 0;
            for (int j = 15; j >= 0; j--)
            {
                recordp[i] += (CA[j][i] * ts);
                ts *= 2;
            }
        }
        if (k == 16)
        {
            for (int i = 0; i < n * m; i++)
            {
                as20[i] = recordp[i];
                as21[i] = recordp[i];
            }
        }
        else if (k == 64)
        {
            for (int i = 0; i < n * m; i++)
            {
                as22[i] = recordp[i];
            }
        }
        else if (k == 128)
        {
            for (int i = 0; i < n * m; i++)
            {
                as23[i] = recordp[i];
            }
        }
    }
    int fg = 0;
    for (int i = 0; i < n * m; i++)
    {
        if (as23[i] != s23[i].id) { fg = 1; break; }
        if (as22[i] != s22[i].id) { fg = 1; break; }
        if (as21[i] != s21[i].id) { fg = 1; break; }
        if (as20[i] != s20[i].id) { fg = 1; break; }
    }
    if (fg == 0) cout << "S2 determined." << endl;
    else cout << "S2 can not be determined." << endl;
    /*
    cout << "T20" << endl;
    for (int i = 0; i < n * m; i++)
    {
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << as20[i] << endl;
        }
    }
    */
}


int avl[MaxN * MaxN];
int avltp[MaxN * MaxN];
int invas20[MaxN * MaxN];
int invas21[MaxN * MaxN];
int invas22[MaxN * MaxN];
int invas23[MaxN * MaxN];
int vis[MaxN * MaxN];
void attvl()
{
    for (int i = 0; i < n * m; i++)
    {
        vis[i] = 0;
        for (int k = 1; k <= 15; k++)
            PA[k][i] = k;
    }
    for (int k = 1; k <= 15; k++)
    {
        encryption(PA[k], CA[k]);
        for (int i = 0; i < n * m; i++)
        {
            CA[k][i] = (CA[k][i] ^ CA0[i]);
            CA[k][i] %= 2;
        }
    }
    for (int k = 1; k <= 15; k++)
    {
        for (int i = 0; i < n * m; i++)
        {
            if (CA[k][i] == ((k % 2) ^ 1) && vis[i] == 0)
            {
                avltp[i] = 16 - k;
                vis[i] = 1;
            }
        }
    }
    for (int i = 0; i < n * m; i++)
    {
        invas20[as20[i]] = i;
        invas21[as21[i]] = i;
        invas22[as22[i]] = i;
        invas23[as23[i]] = i;
    }
    for (int i = 0; i < n * m; i++) avl[i] = avltp[invas20[i]];

    int fg = 0;
    for (int i = 0; i < n * m; i++)
    {
        if (v1[i] % 16 != avl[i]) { fg = 1; break; }
    }
    if (fg == 0) cout << "VL determined." << endl;
    else cout << "VL can not be determined." << endl;

}


int PA1[16][MaxN * MaxN];
int CA1[16][MaxN * MaxN];
void atts1()
{
    for (int k = 1; k <= 8; k *= 2)
    {
        if (k == 2) continue;
        for (int i = 0; i < n * m; i++)
        {
            int ts = i;
            for (int j = 15; j >= 0; j--)
            {
                PA[j][i] = (ts % 2) * k;
                PA1[j][i] = (((PA[j][i] + avl[i]) / 16) * 16);
                ts /= 2;
            }
        }
        for (int j = 0; j < 16; j++)
        {
            encryption(PA[j], CA[j]);
            encryption(PA1[j], CA1[j]);
            for (int i = 0; i < n * m; i++)
            {
                CA[j][i] = (CA[j][i] ^ CA1[j][i]);
                CA[j][i] /= k;
                CA[j][i] %= 2;
            }
        }
        for (int i = 0; i < n * m; i++)
        {
            int ts = 1;
            recordp[i] = 0;
            for (int j = 15; j >= 0; j--)
            {
                recordp[i] += (CA[j][i] * ts);
                ts *= 2;
            }
        }
        if (k == 1)
        {
            for (int i = 0; i < n * m; i++)
            {
                as10[i] = recordp[i];
                as11[i] = recordp[i];
            }
        }
        else if (k == 4)
        {
            for (int i = 0; i < n * m; i++)
            {
                as12[i] = recordp[i];
            }
        }
        else if (k == 8)
        {
            for (int i = 0; i < n * m; i++)
            {
                as13[i] = recordp[i];
            }
        }
    }
    int fg = 0;
    for (int i = 0; i < n * m; i++)
    {
        if (as13[i] != s13[i].id) { fg = 1; break; }
        if (as12[i] != s12[i].id) { fg = 1; break; }
        if (as11[i] != s11[i].id) { fg = 1; break; }
        if (as10[i] != s10[i].id) { fg = 1; break; }
    }
    if (fg == 0) cout << "S1 determined." << endl;
    else cout << "S1 can not be determined." << endl;
    /*
    cout << "T11" << endl;
    for (int i = 0; i < n * m; i++)
    {
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << as11[i] << endl;
        }
    }
    */

}

int avh[MaxN * MaxN];
int av[MaxN * MaxN];
int avhbit[MaxN * MaxN][4];
int avhtpbit[MaxN * MaxN][4];
int L0s[MaxN * MaxN];
int aL1[MaxN * MaxN];
int aH1[MaxN * MaxN];
int ar1[MaxN * MaxN];
void attvh()
{
    for (int i = 0; i < n * m; i++)
        L0s[i] = avl[i];
    for (int k = 1; k <= 3; k++)
    {
        for (int i = 0; i < n * m; i++)
        {
            aL1[i] = ((L0s[i] ^ (1 << (k - 1))) + 16 - avl[i]) % 16;
            ar1[i] = (aL1[i] + avl[i]) / 16;
            aH1[i] = ((1 << (k - 1)) - ar1[i]) % 16;
            PA[k][i] = aH1[i] * 16 + aL1[i];
        }
    }

    for (int k = 1; k <= 3; k++)
    {
        encryption(PA[k], CA[k]);
        for (int i = 0; i < n * m; i++)
        {
            CA[k][i] = (CA[k][i] ^ CA0[i]);
            avhtpbit[i][4 - k] = (CA[k][i] / (1 << k)) % 2;
        }
    }

    for (int i = 0; i < n * m; i++)
    {
        avhbit[i][3] = avhtpbit[invas21[i]][3];
        avhbit[i][2] = avhtpbit[invas22[i]][2];
        avhbit[i][1] = avhtpbit[invas23[i]][1];
        avh[i] = avhbit[i][1] * 4 + avhbit[i][2] * 2 + avhbit[i][3];
        av[i] = avh[i] * 16 + avl[i];
    }

    int fg = 0;
    for (int i = 0; i < n * m; i++)
    {
        if ((v1[i] / 16) % 8 != avh[i]) { fg = 1; break; }
    }
    if (fg == 0) cout << "VH determined." << endl;
    else cout << "VH can not be determined." << endl;
    /*
    cout << "AV" << endl;
    for (int i = 0; i < n * m; i++)
    {
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << av[i] << endl;
        }
    }
    */
}

int invas10[MaxN * MaxN];
int invas11[MaxN * MaxN];
int invas12[MaxN * MaxN];
int invas13[MaxN * MaxN];

int H1s[16][MaxN * MaxN];
int H1sbit[16][MaxN * MaxN][4];
int L1s[16][MaxN * MaxN];
int L1sbit[16][MaxN * MaxN][4];
int swapbit[MaxN * MaxN][4];



void atts4()
{
    for (int i = 0; i < n * m; i++)
    {
        invas10[as10[i]] = i;
        invas11[as11[i]] = i;
        invas12[as12[i]] = i;
        invas13[as13[i]] = i;
    }
    for (int i = 0; i < n * m; i++)
    {
        PA0[i] = (256 - av[i]) % 256;
    }
    encryption(PA0, CA0);
    for (int k = 1; k <= 8; k *= 2)
    {
        if (k == 2) continue;
        for (int i = 0; i < n * m; i++)
        {
            int ts = i;
            for (int j = 15; j >= 0; j--)
            {
                H1s[j][i] = (ts % 2) * k;
                int hspt = H1s[j][i];

                for (int ji = 3; ji >= 0; ji--)
                {
                    H1sbit[j][i][ji] = hspt % 2;
                    hspt /= 2;
                }
                ts /= 2;
            }
        }
        for (int j = 0; j < 16; j++)
        {
            for (int i = 0; i < n * m; i++)
            {
                swapbit[i][0] = H1sbit[j][as23[i]][0];
                swapbit[i][1] = H1sbit[j][as22[i]][1];
                swapbit[i][2] = H1sbit[j][as21[i]][2];
                swapbit[i][3] = H1sbit[j][as20[i]][3];
            }
            for (int i = 0; i < n * m; i++)
            {
                L1sbit[j][i][0] = swapbit[invas13[i]][0];
                L1sbit[j][i][1] = swapbit[invas12[i]][1];
                L1sbit[j][i][2] = swapbit[invas11[i]][2];
                L1sbit[j][i][3] = swapbit[invas10[i]][3];
            }
            for (int i = 0; i < n * m; i++)
            {
                L1s[j][i] = L1sbit[j][i][0] * 8 + L1sbit[j][i][1] * 4 + L1sbit[j][i][2] * 2 + L1sbit[j][i][3];

                PA[j][i] = ((H1s[j][i] * 16 + L1s[j][i]) + 256 - av[i]) % 256;
            }
        }

        for (int j = 0; j < 16; j++)
        {
            encryption(PA[j], CA[j]);
            for (int i = 0; i < n * m; i++)
            {
                CA[j][i] = (CA[j][i] ^ CA0[i]);
                CA[j][i] /= (k * 16);
                CA[j][i] %= 2;
            }
        }
        for (int i = 0; i < n * m; i++)
        {
            int ts = 1;
            recordp[i] = 0;
            for (int j = 15; j >= 0; j--)
            {
                recordp[i] += (CA[j][i] * ts);
                ts *= 2;
            }
        }
        if (k == 1)
        {
            for (int i = 0; i < n * m; i++)
            {
                as40[i] = recordp[i];
                as41[i] = recordp[i];
            }
        }
        else if (k == 4)
        {
            for (int i = 0; i < n * m; i++)
            {
                as42[i] = recordp[i];
            }
        }
        else if (k == 8)
        {
            for (int i = 0; i < n * m; i++)
            {
                as43[i] = recordp[i];
            }
        }
    }
    int fg = 0;
    for (int i = 0; i < n * m; i++)
    {
        if (as43[i] != s43[i].id) { fg = 1; break; }
        if (as42[i] != s42[i].id) { fg = 1; break; }
        if (as41[i] != s41[i].id) { fg = 1; break; }
        if (as40[i] != s40[i].id) { fg = 1; break; }
    }
    if (fg == 0) cout << "S4 determined." << endl;
    else cout << "S4 can not be determined." << endl;
    /*
    cout << "T42" << endl;
    for (int i = 0; i < n * m; i++)
    {
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << as42[i] << endl;
        }
    }
    */

}

int phi[MaxN * MaxN][4];
int L1sr[16][MaxN * MaxN];
int L1srbit[16][MaxN * MaxN][4];

void atts3()
{
    for (int k = 0; k <= 3; k++)
    {
        if (k == 1) continue;
        int k2b = (1 << k);
        for (int i = 0; i < n * m; i++)
        {
            L1sr[k][i] = k2b;
            PA[k][i] = (L1sr[k][i] + 256 - av[i]);
        }
        encryption(PA[k], CA[k]);
        for (int i = 0; i < n * m; i++)
        {
            CA[k][i] ^= CA0[i];
            CA[k][i] /= (16 * k2b);
            CA[k][i] %= 2;
            phi[i][k] = CA[k][i] ^ 1;
        }
    }


    for (int k = 0; k <= 3; k++)
    {
        if (k == 1) continue;
        int k2b = (1 << k);
        for (int i = 0; i < n * m; i++)
        {
            int ts = i;
            for (int j = 15; j >= 0; j--)
            {
                L1sr[j][i] = (ts % 2) * k2b;
                int hspt = L1sr[j][i];

                for (int ji = 3; ji >= 0; ji--)
                {
                    L1srbit[j][i][ji] = hspt % 2;
                    hspt /= 2;
                }
                ts /= 2;
            }
        }

        for (int j = 0; j < 16; j++)
        {
            for (int i = 0; i < n * m; i++)
            {
                L1sbit[j][i][0] = L1srbit[j][invas13[i]][0];
                L1sbit[j][i][1] = L1srbit[j][invas12[i]][1];
                L1sbit[j][i][2] = L1srbit[j][invas11[i]][2];
                L1sbit[j][i][3] = L1srbit[j][invas10[i]][3];
            }
            for (int i = 0; i < n * m; i++)
            {
                L1s[j][i] = L1sbit[j][i][0] * 8 + L1sbit[j][i][1] * 4 + L1sbit[j][i][2] * 2 + L1sbit[j][i][3];

                PA[j][i] = (L1s[j][i] + 256 - av[i]) % 256;
            }
        }

        for (int j = 0; j < 16; j++)
        {
            encryption(PA[j], CA[j]);
            for (int i = 0; i < n * m; i++)
            {
                CA[j][i] = (CA[j][i] ^ CA0[i]);
                CA[j][i] /= (k2b * 16);
                CA[j][i] %= 2;

                if ((L1sr[j][i] / k2b) == 1)
                    CA[j][i] ^= phi[i][k];
            }
        }

        for (int i = 0; i < n * m; i++)
        {
            int ts = 1;
            recordp[i] = 0;
            for (int j = 15; j >= 0; j--)
            {
                recordp[i] += (CA[j][i] * ts);
                ts *= 2;
            }
        }
        if (k == 0)
        {
            for (int i = 0; i < n * m; i++)
            {
                as30[i] = recordp[i];
                as31[i] = recordp[i];
            }
        }
        else if (k == 2)
        {
            for (int i = 0; i < n * m; i++)
            {
                as32[i] = recordp[i];
            }
        }
        else if (k == 3)
        {
            for (int i = 0; i < n * m; i++)
            {
                as33[i] = recordp[i];
            }
        }

    }

    int fg = 0;
    for (int i = 0; i < n * m; i++)
    {
        if (as33[i] != s33[i].id) { fg = 1; break; }
        if (as32[i] != s32[i].id) { fg = 1; break; }
        if (as31[i] != s31[i].id) { fg = 1; break; }
        if (as30[i] != s30[i].id) { fg = 1; break; }
    }
    if (fg == 0) cout << "S3 determined." << endl;
    else cout << "S3 can not be determined." << endl;
    /*
    cout << "T33" << endl;
    for (int i = 0; i < n * m; i++)
    {
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << as33[i] << endl;
        }
    }
    */
}

int Ls[MaxN * MaxN];
int Lsr[MaxN * MaxN];
int Hs[MaxN * MaxN];
int Hsr[MaxN * MaxN];

int Lsbit[MaxN * MaxN][4];
int Lsrbit[MaxN * MaxN][4];
int Hsbit[MaxN * MaxN][4];
int Hsrbit[MaxN * MaxN][4];
int xormap[MaxN * MaxN][256];

int invas40[MaxN * MaxN];
int invas41[MaxN * MaxN];
int invas42[MaxN * MaxN];
int invas43[MaxN * MaxN];

void attxor()
{
    for (int i = 0; i < n * m; i++)
    {
        invas40[as40[i]] = i;
        invas41[as41[i]] = i;
        invas42[as42[i]] = i;
        invas43[as43[i]] = i;
    }
    for (int k = 0; k < 256; k++)
    {

        for (int i = 0; i < m * n; i++)
        {
            int ts = k;
            for (int j = 3; j >= 0; j--)
            {
                Lsrbit[i][j] = ts % 2;
                ts /= 2;
            }
            for (int j = 3; j >= 0; j--)
            {
                Hsrbit[i][j] = ts % 2;
                ts /= 2;
            }
        }
        for (int i = 0; i < m * n; i++)
        {
            for (int j = 3; j >= 0; j--)
            {
                Hsbit[i][j] = Hsrbit[i][j] ^ Lsrbit[i][j];
                Lsbit[i][j] = Hsbit[i][j] ^ Lsrbit[i][j];
            }
            PA0[i] = Hsbit[i][0] * 128 + Hsbit[i][1] * 64 + Hsbit[i][2] * 32 + Hsbit[i][3] * 16
                + Lsbit[i][0] * 8 + Lsbit[i][1] * 4 + Lsbit[i][2] * 2 + Lsbit[i][3] * 1;
            PA0[i] = (PA0[i] + 256 - av[i]) % 256;
        }
        encryption(PA0, CA0);
        for (int i = 0; i < m * n; i++)
        {
            xormap[i][CA0[i]] = k;
        }
    }
}

void CPA(int* C, int* P)
{
    for (int i = 0; i < m * n; i++)
    {
        int ts = xormap[i][C[i]];
        for (int j = 3; j >= 0; j--)
        {
            Lsrbit[i][j] = ts % 2;
            ts /= 2;
        }
        for (int j = 3; j >= 0; j--)
        {
            Hsrbit[i][j] = ts % 2;
            ts /= 2;
        }
    }
    for (int i = 0; i < m * n; i++)
    {
        swapbit[i][0] = Lsrbit[as33[i]][0] ^ Hsrbit[i][0];
        swapbit[i][1] = Lsrbit[as32[i]][1] ^ Hsrbit[i][1];
        swapbit[i][2] = Lsrbit[as31[i]][2] ^ Hsrbit[i][2];
        swapbit[i][3] = Lsrbit[as30[i]][3] ^ Hsrbit[i][3];
    }

    for (int i = 0; i < m * n; i++)
    {
        Hsbit[i][0] = swapbit[invas43[i]][0];
        Hsbit[i][1] = swapbit[invas42[i]][1];
        Hsbit[i][2] = swapbit[invas41[i]][2];
        Hsbit[i][3] = swapbit[invas40[i]][3];
    }

    for (int i = 0; i < m * n; i++)
    {
        swapbit[i][0] = Hsbit[as23[i]][0] ^ Lsrbit[i][0];
        swapbit[i][1] = Hsbit[as22[i]][1] ^ Lsrbit[i][1];
        swapbit[i][2] = Hsbit[as21[i]][2] ^ Lsrbit[i][2];
        swapbit[i][3] = Hsbit[as20[i]][3] ^ Lsrbit[i][3];
    }

    for (int i = 0; i < m * n; i++)
    {
        Lsbit[i][0] = swapbit[invas13[i]][0];
        Lsbit[i][1] = swapbit[invas12[i]][1];
        Lsbit[i][2] = swapbit[invas11[i]][2];
        Lsbit[i][3] = swapbit[invas10[i]][3];
    }
    
    //cout << "P* P" << endl;
    for (int i = 0; i < m * n; i++)
    {
        int resimg = Hsbit[i][0] * 128 + Hsbit[i][1] * 64 + Hsbit[i][2] * 32 + Hsbit[i][3] * 16
            + Lsbit[i][0] * 8 + Lsbit[i][1] * 4 + Lsbit[i][2] * 2 + Lsbit[i][3] * 1;

        P[i] = (resimg + 256 - av[i]) % 256;
        /*
        if (i == 1 || i == 2 || i == 4 || i == 8 || i == 16 || i == 32 || i == 64 || i == 128 || i == 256 || i == 512 || i == 1024)
        {
            cout << i << " " << resimg << " " << P[i] << endl;
        }
        */
    }
    
}

//加密解密主函数
int main()
{
    /*
    Mat image = imread("./image/Lena.jpg");
    //Mat image = imread("./image/test1.jpeg");
    //Mat image = imread("./image/Peppers.jpg");
    n = image.rows;
    m = image.cols;
    int cns23 = 0, cns22 = 0, cns21 = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            mr[cns23++] = image.at<Vec3b>(i, j)[2];
            mg[cns22++] = image.at<Vec3b>(i, j)[1];
            mb[cns21++] = image.at<Vec3b>(i, j)[0];
        }
    }

    //设置密钥(初始值)
    SumR = 29676;
    SumG = 9202;
    SumB = 62299;

    init();

    encryption(mr, cr);
    encryption(mg, cg);
    encryption(mb, cb);

    cns23 = 0; cns22 = 0; cns21 = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            image.at<Vec3b>(i, j)[2] = cr[cns23++];
            image.at<Vec3b>(i, j)[1] = cg[cns22++];
            image.at<Vec3b>(i, j)[0] = cb[cns21++];
        }
    }
    waitKey(0);
    imwrite("./image/Tempchiper.bmp", image);
    //printf("%d %d %d %d %d\n",s21[10].id, s21[558].id, s21[640].id, s21[780].id, s21[110].id);
    */
    
    //设置密钥(初始值)
    SumR = 29676;
    SumG = 9202;
    SumB = 62299;

    Mat image = imread("./image/Tempchiper.bmp");

    n = image.rows;
    m = image.cols;
    int cns23 = 0, cns22 = 0, cns21 = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cr[cns23++] = image.at<Vec3b>(i, j)[2];
            cg[cns22++] = image.at<Vec3b>(i, j)[1];
            cb[cns21++] = image.at<Vec3b>(i, j)[0];
        }
    }

    init();
    
    CPAinit();
    atts2();
    attvl();
    atts1();
    attvh();
    atts4();
    atts3();
    attxor();

    CPA(cr, ar);
    CPA(cg, ag);
    CPA(cb, ab);

    cns23 = 0; cns22 = 0; cns21 = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            image.at<Vec3b>(i, j)[2] = ar[cns23++];
            image.at<Vec3b>(i, j)[1] = ag[cns22++];
            image.at<Vec3b>(i, j)[0] = ab[cns21++];
        }
    }
    imwrite("./image/Tempattack.jpg", image);
    waitKey(0);

}


//lena 29676 9202 62299
//peppers 29232 54749 57603