// ConsoleApplication1.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include <stdio.h>
#include <math.h>
#include <iostream> //入出力ライブラリ
#include <fstream> //iostreamのファイル入出力をサポート

using namespace std;

//流体の条件
double  re = 70.0;
double  cfl = 0.2;
double  omegap = 1.00;
int     maxitp = 100;
double  errorp = 1.0e-4;
int     nlast = 5000;
int     nlp = 10;

//計算格子の設定
#define mx  401
int     i_1 = 96;
int     i_2 = 106;

#define my  201
int     j_1 = 96;
int     j_2 = 106;

//関数の宣言
void intcnd(double u_dummy[][my + 2], double v_dummy[][my + 2], double p_dummy[][my + 2]);
void bcforp(double p_dummy[][my + 2]);
void bcforv(double u_dummy[][my + 2], double v_dummy[][my + 2]);
int poiseq(double u_dummy[][my + 2], double v_dummy[][my + 2], double p_dummy[][my + 2], double dx_dummy, double dy_dummy, double dt_dummy);
void veloeq(double u_dummy[][my + 2], double v_dummy[][my + 2], double p_dummy[][my + 2], double dx_dummy, double dy_dummy, double dt_dummy);

//template <typename T>
//void p(T s) {
//    std::cout << s << std::endl;
//}

int main()
{
    //計算格子の設定（物体）
    double  dx = 1.0 / (double)(i_2 - i_1);
    double  dy = 1.0 / (double)(j_2 - j_1);
    int     icent = (i_1 + i_2) / 2;
    int     jcent = (j_1 + j_2) / 2;

    double  x[mx + 1][my + 1];
    double  y[mx + 1][my + 1];

    {
        for (int i = 1; i < mx + 1; i++) {
            for (int j = 1; j < my + 1; j++) {
                x[i][j] = dx*(double)(i - icent);
                y[i][j] = dy*(double)(j - jcent);
            }
        }
    }

    //{
    //    for (int i; i < mx; i++) {
    //        x[i] = dx*(double)(i - icent);
    //    }
    //    for (int j; j < my; j++) {
    //        y[j] = dy*(double)(j - jcent);
    //    }
    //}
    //計算格子の設定終了

    //クーラン条件
    double dt = cfl*fmin(dx, dy);

    ////////////////////////////////////////////////////////////////////////////
    //流れを解く
    ////////////////////////////////////////////////////////////////////////////
    double  u[mx + 2][my + 2];
    double  v[mx + 2][my + 2];
    double  p[mx + 2][my + 2];

    //初期条件の設定
    //int     n = 0;
    //double  time = 0.0;
    intcnd(u, v, p);
    bcforp(p);
    bcforv(u, v);

    //初期化成功か見る
    /*for (int i = 0; i < mx; i++) {
    for (int j = 0; j < my; j++) {
    printf("%f", u[i][j]);
    }
    putchar('\n');
    }*/

    int itr_num;

    //サイクル計算
    for (int n = 1; n < nlast + 1; n++) {

        itr_num = poiseq(u, v, p, dx, dy, dt);
        bcforp(p);

        veloeq(u, v, p, dx, dy, dt);
        bcforv(u, v);

        if (n % 10 == 0) {
            printf("%d,%d\n", n, itr_num);
        }
    }

    //圧力を正しい値に変換する
    for (int i = 1; i < mx + 1; i++) {
        for (int j = 1; j < my + 1; j++) {
            p[i][j] = 2.0*p[i][j];
        }
    }

    //.dataに書き出し
    {
        ofstream ofs("Result.data");
        for (int i = 1; i < mx + 1; i++) {
            for (int j = 1; j < my + 1; j++) {
                ofs << x[i][j] << " " << y[i][j] << " " << 2.0*p[i][j] << endl;
            }
            ofs << endl;
        }
        ofs << endl;
    }
    //csvに書き出し
    {
        ofstream ofs("Result_x.csv");
        for (int i = 1; i < mx + 1; i++) {
            ofs << x[i][1] << ",";
        }
        ofs << endl;
    }
    {
        ofstream ofs("Result_y.csv");
        for (int j = 1; j < my + 1; j++) {
            ofs << y[1][j] << ",";
        }
        ofs << endl;
    }
    {
        ofstream ofs("Result_p.csv");
        for (int i = 1; i < mx + 1; i++) {
            for (int j = 1; j < my + 1; j++) {
                ofs << 2.0*p[i][j] << ",";
            }
            ofs << endl;
        }
        ofs << endl;
    }

    double cd = 0.0;
    {
        for (int j = j_1; j < j_2; j++) {
            static double cpfore = (2.0*p[i_1][j] + 2.0*p[i_1][j + 1]) / 2.0;
            static double cpback = (2.0*p[i_2][j] + 2.0*p[i_2][j + 1]) / 2.0;
            cd = cd + (cpfore - cpback)*dy;
        }
    }
    printf("%f", cd);

    double cl = 0.0;
    {
        for (int i = i_1; i < i_2; i++) {
            static double cpbtm = (2.0*p[i][j_1] + 2.0*p[i + 1][j_1]) / 2.0;
            static double cptop = (2.0*p[i][j_2] + 2.0*p[i + 1][j_2]) / 2.0;
            cl = cl + (cpbtm - cptop)*dx;
        }
    }
    printf("%f", cl);

    printf("成功\n");
    return 0;

}

/*初期条件の設定*/
void intcnd(double u_dummy[][my + 2], double v_dummy[][my + 2], double p_dummy[][my + 2])
{
    for (int i = 0; i < mx + 2; i++) {
        for (int j = 0; j < my + 2; j++) {
            u_dummy[i][j] = 1.0;
            v_dummy[i][j] = 0.0;
            p_dummy[i][j] = 0.0;
        }
    }
}

/*圧力の境界条件の設定*/
void bcforp(double p_dummy[][my + 2])
{
    {
        int i = 1;
        for (int j = 1; j < my + 1; j++) {
            p_dummy[i][j] = 0.0;
        }
    }
    {
        int i = mx;
        for (int j = 1; j < my + 1; j++) {
            p_dummy[i][j] = 0.0;
        }
    }
    {
        int j = 1;
        for (int i = 1; i < mx + 1; i++) {
            p_dummy[i][j] = 0.0;
        }
    }
    {
        int j = my;
        for (int i = 1; i < mx + 1; i++) {
            p_dummy[i][j] = 0.0;
        }
    }

    //壁周りの境界条件
    p_dummy[i_1][j_1] = p_dummy[i_1 - 1][j_1 - 1];
    p_dummy[i_1][j_2] = p_dummy[i_1 - 1][j_2 + 1];
    p_dummy[i_2][j_1] = p_dummy[i_2 + 1][j_1 - 1];
    p_dummy[i_2][j_2] = p_dummy[i_2 + 1][j_2 + 1];
    {
        int i = i_1;
        for (int j = j_1 + 1; j < j_2; j++) {
            p_dummy[i][j] = p_dummy[i - 1][j];
        }
    }
    {
        int i = i_2;
        for (int j = j_1 + 1; j < j_2; j++) {
            p_dummy[i][j] = p_dummy[i + 1][j];
        }
    }
    {
        int j = j_1;
        for (int i = i_1 + 1; i < i_2; i++) {
            p_dummy[i][j] = p_dummy[i][j - 1];
        }
    }
    {
        int j = j_2;
        for (int i = i_1 + 1; i < i_2; i++) {
            p_dummy[i][j] = p_dummy[i][j + 1];
        }
    }
}

/*速度の境界条件の設定*/
void bcforv(double u_dummy[][my + 2], double v_dummy[][my + 2])
{
    {
        int i = 1;
        for (int j = 1; j < my + 1; j++) {
            u_dummy[i][j] = 1.0;
            v_dummy[i][j] = 0.0;
            u_dummy[i - 1][j] = 1.0;
            v_dummy[i - 1][j] = 0.0;
        }
    }

    {
        int i = mx;
        for (int j = 1; j < my + 1; j++) {
            u_dummy[i][j] = 2.0*u_dummy[i - 1][j] - u_dummy[i - 2][j];
            v_dummy[i][j] = 2.0*v_dummy[i - 1][j] - v_dummy[i - 2][j];
            u_dummy[i + 1][j] = 2.0 * u_dummy[i][j] - u_dummy[i - 1][j];
            v_dummy[i + 1][j] = 2.0 * v_dummy[i][j] - v_dummy[i - 1][j];
        }
    }

    {
        int j = 1;
        for (int i = 1; i < mx + 1; i++) {
            u_dummy[i][j] = 2.0*u_dummy[i][j + 1] - u_dummy[i][j + 2];
            v_dummy[i][j] = 2.0*v_dummy[i][j + 1] - v_dummy[i][j + 2];
            u_dummy[i][j - 1] = 2.0*u_dummy[i][j] - u_dummy[i][j + 1];
            v_dummy[i][j - 1] = 2.0*v_dummy[i][j] - v_dummy[i][j + 1];
        }
    }
    {
        int j = my;
        for (int i = 1; i < mx + 1; i++) {
            u_dummy[i][j] = 2.0*u_dummy[i][j - 1] - u_dummy[i][j - 2];
            v_dummy[i][j] = 2.0*v_dummy[i][j - 1] - v_dummy[i][j - 2];
            u_dummy[i][j + 1] = 2.0*u_dummy[i][j] - u_dummy[i][j - 1];
            v_dummy[i][j + 1] = 2.0*v_dummy[i][j] - v_dummy[i][j - 1];
        }
    }
    for (int i = i_1; i < i_2 + 1; i++) {
        for (int j = j_1; j < j_2 + 1; j++) {
            u_dummy[i][j] = 0.0;
            v_dummy[i][j] = 0.0;
        }
    }

}

/*圧力場を解く*/
int poiseq(double u_dummy[][my + 2], double v_dummy[][my + 2], double p_dummy[][my + 2], double dx_dummy, double dy_dummy, double dt_dummy)
{
    static double rhs[mx + 2][my + 2];

    //RHSを計算する
    {
        for (int i = 2; i < mx; i++) {
            for (int j = 2; j < my; j++) {
                if ((i_1 <= i && i <= i_2) && (j_1 <= j && j <= j_2)) {
                    continue;
                }
                else {
                    double ux = (u_dummy[i + 1][j] - u_dummy[i - 1][j]) / (2.0*dx_dummy);
                    double uy = (u_dummy[i][j + 1] - u_dummy[i][j - 1]) / (2.0*dy_dummy);
                    double vx = (v_dummy[i + 1][j] - v_dummy[i - 1][j]) / (2.0*dx_dummy);
                    double vy = (v_dummy[i][j + 1] - v_dummy[i][j - 1]) / (2.0*dy_dummy);

                    rhs[i][j] = (ux + vy) / dt_dummy - (pow(ux, 2.0) + 2.0*uy*vx + pow(vy, 2.0));
                }
            }
        }
    }
    //iteraionを回す
    double res = 0.0;
    int itr_final = 0;
    {

        for (int itr = 1; itr < maxitp + 1; itr++) {
            res = 0.0;

            for (int i = 2; i < mx; i++) {
                for (int j = 2; j < my; j++) {
                    //物体内部の時は計算しない
                    if ((i_1 <= i && i <= i_2) && (j_1 <= j && j <= j_2)) {
                        //何もしない
                    }
                        //それ以外の時
                    else {
                        double dp;
                        dp = (p_dummy[i + 1][j] + p_dummy[i - 1][j]) / pow(dx_dummy, 2.0)
                             + (p_dummy[i][j + 1] + p_dummy[i][j - 1]) / pow(dy_dummy, 2.0)
                             - rhs[i][j];
                        dp = dp / (2.0 / pow(dx_dummy, 2.0) + 2.0 / pow(dy_dummy, 2.0)) - p_dummy[i][j];
                        res = res + pow(dp, 2.0);
                        p_dummy[i][j] = p_dummy[i][j] + omegap*dp;
                    }

                }

            }
            //境界条件の設定
            bcforp(p_dummy);

            res = sqrt(res / ((double)(mx*my)));
            itr_final = itr;
            //printf("%f", res);

            //残差が十分小さいときループから脱出
            if (res < errorp) {
                //printf("脱出します");
                break;
            }
            itr_final = itr;
        }

    }
    return itr_final;
}

void veloeq(double u_dummy[][my + 2], double v_dummy[][my + 2], double p_dummy[][my + 2], double dx_dummy, double dy_dummy, double dt_dummy)
{
    static double urhs[mx + 2][my + 2];
    static double vrhs[mx + 2][my + 2];

    //圧力勾配
    {
        for (int i = 2; i < mx; i++) {
            for (int j = 2; j < my; j++) {
                if ((i_1 <= i && i <= i_2) && (j_1 <= j && j <= j_2)) {
                    continue;
                }
                else {
                    urhs[i][j] = -(p_dummy[i + 1][j] - p_dummy[i - 1][j]) / (2.0*dx_dummy);
                    vrhs[i][j] = -(p_dummy[i][j + 1] - p_dummy[i][j - 1]) / (2.0*dy_dummy);
                }
            }
        }
    }

    //粘性項
    {
        for (int i = 2; i < mx; i++) {
            for (int j = 2; j < my; j++) {
                if ((i_1 <= i && i <= i_2) && (j_1 <= j && j <= j_2)) {
                    continue;
                }
                else {
                    urhs[i][j] = urhs[i][j]
                                 + (u_dummy[i + 1][j] - 2.0*u_dummy[i][j] + u_dummy[i - 1][j]) / (re*pow(dx_dummy, 2.0))
                                 + (u_dummy[i][j + 1] - 2.0*u_dummy[i][j] + u_dummy[i][j - 1]) / (re*pow(dy_dummy, 2.0));
                    vrhs[i][j] = vrhs[i][j]
                                 + (v_dummy[i + 1][j] - 2.0*v_dummy[i][j] + v_dummy[i - 1][j]) / (re*pow(dx_dummy, 2.0))
                                 + (v_dummy[i][j + 1] - 2.0*v_dummy[i][j] + v_dummy[i][j - 1]) / (re*pow(dy_dummy, 2.0));
                }
            }
        }
    }

    //移流項（x軸について)
    {
        for (int j = j_1 + 1; j < j_2; j++) {
            u_dummy[i_1 + 1][j] = 2.0*u_dummy[i_1][j] - u_dummy[i_1 - 1][j];
            u_dummy[i_2 - 1][j] = 2.0*u_dummy[i_2][j] - u_dummy[i_2 + 1][j];
            v_dummy[i_1 + 1][j] = 2.0*v_dummy[i_1][j] - v_dummy[i_1 - 1][j];
            v_dummy[i_2 - 1][j] = 2.0*v_dummy[i_2][j] - v_dummy[i_2 + 1][j];
        }

    }
    {
        for (int i = 2; i < mx; i++) {
            for (int j = 2; j < my; j++) {
                if ((i_1 <= i && i <= i_2) && (j_1 <= j && j <= j_2)) {
                    continue;
                }
                else {
                    urhs[i][j] = urhs[i][j]
                                 - u_dummy[i][j]
                                   * (-u_dummy[i + 2][j] + 8.0*(u_dummy[i + 1][j] - u_dummy[i - 1][j]) + u_dummy[i - 2][j])
                                   / (12.0*dx_dummy)
                                 - fabs(u_dummy[i][j])
                                   *(u_dummy[i + 2][j] - 4.0*u_dummy[i + 1][j] + 6.0*u_dummy[i][j] - 4.0*u_dummy[i - 1][j] + u_dummy[i - 2][j])
                                   / (4.0*dx_dummy);

                    vrhs[i][j] = vrhs[i][j]
                                 - u_dummy[i][j]
                                   * (-v_dummy[i + 2][j] + 8.0*(v_dummy[i + 1][j] - v_dummy[i - 1][j]) + v_dummy[i - 2][j])
                                   / (12.0*dx_dummy)
                                 - fabs(u_dummy[i][j])
                                   *(v_dummy[i + 2][j] - 4.0*v_dummy[i + 1][j] + 6.0*v_dummy[i][j] - 4.0*v_dummy[i - 1][j] + v_dummy[i - 2][j])
                                   / (4.0*dx_dummy);
                }
            }
        }
    }

    //移流項（y軸）
    {
        for (int i = i_1 + 1; i < i_2; i++) {
            u_dummy[i][j_1 + 1] = 2.0*u_dummy[i][j_1] - u_dummy[i][j_1 - 1];
            u_dummy[i][j_2 - 1] = 2.0*u_dummy[i][j_2] - u_dummy[i][j_2 + 1];
            v_dummy[i][j_1 + 1] = 2.0*v_dummy[i][j_1] - v_dummy[i][j_1 - 1];
            v_dummy[i][j_2 - 1] = 2.0*v_dummy[i][j_2] - v_dummy[i][j_2 + 1];
        }

    }
    {
        for (int i = 2; i < mx; i++) {
            for (int j = 2; j < my; j++) {
                if ((i_1 <= i && i <= i_2) && (j_1 <= j && j <= j_2)) {
                    continue;
                }
                else {
                    urhs[i][j] = urhs[i][j]
                                 - v_dummy[i][j]
                                   * (-u_dummy[i][j + 2] + 8.0*(u_dummy[i][j + 1] - u_dummy[i][j - 1]) + u_dummy[i][j - 2])
                                   / (12.0*dy_dummy)
                                 - fabs(v_dummy[i][j])
                                   *(u_dummy[i][j + 2] - 4.0*u_dummy[i][j + 1] + 6.0*u_dummy[i][j] - 4.0*u_dummy[i][j - 1] + u_dummy[i][j - 2])
                                   / (4.0*dy_dummy);

                    vrhs[i][j] = vrhs[i][j]
                                 - v_dummy[i][j]
                                   * (-v_dummy[i][j + 2] + 8.0*(v_dummy[i][j + 1] - v_dummy[i][j - 1]) + v_dummy[i][j - 2])
                                   / (12.0*dy_dummy)
                                 - fabs(v_dummy[i][j])
                                   *(v_dummy[i][j + 2] - 4.0*v_dummy[i][j + 1] + 6.0*v_dummy[i][j] - 4.0*v_dummy[i][j - 1] + v_dummy[i][j - 2])
                                   / (4.0*dy_dummy);
                }
            }
        }
    }

    //update
    {
        for (int i = 2; i < mx; i++) {
            for (int j = 2; j < my; j++) {
                if ((i_1 <= i && i <= i_2) && (j_1 <= j && j <= j_2)) {
                    continue;
                }
                else {
                    u_dummy[i][j] = u_dummy[i][j] + dt_dummy*urhs[i][j];
                    v_dummy[i][j] = v_dummy[i][j] + dt_dummy*vrhs[i][j];
                }
            }
        }

    }

}