#include <bits/stdc++.h>
#ifdef _WIN32
#include <windows.h>
#endif
using namespace std;

#define RST "\033[0m"
#define BLD "\033[1m"
#define CYN "\033[36m"
#define GRN "\033[32m"
#define YLW "\033[33m"
#define RED "\033[31m"
#define GRY "\033[90m"

const double PI  = acos(-1.0);
const double EPS = 1e-9;

struct Pt { double x, y; };
struct Ln { double A, B, C; };
struct Cir { Pt O; double r; };

// ==================== 基础函数 ====================
Ln linePP(Pt p, Pt q) {
    return {p.y - q.y, q.x - p.x, (q.y - p.y) * p.x + (p.x - q.x) * p.y};
}
Pt interLL(Ln L, Ln M) {
    double D = L.A * M.B - L.B * M.A;
    if (fabs(D) < EPS) return {NAN, NAN};
    return {(L.B * M.C - L.C * M.B) / D, (L.C * M.A - L.A * M.C) / D};
}
double distPL(Pt p, Ln L) {
    return fabs(L.A * p.x + L.B * p.y + L.C) / sqrt(L.A * L.A + L.B * L.B);
}
Ln perpBisector(Pt p, Pt q) {
    Pt m{(p.x + q.x) / 2, (p.y + q.y) / 2};
    return {q.x - p.x, q.y - p.y, -((q.x - p.x) * m.x + (q.y - p.y) * m.y)};
}
double dist2(Pt a, Pt b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

// ==================== 圆相关 ====================
bool circumcircle(Pt A, Pt B, Pt C, Cir &res) {
    Pt O = interLL(perpBisector(A, B), perpBisector(A, C));
    if (isnan(O.x)) return false;
    res.O = O;
    res.r = sqrt(dist2(O, A));
    return true;
}
int interLC(Ln L, Cir C, Pt pts[2]) {
    double d = distPL(C.O, L);
    if (d > C.r + EPS) return 0;
    double den = L.A * L.A + L.B * L.B;
    double t = (L.A * C.O.x + L.B * C.O.y + L.C) / den;
    Pt P0 = {C.O.x - t * L.A, C.O.y - t * L.B};
    if (fabs(d - C.r) < EPS) { pts[0] = P0; return 1; }
    double h = sqrt(max(0.0, C.r * C.r - d * d));
    double len = sqrt(den);
    double ux = L.B / len, uy = -L.A / len;
    pts[0] = {P0.x + h * ux, P0.y + h * uy};
    pts[1] = {P0.x - h * ux, P0.y - h * uy};
    return 2;
}
string circleRelation(Cir C1, Cir C2) {
    double d = sqrt(dist2(C1.O, C2.O));
    if (d < EPS) return fabs(C1.r - C2.r) < EPS ? "重合" : "同心";
    double sum = C1.r + C2.r, diff = fabs(C1.r - C2.r);
    if (fabs(d - sum) < EPS) return "外切";
    if (fabs(d - diff) < EPS) return "内切";
    if (d > sum) return "外离";
    if (d < diff) return "内含";
    return "相交";
}
int interCC(Cir C1, Cir C2, Pt pts[2]) {
    if (dist2(C1.O, C2.O) < EPS) return 0;
    Ln rad;
    rad.A = 2 * (C2.O.x - C1.O.x);
    rad.B = 2 * (C2.O.y - C1.O.y);
    rad.C = (C1.O.x * C1.O.x + C1.O.y * C1.O.y - C1.r * C1.r)
          - (C2.O.x * C2.O.x + C2.O.y * C2.O.y - C2.r * C2.r);
    return interLC(rad, C1, pts);
}
int tangentsFromPoint(Pt P, Cir C, Ln tans[2]) {
    double d2 = dist2(P, C.O);
    if (d2 < C.r * C.r - EPS) return 0;
    if (fabs(d2 - C.r * C.r) < EPS) {
        double dx = P.x - C.O.x, dy = P.y - C.O.y;
        tans[0] = {dx, dy, -(dx * P.x + dy * P.y)};
        return 1;
    }
    Cir aux{P, sqrt(d2 - C.r * C.r)};
    Pt pts[2];
    int cnt = interCC(C, aux, pts);
    for (int i = 0; i < cnt; i++) tans[i] = linePP(P, pts[i]);
    return cnt;
}

// ==================== 其他功能 ====================
double distLL(Ln L, Ln M) {
    double len1 = sqrt(L.A * L.A + L.B * L.B);
    double len2 = sqrt(M.A * M.A + M.B * M.B);
    double s = (L.A * M.A + L.B * M.B) >= 0 ? 1 : -1;
    return fabs(L.C / len1 - s * M.C / len2);
}
Pt reflectPt(Pt P, Ln L) {
    double d = L.A * L.A + L.B * L.B;
    double t = (L.A * P.x + L.B * P.y + L.C) / d;
    return {P.x - 2 * t * L.A, P.y - 2 * t * L.B};
}
double angleLL(Ln L, Ln M) {
    double c = fabs(L.A * M.A + L.B * M.B)
             / (sqrt(L.A * L.A + L.B * L.B) * sqrt(M.A * M.A + M.B * M.B));
    if (c > 1) c = 1;
    return acos(c) * 180 / PI;
}
double chordLC(Ln L, Cir C) {
    double d = distPL(C.O, L);
    if (d > C.r + EPS) return -1;
    return 2 * sqrt(max(0.0, C.r * C.r - d * d));
}
double triArea(Pt A, Pt B, Pt C) {
    return fabs((B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y)) / 2;
}
Pt centroid(Pt A, Pt B, Pt C) {
    return {(A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3};
}
Pt orthocenter(Pt A, Pt B, Pt C) {
    Cir c;
    if (!circumcircle(A, B, C, c)) return {NAN, NAN};
    return {A.x + B.x + C.x - 2 * c.O.x, A.y + B.y + C.y - 2 * c.O.y};
}
Ln commonChord(Cir C1, Cir C2) {
    Ln L;
    L.A = 2 * (C2.O.x - C1.O.x);
    L.B = 2 * (C2.O.y - C1.O.y);
    L.C = (C1.O.x * C1.O.x + C1.O.y * C1.O.y - C1.r * C1.r)
        - (C2.O.x * C2.O.x + C2.O.y * C2.O.y - C2.r * C2.r);
    return L;
}
void distPR(Pt P, Cir C, double &nd, double &fd) {
    double d = sqrt(dist2(P, C.O));
    nd = fabs(d - C.r);
    fd = d + C.r;
}

// ==================== 界面（方法一：> 和 :）====================
void initConsole() {
#ifdef _WIN32
    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode;
    GetConsoleMode(h, &mode);
    SetConsoleMode(h, mode | 0x0004);
#endif
    cout << fixed << setprecision(4);
}
void cls() {
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}
void pauseScreen() {
    cout << GRY "\n  按回车返回菜单..." RST;
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cin.get();
}
void groupTitle(const string &s) {
    cout << "\n  " CYN "* " BLD << s << RST "\n";
}
void menuItem(int n, const string &s) {
    cout << "     " YLW "[" << setw(2) << n << "]" RST "  " << s << "\n";
}
void section(const string &s) {
    cout << "\n  " GRN "> " BLD << s << RST "\n";
    cout << GRY "  --------------------------------------------\n" RST "\n";
}

bool readVal(double &v) {
    while (true) {
        string s;
        if (!(cin >> s)) return false;
        if (s == "q" || s == "Q") {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            return false;
        }
        try { v = stod(s); return true; }
        catch (...) { cout << "  " RED "无效，请重输或 q: " RST; }
    }
}
bool ip(Pt &p, const string &n) {
    cout << "  " GRY << n << RST " x y: " GRY "(q 返回)" RST "  ";
    return readVal(p.x) && readVal(p.y);
}
bool ic(Cir &c, const string &n) {
    cout << "  " GRY << n << RST " 圆心x 圆心y r: " GRY "(q 返回)" RST "  ";
    return readVal(c.O.x) && readVal(c.O.y) && readVal(c.r);
}
bool il(Ln &L, const string &n) {
    cout << "  " GRY << n << RST " A B C: " GRY "(q 返回)" RST "  ";
    return readVal(L.A) && readVal(L.B) && readVal(L.C);
}

void showHelp() {
    cls();
    cout << "\n";
    cout << CYN BLD "  +--------------------------------------------+\n" RST;
    cout << CYN BLD "  |            使 用 说 明                    |\n" RST;
    cout << CYN BLD "  +--------------------------------------------+\n" RST;
    cout << "\n";
    cout << "  " GRN "> 基本操作\n" RST;
    cout << "      - 菜单输入数字选择功能，输 0 退出\n";
    cout << "      - 输入数据时用空格分隔，例如：1 2 表示 (1, 2)\n";
    cout << "      - 输入过程中随时输入 " YLW "q" RST " 可返回菜单\n";
    cout << "      - 每个功能算完按回车返回菜单\n\n";
    cout << "  " GRN "> 输入格式\n" RST;
    cout << "      - 点       " GRY "x y" RST "          例:  1 2\n";
    cout << "      - 直线     " GRY "A B C" RST "        例:  1 -1 0   即 x - y = 0\n";
    cout << "      - 圆       " GRY "圆心x 圆心y r" RST "  例:  0 0 5   即 x2+y2=25\n\n";
    cout << "  " GRN "> 功能速查\n" RST;
    cout << "      - 直线      求方程 交点 距离 中垂线 对称点 夹角\n";
    cout << "      - 圆        外接圆 直线与圆 两圆关系 切线 弦长 公共弦\n";
    cout << "      - 三角形    面积 重心 垂心\n\n";
    cout << "  " GRY "按回车开始使用..." RST;
    cin.get();
    cls();
}

void menu() {
    cout << "\n";
    cout << CYN BLD "  ========== 解析几何求解器 v3.0 ==========\n" RST;

    groupTitle("直线");
    menuItem(1,  "两点求直线");
    menuItem(2,  "两直线交点");
    menuItem(3,  "点到直线距离");
    menuItem(4,  "两点中垂线");
    menuItem(9,  "两平行线距离");
    menuItem(10, "点关于直线对称");
    menuItem(11, "两直线夹角");

    groupTitle("圆");
    menuItem(5,  "三点求外接圆");
    menuItem(6,  "直线与圆");
    menuItem(7,  "两圆关系");
    menuItem(8,  "点到圆的切线");
    menuItem(12, "直线截圆弦长");
    menuItem(16, "两圆公共弦");
    menuItem(17, "点到圆最远/最近距离");

    groupTitle("三角形");
    menuItem(13, "三角形面积");
    menuItem(14, "三角形重心");
    menuItem(15, "三角形垂心");

    cout << "\n";
    menuItem(99, "使用说明");
    menuItem(0,  "退出");
}

// ==================== 主程序 ====================
int main() {
    initConsole();
    showHelp();

    const char *FN[18] = {
        "",
        "两点求直线", "两直线交点", "点到直线距离", "两点中垂线",
        "三点求外接圆", "直线与圆", "两圆关系", "点到圆的切线",
        "两平行线距离", "点关于直线对称", "两直线夹角", "直线截圆弦长",
        "三角形面积", "三角形重心", "三角形垂心", "两圆公共弦",
        "点到圆最远/最近距离"
    };

    while (true) {
        cls();
        menu();
        cout << "\n  " GRN "请选择: " RST;

        int c;
        if (!(cin >> c)) {
            cin.clear();
            cin.ignore(1000, '\n');
            continue;
        }
        cin.ignore(1000, '\n');

        if (c == 0) break;
        if (c == 99) { showHelp(); continue; }
        if (c < 1 || c > 17) {
            cout << "\n  " RED "无效选项" RST "\n";
            pauseScreen();
            continue;
        }

        cls();
        section(FN[c]);
        cout << "  " GRY "输入 " YLW "q" GRY " 可随时返回菜单" RST "\n\n";

        switch (c) {
        case 1: {
            Pt p, q;
            if (!ip(p, "P1")) continue;
            if (!ip(q, "P2")) continue;
            Ln L = linePP(p, q);
            cout << "\n  " GRN "> 直线: " RST
                 << L.A << "x + " << L.B << "y + " << L.C << " = 0\n";
            break;
        }
        case 2: {
            Ln L, M;
            if (!il(L, "L1")) continue;
            if (!il(M, "L2")) continue;
            Pt p = interLL(L, M);
            if (isnan(p.x)) cout << "\n  " GRN "> 平行" RST "\n";
            else cout << "\n  " GRN "> 交点: (" RST
                      << p.x << ", " << p.y << ")\n";
            break;
        }
        case 3: {
            Pt p; Ln L;
            if (!ip(p, "P")) continue;
            if (!il(L, "直线")) continue;
            cout << "\n  " GRN "> 距离: " RST << distPL(p, L) << "\n";
            break;
        }
        case 4: {
            Pt p, q;
            if (!ip(p, "P1")) continue;
            if (!ip(q, "P2")) continue;
            Ln L = perpBisector(p, q);
            cout << "\n  " GRN "> 中垂线: " RST
                 << L.A << "x + " << L.B << "y + " << L.C << " = 0\n";
            break;
        }
        case 5: {
            Pt A, B, C;
            if (!ip(A, "A")) continue;
            if (!ip(B, "B")) continue;
            if (!ip(C, "C")) continue;
            Cir r;
            if (!circumcircle(A, B, C, r))
                cout << "\n  " RED "> 三点共线" RST "\n";
            else
                cout << "\n  " GRN "> 圆心: (" RST
                     << r.O.x << ", " << r.O.y << ")\n"
                     << GRN "    半径: " RST << r.r << "\n";
            break;
        }
        case 6: {
            Pt p, q; Cir C;
            if (!ip(p, "直线上点P1")) continue;
            if (!ip(q, "直线上点P2")) continue;
            Ln L = linePP(p, q);
            if (!ic(C, "圆")) continue;
            double d = distPL(C.O, L);
            cout << "\n  " GRY "d = " RST << d << GRY ",  r = " RST << C.r << "\n";
            if (d > C.r + EPS) cout << "  " GRN "> 相离" RST "\n";
            else if (fabs(d - C.r) < EPS) {
                Pt pts[2]; interLC(L, C, pts);
                cout << "  " GRN "> 相切，切点: " RST
                     << "(" << pts[0].x << ", " << pts[0].y << ")\n";
            } else {
                Pt pts[2]; interLC(L, C, pts);
                cout << "  " GRN "> 相交" RST "\n";
                cout << "    交点1: (" << pts[0].x << ", " << pts[0].y << ")\n";
                cout << "    交点2: (" << pts[1].x << ", " << pts[1].y << ")\n";
            }
            break;
        }
        case 7: {
            Cir C1, C2;
            if (!ic(C1, "圆1")) continue;
            if (!ic(C2, "圆2")) continue;
            cout << "\n  " GRN "> 位置关系: " RST << circleRelation(C1, C2) << "\n";
            Pt pts[2];
            int cnt = interCC(C1, C2, pts);
            for (int i = 0; i < cnt; i++)
                cout << "    交点" << i + 1 << ": ("
                     << pts[i].x << ", " << pts[i].y << ")\n";
            break;
        }
        case 8: {
            Pt P; Cir C;
            if (!ip(P, "P")) continue;
            if (!ic(C, "圆")) continue;
            Ln t[2];
            int cnt = tangentsFromPoint(P, C, t);
            cout << "\n";
            if (cnt == 0) cout << "  " GRN "> 点在圆内，无切线" RST "\n";
            else for (int i = 0; i < cnt; i++)
                cout << "  " GRN "> 切线" << (cnt > 1 ? to_string(i + 1) : "")
                     << ": " RST << t[i].A << "x + " << t[i].B
                     << "y + " << t[i].C << " = 0\n";
            break;
        }
        case 9: {
            Ln L, M;
            if (!il(L, "L1")) continue;
            if (!il(M, "L2")) continue;
            if (fabs(L.A * M.B - L.B * M.A) > EPS)
                cout << "\n  " RED "> 两直线不平行" RST "\n";
            else
                cout << "\n  " GRN "> 距离: " RST << distLL(L, M) << "\n";
            break;
        }
        case 10: {
            Pt p; Ln L;
            if (!ip(p, "P")) continue;
            if (!il(L, "直线")) continue;
            Pt Q = reflectPt(p, L);
            cout << "\n  " GRN "> 对称点: (" RST << Q.x << ", " << Q.y << ")\n";
            break;
        }
        case 11: {
            Ln L, M;
            if (!il(L, "L1")) continue;
            if (!il(M, "L2")) continue;
            cout << "\n  " GRN "> 夹角: " RST << angleLL(L, M) << "°\n";
            break;
        }
        case 12: {
            Pt p, q; Cir C;
            if (!ip(p, "直线上点P1")) continue;
            if (!ip(q, "直线上点P2")) continue;
            Ln L = linePP(p, q);
            if (!ic(C, "圆")) continue;
            double ch = chordLC(L, C);
            if (ch < 0) cout << "\n  " GRN "> 相离，无弦" RST "\n";
            else cout << "\n  " GRN "> 弦长: " RST << ch << "\n";
            break;
        }
        case 13: {
            Pt A, B, C;
            if (!ip(A, "A")) continue;
            if (!ip(B, "B")) continue;
            if (!ip(C, "C")) continue;
            cout << "\n  " GRN "> 面积: " RST << triArea(A, B, C) << "\n";
            break;
        }
        case 14: {
            Pt A, B, C;
            if (!ip(A, "A")) continue;
            if (!ip(B, "B")) continue;
            if (!ip(C, "C")) continue;
            Pt G = centroid(A, B, C);
            cout << "\n  " GRN "> 重心: (" RST << G.x << ", " << G.y << ")\n";
            break;
        }
        case 15: {
            Pt A, B, C;
            if (!ip(A, "A")) continue;
            if (!ip(B, "B")) continue;
            if (!ip(C, "C")) continue;
            Pt H = orthocenter(A, B, C);
            if (isnan(H.x)) cout << "\n  " RED "> 三点共线" RST "\n";
            else cout << "\n  " GRN "> 垂心: (" RST << H.x << ", " << H.y << ")\n";
            break;
        }
        case 16: {
            Cir C1, C2;
            if (!ic(C1, "圆1")) continue;
            if (!ic(C2, "圆2")) continue;
            string rel = circleRelation(C1, C2);
            cout << "\n  " GRN "> 位置关系: " RST << rel << "\n";
            if (rel == "相交" || rel == "外切" || rel == "内切") {
                Ln L = commonChord(C1, C2);
                cout << "    公共弦: " << L.A << "x + " << L.B
                     << "y + " << L.C << " = 0\n";
            } else {
                cout << "    无公共弦\n";
            }
            break;
        }
        case 17: {
            Pt P; Cir C;
            if (!ip(P, "P")) continue;
            if (!ic(C, "圆")) continue;
            double nd, fd;
            distPR(P, C, nd, fd);
            cout << "\n  " GRN "> 最近距离: " RST << nd << "\n";
            cout << "  " GRN "> 最远距离: " RST << fd << "\n";
            break;
        }
        }
        pauseScreen();
    }

    cls();
    cout << "\n\n  " CYN BLD "  感谢使用，再见！" RST "\n\n";
    return 0;
}
