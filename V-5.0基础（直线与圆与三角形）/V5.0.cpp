#include <bits/stdc++.h>
using namespace std;

const double PI = acos(-1.0);

struct Pt { double x, y; };
struct Ln { double A, B, C; };   // Ax + By + C = 0
struct Cir { Pt O; double r; };

// ==================== 基础函数 ====================
Ln linePP(Pt p, Pt q) {
    return {p.y - q.y, q.x - p.x, (q.y - p.y) * p.x + (p.x - q.x) * p.y};
}

Pt interLL(Ln L, Ln M) {
    double D = L.A * M.B - L.B * M.A;
    if (fabs(D) < 1e-9) return {NAN, NAN};
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
    Ln L1 = perpBisector(A, B);
    Ln L2 = perpBisector(A, C);
    Pt O = interLL(L1, L2);
    if (isnan(O.x)) return false;
    res.O = O;
    res.r = sqrt(dist2(O, A));
    return true;
}

int interLC(Ln L, Cir C, Pt pts[2]) {
    double d = distPL(C.O, L);
    if (d > C.r + 1e-9) return 0;
    double den = L.A * L.A + L.B * L.B;
    double t = (L.A * C.O.x + L.B * C.O.y + L.C) / den;
    Pt P0 = {C.O.x - t * L.A, C.O.y - t * L.B};
    if (fabs(d - C.r) < 1e-9) {
        pts[0] = P0;
        return 1;
    }
    double h = sqrt(max(0.0, C.r * C.r - d * d));
    double len = sqrt(den);
    double ux = L.B / len, uy = -L.A / len;
    pts[0] = {P0.x + h * ux, P0.y + h * uy};
    pts[1] = {P0.x - h * ux, P0.y - h * uy};
    return 2;
}

string circleRelation(Cir C1, Cir C2) {
    double d = sqrt(dist2(C1.O, C2.O));
    if (d < 1e-9)
        return fabs(C1.r - C2.r) < 1e-9 ? "重合" : "同心";
    double sum = C1.r + C2.r, diff = fabs(C1.r - C2.r);
    if (fabs(d - sum) < 1e-9) return "外切";
    if (fabs(d - diff) < 1e-9) return "内切";
    if (d > sum) return "外离";
    if (d < diff) return "内含";
    return "相交";
}

int interCC(Cir C1, Cir C2, Pt pts[2]) {
    if (dist2(C1.O, C2.O) < 1e-9) return 0;
    double A = 2 * (C2.O.x - C1.O.x);
    double B = 2 * (C2.O.y - C1.O.y);
    double Cc = (C1.O.x * C1.O.x + C1.O.y * C1.O.y - C1.r * C1.r)
              - (C2.O.x * C2.O.x + C2.O.y * C2.O.y - C2.r * C2.r);
    Ln rad = {A, B, Cc};
    return interLC(rad, C1, pts);
}

int tangentsFromPoint(Pt P, Cir C, Ln tans[2]) {
    double d2 = dist2(P, C.O);
    if (d2 < C.r * C.r - 1e-9) return 0;
    if (fabs(d2 - C.r * C.r) < 1e-9) {
        double dx = P.x - C.O.x, dy = P.y - C.O.y;
        tans[0] = {dx, dy, -(dx * P.x + dy * P.y)};
        return 1;
    }
    Cir aux;
    aux.O = P;
    aux.r = sqrt(d2 - C.r * C.r);
    Pt pts[2];
    int cnt = interCC(C, aux, pts);
    for (int i = 0; i < cnt; i++)
        tans[i] = linePP(P, pts[i]);
    return cnt;
}

// ==================== 新增：其他常用功能 ====================

// 9. 两平行线距离（只取法向量同向的偏移）
double distLL(Ln L, Ln M) {
    double len1 = sqrt(L.A * L.A + L.B * L.B);
    double len2 = sqrt(M.A * M.A + M.B * M.B);
    double s = (L.A * M.A + L.B * M.B) >= 0 ? 1 : -1;
    double c1 = L.C / len1;
    double c2 = s * M.C / len2;
    return fabs(c1 - c2);
}

// 10. 点关于直线的对称点
Pt reflectPt(Pt P, Ln L) {
    double d = L.A * L.A + L.B * L.B;
    double t = (L.A * P.x + L.B * P.y + L.C) / d;
    return {P.x - 2 * t * L.A, P.y - 2 * t * L.B};
}

// 11. 两直线夹角（0~90 度）
double angleLL(Ln L, Ln M) {
    double dot = L.A * M.A + L.B * M.B;
    double len1 = sqrt(L.A * L.A + L.B * L.B);
    double len2 = sqrt(M.A * M.A + M.B * M.B);
    double c = fabs(dot) / (len1 * len2);
    if (c > 1) c = 1;
    return acos(c) * 180 / PI;
}

// 12. 直线截圆弦长（相离返回 -1）
double chordLC(Ln L, Cir C) {
    double d = distPL(C.O, L);
    if (d > C.r + 1e-9) return -1;
    return 2 * sqrt(max(0.0, C.r * C.r - d * d));
}

// 13. 三角形面积（向量叉积的一半）
double triArea(Pt A, Pt B, Pt C) {
    return fabs((B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y)) / 2;
}

// 14. 三角形重心
Pt centroid(Pt A, Pt B, Pt C) {
    return {(A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3};
}

// 15. 三角形垂心（用欧拉线：H = A+B+C - 2O）
Pt orthocenter(Pt A, Pt B, Pt C) {
    Cir cir;
    if (!circumcircle(A, B, C, cir)) return {NAN, NAN};
    return {A.x + B.x + C.x - 2 * cir.O.x,
            A.y + B.y + C.y - 2 * cir.O.y};
}

// 16. 两圆公共弦（根轴直线），注意只在相交/相切时有意义
Ln commonChord(Cir C1, Cir C2) {
    double A = 2 * (C2.O.x - C1.O.x);
    double B = 2 * (C2.O.y - C1.O.y);
    double Cc = (C1.O.x * C1.O.x + C1.O.y * C1.O.y - C1.r * C1.r)
              - (C2.O.x * C2.O.x + C2.O.y * C2.O.y - C2.r * C2.r);
    return {A, B, Cc};
}

// 17. 点到圆的最远、最近距离
void distPR(Pt P, Cir C, double &nearD, double &farD) {
    double d = sqrt(dist2(P, C.O));
    nearD = fabs(d - C.r);
    farD = d + C.r;
}

// ==================== 输入辅助 ====================
void ip(Pt &p, string n) {
    cout << n << " x y: ";
    cin >> p.x >> p.y;
}
void ic(Cir &c, string n) {
    cout << n << " 圆心 x y r: ";
    cin >> c.O.x >> c.O.y >> c.r;
}
void il(Ln &L, string n) {
    cout << n << " A B C: ";
    cin >> L.A >> L.B >> L.C;
}

// ==================== 主程序 ====================
int main() {
    cout << "===== 解析几何求解器 =====\n";
    for (;;) {
        cout << "\n============ 菜单 ============\n";
        cout << "  --- 直线 ---\n";
        cout << "  1. 两点求直线\n";
        cout << "  2. 两直线交点\n";
        cout << "  3. 点到直线距离\n";
        cout << "  4. 两点中垂线\n";
        cout << "  9. 两平行线距离\n";
        cout << " 10. 点关于直线对称\n";
        cout << " 11. 两直线夹角\n";
        cout << "  --- 圆 ---\n";
        cout << "  5. 三点求外接圆\n";
        cout << "  6. 直线与圆\n";
        cout << "  7. 两圆关系\n";
        cout << "  8. 点到圆的切线\n";
        cout << " 12. 直线截圆弦长\n";
        cout << " 16. 两圆公共弦\n";
        cout << " 17. 点到圆最远/最近距离\n";
        cout << "  --- 三角形 ---\n";
        cout << " 13. 三角形面积\n";
        cout << " 14. 三角形重心\n";
        cout << " 15. 三角形垂心\n";
        cout << "  0. 退出\n";
        cout << "请选择: ";
        int c; cin >> c;
        if (c == 0) break;

        if (c == 1) {
            Pt p, q; ip(p, "P1"); ip(q, "P2");
            Ln L = linePP(p, q);
            cout << "直线: " << L.A << "x + " << L.B << "y + " << L.C << " = 0\n";
        }
        else if (c == 2) {
            Ln L, M; il(L, "L1"); il(M, "L2");
            Pt p = interLL(L, M);
            if (isnan(p.x)) cout << "平行\n";
            else cout << "交点: (" << p.x << ", " << p.y << ")\n";
        }
        else if (c == 3) {
            Pt p; ip(p, "P");
            Ln L; il(L, "直线");
            cout << "距离: " << distPL(p, L) << "\n";
        }
        else if (c == 4) {
            Pt p, q; ip(p, "P1"); ip(q, "P2");
            Ln L = perpBisector(p, q);
            cout << "中垂线: " << L.A << "x + " << L.B << "y + " << L.C << " = 0\n";
        }
        else if (c == 5) {
            Pt A, B, C; ip(A, "A"); ip(B, "B"); ip(C, "C");
            Cir res;
            if (!circumcircle(A, B, C, res)) cout << "三点共线\n";
            else cout << "圆心: (" << res.O.x << ", " << res.O.y
                      << ")  半径: " << res.r << "\n";
        }
        else if (c == 6) {
            Pt p, q; Cir C;
            ip(p, "直线上点P1"); ip(q, "直线上点P2");
            Ln L = linePP(p, q);
            ic(C, "圆");
            double d = distPL(C.O, L);
            cout << "d = " << d << ", r = " << C.r << "\n";
            if (d > C.r + 1e-9) cout << "相离\n";
            else if (fabs(d - C.r) < 1e-9) {
                Pt pts[2]; interLC(L, C, pts);
                cout << "相切, 切点: (" << pts[0].x << ", " << pts[0].y << ")\n";
            } else {
                Pt pts[2]; interLC(L, C, pts);
                cout << "相交\n";
                cout << "交点1: (" << pts[0].x << ", " << pts[0].y << ")\n";
                cout << "交点2: (" << pts[1].x << ", " << pts[1].y << ")\n";
            }
        }
        else if (c == 7) {
            Cir C1, C2; ic(C1, "圆1"); ic(C2, "圆2");
            cout << "位置关系: " << circleRelation(C1, C2) << "\n";
            Pt pts[2];
            int cnt = interCC(C1, C2, pts);
            for (int i = 0; i < cnt; i++)
                cout << "交点" << i + 1 << ": (" << pts[i].x << ", " << pts[i].y << ")\n";
        }
        else if (c == 8) {
            Pt P; Cir C; ip(P, "P"); ic(C, "圆");
            Ln tans[2];
            int cnt = tangentsFromPoint(P, C, tans);
            if (cnt == 0) cout << "点在圆内，无切线\n";
            else for (int i = 0; i < cnt; i++)
                cout << "切线" << (cnt > 1 ? to_string(i + 1) : "") << ": "
                     << tans[i].A << "x + " << tans[i].B << "y + "
                     << tans[i].C << " = 0\n";
        }
        else if (c == 9) {
            Ln L, M; il(L, "L1"); il(M, "L2");
            if (fabs(L.A * M.B - L.B * M.A) > 1e-9) {
                cout << "两直线不平行，距离无意义\n";
            } else {
                cout << "两平行线距离: " << distLL(L, M) << "\n";
            }
        }
        else if (c == 10) {
            Pt p; ip(p, "P"); Ln L; il(L, "直线");
            Pt Q = reflectPt(p, L);
            cout << "对称点: (" << Q.x << ", " << Q.y << ")\n";
        }
        else if (c == 11) {
            Ln L, M; il(L, "L1"); il(M, "L2");
            cout << "夹角: " << angleLL(L, M) << " 度\n";
        }
        else if (c == 12) {
            Pt p, q; Cir C;
            ip(p, "直线上点P1"); ip(q, "直线上点P2");
            Ln L = linePP(p, q);
            ic(C, "圆");
            double ch = chordLC(L, C);
            if (ch < 0) cout << "相离，无弦\n";
            else cout << "弦长: " << ch << "\n";
        }
        else if (c == 13) {
            Pt A, B, C; ip(A, "A"); ip(B, "B"); ip(C, "C");
            cout << "面积: " << triArea(A, B, C) << "\n";
        }
        else if (c == 14) {
            Pt A, B, C; ip(A, "A"); ip(B, "B"); ip(C, "C");
            Pt G = centroid(A, B, C);
            cout << "重心: (" << G.x << ", " << G.y << ")\n";
        }
        else if (c == 15) {
            Pt A, B, C; ip(A, "A"); ip(B, "B"); ip(C, "C");
            Pt H = orthocenter(A, B, C);
            if (isnan(H.x)) cout << "三点共线\n";
            else cout << "垂心: (" << H.x << ", " << H.y << ")\n";
        }
        else if (c == 16) {
            Cir C1, C2; ic(C1, "圆1"); ic(C2, "圆2");
            string rel = circleRelation(C1, C2);
            cout << "位置关系: " << rel << "\n";
            if (rel == "相交" || rel == "外切" || rel == "内切") {
                Ln L = commonChord(C1, C2);
                cout << "公共弦: " << L.A << "x + " << L.B << "y + "
                     << L.C << " = 0\n";
            } else {
                cout << "无公共弦\n";
            }
        }
        else if (c == 17) {
            Pt P; Cir C; ip(P, "P"); ic(C, "圆");
            double nd, fd;
            distPR(P, C, nd, fd);
            cout << "最近距离: " << nd << "\n";
            cout << "最远距离: " << fd << "\n";
        }
        else {
            cout << "无效选项\n";
        }
    }
    cout << "再见！\n";
    return 0;
}
