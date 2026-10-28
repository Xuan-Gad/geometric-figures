#include <bits/stdc++.h>
using namespace std;

struct Pt { double x, y; };
struct Ln { double A, B, C; };

Ln linePP(Pt p, Pt q) {
    return {p.y - q.y, q.x - p.x, (q.y - p.y) * p.x + (p.x - q.x) * p.y};
}

Pt inter(Ln L, Ln M) {
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

void ip(Pt &p, string n) {
    cout << n << " x y: ";
    cin >> p.x >> p.y;
}

int main() {
    cout << "===== 解析几何求解器 =====\n";
    for (;;) {
        cout << "\n========== 菜单 ==========\n";
        cout << "  1. 两点求直线\n";
        cout << "  2. 两直线交点\n";
        cout << "  3. 点到直线距离\n";
        cout << "  4. 两点中垂线\n";
        cout << "  0. 退出\n";
        cout << "==========================\n";
        cout << "请选择: ";
        int c;
        cin >> c;

        if (c == 0) break;

        if (c == 1) {
            Pt p, q;
            ip(p, "P1");
            ip(q, "P2");
            Ln L = linePP(p, q);
            cout << "直线: " << L.A << "x + " << L.B << "y + " << L.C << " = 0\n";
        }
        else if (c == 2) {
            double a1, b1, c1, a2, b2, c2;
            cout << "请按 Ax+By+C=0 的形式输入系数\n";
            cout << "L1 A1 B1 C1: "; cin >> a1 >> b1 >> c1;
            cout << "L2 A2 B2 C2: "; cin >> a2 >> b2 >> c2;
            Pt p = inter({a1, b1, c1}, {a2, b2, c2});
            if (isnan(p.x)) cout << "平行（或无交点）\n";
            else cout << "交点: (" << p.x << ", " << p.y << ")\n";
        }
        else if (c == 3) {
            Pt p;
            ip(p, "P");
            double a, b, c;
            cout << "直线 A B C: "; cin >> a >> b >> c;
            cout << "距离: " << distPL(p, {a, b, c}) << "\n";
        }
        else if (c == 4) {
            Pt p, q;
            ip(p, "P1");
            ip(q, "P2");
            Ln L = perpBisector(p, q);
            cout << "中垂线: " << L.A << "x + " << L.B << "y + " << L.C << " = 0\n";
        }
        else {
            cout << "无效选项，请重新输入\n";
        }
    }
    cout << "再见！\n";
    return 0;
}
