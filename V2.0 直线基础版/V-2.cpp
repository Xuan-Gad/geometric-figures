#include <bits/stdc++.h>
using namespace std;

struct Frac {
    long long n, d;
    static long long gcd(long long a, long long b) {
        a = a < 0 ? -a : a;
        b = b < 0 ? -b : b;
        while (b) { long long t = a % b; a = b; b = t; }
        return a;
    }
    Frac(long long num = 0, long long den = 1) : n(num), d(den) {
        if (d == 0) { n = 0; d = 1; return; }
        if (d < 0) { n = -n; d = -d; }
        long long g = gcd(n < 0 ? -n : n, d);
        n /= g; d /= g;
    }
    Frac operator+(const Frac& o) const { return {n * o.d + o.n * d, d * o.d}; }
    Frac operator-(const Frac& o) const { return {n * o.d - o.n * d, d * o.d}; }
    Frac operator*(const Frac& o) const { return {n * o.n, d * o.d}; }
    Frac operator/(const Frac& o) const { return {n * o.d, d * o.n}; }
    Frac operator-() const { return {-n, d}; }
    bool operator==(const Frac& o) const { return n == o.n && d == o.d; }
    friend ostream& operator<<(ostream& os, const Frac& f) {
        if (f.d == 1) return os << f.n;
        return os << f.n << '/' << f.d;
    }
};

istream& operator>>(istream& is, Frac& f) {
    long long a, b;
    char c;
    if (is >> a && is >> c && c == '/' && is >> b)
        f = Frac(a, b);
    else
        is.setstate(ios::failbit);
    return is;
}

struct Pt { Frac x, y; };
struct Ln { Frac A, B, C; };

Ln linePP(Pt p, Pt q) {
    Frac A = p.y - q.y;
    Frac B = q.x - p.x;
    Frac C = (q.y - p.y) * p.x + (p.x - q.x) * p.y;
    return {A, B, C};
}

optional<Pt> inter(Ln L, Ln M) {
    Frac D = L.A * M.B - L.B * M.A;
    if (D.n == 0) return {};
    Frac x = (L.B * M.C - L.C * M.B) / D;
    Frac y = (L.C * M.A - L.A * M.C) / D;
    return Pt{x, y};
}

void distPL(Pt p, Ln L) {
    Frac num = L.A * p.x + L.B * p.y + L.C;
    if (num.n < 0) num.n = -num.n;
    long long a2 = L.A.n * L.A.n, b2 = L.B.n * L.B.n;
    long long denSq = a2 * L.B.d * L.B.d + b2 * L.A.d * L.A.d;
    if (denSq == 0) { cout << "0\n"; return; }
    double val = (double)num.n / num.d / sqrt((double)denSq);
    cout << num << " / sqrt(" << denSq << ")";
    if (num.d != 1) cout << "  (约 " << fixed << setprecision(4) << val << ")";
    cout << "\n";
}

Ln perpBisector(Pt p, Pt q) {
    Pt m{(p.x + q.x) / Frac(2), (p.y + q.y) / Frac(2)};
    Frac A = q.x - p.x;
    Frac B = q.y - p.y;
    Frac C = -(A * m.x + B * m.y);
    return {A, B, C};
}

void ip(Pt& p, string n) {
    cout << n << " x y: ";
    cin >> p.x >> p.y;
}

int main() {
    cout << "===== 解析几何 v2.0（分数版）=====\n";
    for (;;) {
        cout << "\n1.两点求直线 2.交点 3.点到直线距离 4.中垂线 0.退出\n> ";
        int c; cin >> c;
        if (c == 0) break;
        if (c == 1) {
            Pt p, q; ip(p, "P1"); ip(q, "P2");
            Ln L = linePP(p, q);
            cout << "直线: " << L.A << "x + " << L.B << "y + " << L.C << " = 0\n";
        } else if (c == 2) {
            Frac a1, b1, c1, a2, b2, c2;
            cout << "L1 A B C: "; cin >> a1 >> b1 >> c1;
            cout << "L2 A B C: "; cin >> a2 >> b2 >> c2;
            auto res = inter({a1, b1, c1}, {a2, b2, c2});
            if (!res) cout << "平行\n";
            else cout << "交点: (" << res->x << ", " << res->y << ")\n";
        } else if (c == 3) {
            Pt p; ip(p, "P");
            Frac a, b, c;
            cout << "A B C: "; cin >> a >> b >> c;
            distPL(p, {a, b, c});
        } else if (c == 4) {
            Pt p, q; ip(p, "P1"); ip(q, "P2");
            Ln L = perpBisector(p, q);
            cout << "中垂线: " << L.A << "x + " << L.B << "y + " << L.C << " = 0\n";
        }
    }
}
