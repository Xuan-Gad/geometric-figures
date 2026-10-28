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

// [修复] 同时支持 "3"、"1/2"、"1 / 2"、"1 /2" 等写法
istream& operator>>(istream& is, Frac& f) {
    long long a, b;
    if (!(is >> a)) return is;
    // 跳过 a 和 '/' 之间的空格/制表符
    while (is.peek() == ' ' || is.peek() == '\t') is.get();
    if (is.peek() == '/') {
        is.get();                        // 吃掉 '/'
        if (is >> b) { f = Frac(a, b); return is; }
        is.setstate(ios::failbit);
        return is;
    }
    f = Frac(a, 1);
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

// 返回 false 表示平行/重合，无唯一交点
bool inter(Ln L, Ln M, Pt& out) {
    Frac D = L.A * M.B - L.B * M.A;
    if (D.n == 0) return false;
    out.x = (L.B * M.C - L.C * M.B) / D;
    out.y = (L.C * M.A - L.A * M.C) / D;
    return true;
}

// [修复] 补上 A.d * B.d 这个因子
void distPL(Pt p, Ln L) {
    Frac num = L.A * p.x + L.B * p.y + L.C;
    if (num.n < 0) num.n = -num.n;

    // denSq = A.n2·B.d2 + B.n2·A.d2，使 sqrt(A2+B2) = sqrt(denSq)/(A.d·B.d)
    long long denSq = (L.A.n * L.A.n) * (L.B.d * L.B.d)
                    + (L.B.n * L.B.n) * (L.A.d * L.A.d);
    if (denSq == 0) { cout << "0\n"; return; }   // 退化直线

    long long scale = L.A.d * L.B.d;
    long long numer = num.n * scale;    // |num| * A.d * B.d
    long long denom = num.d;
    long long g = Frac::gcd(numer, denom);
    numer /= g; denom /= g;

    double val = (double)numer / denom / sqrt((double)denSq);

    cout << numer;
    if (denom != 1) cout << "/" << denom;
    cout << " / sqrt(" << denSq << ")";
    if (denom != 1)
        cout << "  (约 " << fixed << setprecision(4) << val << ")";
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
        }  else if (c == 2) {
    Frac a1, b1, c1, a2, b2, c2;
    cout << "L1 A B C: "; cin >> a1 >> b1 >> c1;
    cout << "L2 A B C: "; cin >> a2 >> b2 >> c2;
    Pt res;
    if (!inter({a1, b1, c1}, {a2, b2, c2}, res))
        cout << "平行\n";
    else
        cout << "交点: (" << res.x << ", " << res.y << ")\n";
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
