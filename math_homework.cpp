#include <iostream>
#include <vector>
#include <windows.h>
#include <cmath>
#include <algorithm>

using namespace std;

vector<vector<double>> small_matrix(const vector<vector<double>>& m, int row, int col) {
    int n = m.size();
    vector<vector<double>> m2;

    for (int i=0; i<n; i++) {
        if (i == row) continue;
        vector<double> temp;
        for (int j=0; j<n; j++) {
            if (j == col) continue;
            temp.push_back(m[i][j]);
        }
        m2.push_back(temp);
    }
    return m2;
}

double det(const vector<vector<double>>& m) {
    int n = m.size();

    if (n == 1) return m[0][0];
    if (n == 2) return m[0][0]*m[1][1] - m[0][1]*m[1][0];

    double result = 0;
    for (int i=0; i<n; i++) {
        double sign = (i % 2 == 0) ? 1 : -1;
        result += sign * m[0][i] * det(small_matrix(m, 0, i));
    }
    return result;
}

bool inverse(const vector<vector<double>>& m, vector<vector<double>>& result) {
    int n = m.size();
    double d = det(m);

    if (fabs(d) < 1e-9) return false;

    result = vector<vector<double>> (n, vector<double>(n));

    if (n == 1) {
        result[0][0] = 1.0 / m[0][0];
        return true;
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            double sign = ((i + j) % 2 == 0) ? 1 : -1;
            double cofactor = sign * det(small_matrix(m, i, j));
            result[j][i] = cofactor / d;
        }
    }
    return true;
}

bool gauss(const vector<vector<double>>& m, vector<vector<double>>& result) {
    int n = m.size();

    vector<vector<double>> aug(n, vector<double>(2*n, 0.0));
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) aug[i][j] = m[i][j];
        aug[i][n+i] = 1.0;
    }

    for (int col=0; col<n; col++) {
        int pivot = col;
        for (int r=col+1; r<n; r++) {
            if (fabs(aug[r][col]) > fabs(aug[pivot][col])) pivot = r;
        }

        if (fabs(aug[pivot][col]) < 1e-9) return false;
        swap(aug[col], aug[pivot]);

        double p = aug[col][col];
        for (int j=0; j<2*n; j++) aug[col][j] /= p;

        for (int r=0; r<n; r++) {
            if (r == col) continue;
            double factor = aug[r][col];
            for (int j=0; j<2*n; j++) {
                aug[r][j] -= factor * aug[col][j];
            }
        }
    }

    result = vector<vector<double>> (n, vector<double>(n));
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            result[i][j] = aug[i][n+j];
        }
    }

    return true;
}

void print_matrix(const vector<vector<double>>& m) {
    for (const auto& row : m) {
        for (double v : row) cout << v << " ";
        cout << endl;
    }
}

bool same(const vector<vector<double>>& a, const vector<vector<double>>& b) {
    int n = a.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (fabs(a[i][j] - b[i][j]) > 1e-6) return false; 
        }
    }
    return true;
}

vector<vector<double>> multiply(const vector<vector<double>>& a, const vector<vector<double>>& b) {
    int n = a.size();
    vector<vector<double>> c(n, vector<double>(n, 0.0));

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int k=0; k<n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
            if (fabs(c[i][j]) < 1e-9) c[i][j] = 0;
        }
    }
    return c;
}

bool is_identity(const vector<vector<double>>& m) {
    int n = m.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (fabs(m[i][j] - expected) > 1e-6) return false;
        }
    }
    return true;
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    
    int n;
    cout << "정방행렬의 차수를 입력하시오 : ";
    cin >> n;

    vector<vector<double>> matrix(n, vector<double>(n));

    for (int i=0; i<n; i++) {
        cout << i+1 << "행 : ";
        for (int j=0; j<n; j++) {
            cin >> matrix[i][j];
        }
    }

    vector<vector<double>> inv1, inv2;
    bool ok1 = inverse(matrix, inv1);
    bool ok2 = gauss(matrix, inv2);

    if (!ok1 || !ok2) {
        cout << "\n오류: 행렬식이 0이므로 역행렬이 존재하지 않습니다." << endl;
        return 0;
    }

    cout << "\n행렬식으로 구한 역행렬:" << endl;
    print_matrix(inv1);
    
    cout << "\n가우스-조던 소거법으로 구한 역행렬:" << endl;
    print_matrix(inv2);

    if (same(inv1, inv2))
        cout << "\n두 방법의 결과가 동일합니다." << endl;
    else
        cout << "\n두 방법의 결과가 다릅니다." << endl;

    vector<vector<double>> check = multiply(matrix, inv1);
    cout << "\n검산 (A x 역행렬):" << endl;
    print_matrix(check);

    if (is_identity(check))
        cout << "\n검산 결과: 단위행렬이 됩니다. 역행렬이 맞습니다" << endl;
    else
        cout << "\n검산 결과: 단위행렬이 아닙니다." << endl;
    
    return 0;
}