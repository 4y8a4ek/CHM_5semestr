#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include <algorithm>
using namespace std;

class CubicSpline {
private:
    vector<double> x, y; // Узлы интерполяции и значения функции в узлах
    vector<double> a, b, c, d; // Коэффициенты для каждого сегмента
    vector<double> h; // Шаги между узлами
    int n; // Количество узлов

public:
    // Конструктор, принимающий узлы интерполяции и значения
    CubicSpline(const vector<double>& x_vals, const vector<double>& y_vals) {
        if (x_vals.size() != y_vals.size() || x_vals.size() < 2) {
            throw invalid_argument("Некорректные данные для интерполяции.");
        }

        n = x_vals.size() - 1; // Количество сегментов n = количество узлов - 1
        x = x_vals;
        y = y_vals;

        a = y_vals; // Коэффициенты a_i равны y_i
        b.resize(n);
        c.resize(n + 1);
        d.resize(n);
        h.resize(n);

        // Вычисление шагов между узлами
        for (int i = 0; i < n; ++i) {
            h[i] = x[i + 1] - x[i];
            if (h[i] == 0) {
                throw invalid_argument("Некорректные узлы интерполяции: совпадающие значения.");
            }
        }

        computeSpline();
    }

    // Вычисление коэффициентов сплайна
    void computeSpline() {
        vector<double> alpha(n), l(n + 1), mu(n), z(n + 1);

        // Вычисление вспомогательных коэффициентов
        for (int i = 1; i < n; ++i) {
            alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1]);
        }

        // Граничные условия для натурального сплайна
        l[0] = 1;
        mu[0] = 0;
        z[0] = 0;

        for (int i = 1; i < n; ++i) {
            l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n] = 1;
        z[n] = 0;
        c[n] = 0;

        // Вычисление коэффициентов b, c, d
        for (int j = n - 1; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
            d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
        }
    }


    // Метод интерполяции значения в произвольной точке x
    double interpolate(double x_val){
        if (x_val < x[0] || x_val > x[n]) {
            throw out_of_range("Значение выходит за пределы интервала интерполяции.");
        }

        int i = n - 1;
        for (int j = 0; j < n; ++j) {
            if (x_val >= x[j] && x_val <= x[j + 1]) {
                i = j;
                break;
            }
        }
        double dx = x_val - x[i];
        return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
    }

    double first_derivative(double x_val){
        if (x_val < x[0] || x_val > x[n]) {
            throw out_of_range("Значение выходит за пределы интервала интерполяции.");
        }
        int i = n - 1;
        for (int j = 0; j < n; ++j) {
            if (x_val >= x[j] && x_val <= x[j + 1]) {
                i = j;
                break;
            }
        }
        double dx = x_val - x[i];
        return b[i] + 2 * c[i] * dx + 3 * d[i] * dx * dx;
    }

    double second_derivative(double x_val){
        if (x_val < x[0] || x_val > x[n]) {
            throw out_of_range("Значение выходит за пределы интервала интерполяции.");
        }
        int i = n - 1;
        for (int j = 0; j < n; ++j) {
            if (x_val >= x[j] && x_val <= x[j + 1]) {
                i = j;
                break;
            }
        }
        double dx = x_val -x[i];
        return 2 * c[i] + 6 * d[i] * dx;
    }
};

vector<double> function(const vector<double>& x){
    vector<double> y;
    for(int i = 0; i < x.size(); i++){
        y.push_back(pow((x[i] - 1), 2));
    }
    return y;
}
double f(double x){
    return pow((x-1),2) ;
}
double f_f_der(double x){
    return 2*x - 2;
}
double f_s_der(double x){
    return 2;
}
// Функция для регулярного разбиения отрезка [a, b] на n частей
vector<double> regular_grid(double a, double b, int n) {
    vector<double> grid(n + 1);
    double step = (b - a) / n;
    
    for (int i = 0; i <= n; ++i) {
        grid[i] = a + i * step;
    }
    
    return grid;
}

// Функция для адаптивного разбиения отрезка [a, b] с коэффициентом r
vector<double> adaptive_grid(double a, double b, int n, double r) {
    vector<double> grid(n + 1);
    double h1 = (b - a) * (r - 1) / (pow(r, n) - 1);
    grid[0] = a;
    for (int i = 1; i <= n; ++i) {
        grid[i] = grid[i - 1] + h1 * pow(r, i - 1);
    }
    return grid;
}

// Функция для вывода сетки
void print_grid(const vector<double>& grid) {
    for (double x : grid) {
        cout << x << " ";
    }
    cout << endl;
}
vector<vector<double>> create_Matrix(CubicSpline spline, vector<double> x){
    vector<vector<double>> Matrix;
    for( int i = 0; i < 10; i++){
        vector<double> cur_string;
        cur_string.push_back(x[i]);
        cur_string.push_back(spline.interpolate(x[i]));
        cur_string.push_back(spline.first_derivative(x[i]));
        cur_string.push_back(spline.second_derivative(x[i]));
        Matrix.push_back(cur_string);
    }
    return Matrix;
}
void print_Matrix(const vector<vector<double>>& matrix){
    for(vector<double> string : matrix){
        for(double x : string){
            cout << x << " ";
        }
        cout << endl;
    }
    cout << endl;
}
void print_Matrix_exp(const vector<vector<double>>& matrix){
    for(vector<double> string : matrix){
        for(double x : string){
            printf("%e", x);
            cout << "   |   " ;
        }
        cout << endl;
    }
    cout << endl;
}

double random_number(double lower_border, double higher_border){
    static std::mt19937 generator(static_cast<unsigned int>(std::time(0)));
    static std::uniform_real_distribution<> dis(lower_border, higher_border);
	return dis(generator);
}

vector<double> generate_random_x(double a, double b, vector<vector<double>> x_vals){
    vector<double> x;
    while(x.size() != 10 ){
        double curr = random_number(a,b);
        for(int i = 0; i < 3; i++){
        for(double element : x_vals[i]){
            if(element == curr){
                continue;
            }
        }
        }
        x.push_back(curr);

    }
    return x;
}

vector<vector<double>> approximate(vector<double> x, CubicSpline spline){
    double curr_aprox1 = 0, curr_aprox2 = 0, curr_aprox3 = 0;
    vector<vector<double>> result;
    vector<double> string;
    for(double element : x){
        string.push_back(fabs((f(element)) - (spline.interpolate(element))));
        string.push_back(fabs((f_f_der(element)) - (spline.first_derivative(element))));
        result.push_back(string);
        string.clear();
    }
    return result;
}

vector<double> first_der(vector<double> dots){
    vector<double> result;
    double der1 = (f(dots[1]) - f(dots[0]))/(dots[1] - dots[0]);
    double der2 = (f(dots[2]) - f(dots[1]))/(dots[2] - dots[1]);
    double der3 = (f(dots[2]) - f(dots[1]))/(dots[2] - dots[1]);
    result.push_back(der1); result.push_back(der2); result.push_back(der3);
    return result;
}
int main() {
    double a, b, r;
    int n;
    vector<double> dots = {-1, 0, 1};
    cout << "Введите a (начало отрезка): ";
    cin >> a;
    
    cout << "Введите b (конец отрезка): ";
    cin >> b;
    
    cout << "Введите n (количество сегментов): ";
    cin >> n;
    vector<vector<double>> x_valsh;
    for(int i = 0; i < 3; i++){
        x_valsh.push_back(regular_grid(a, b, n * pow(2,i)));
    }
    vector<double> x = generate_random_x(a + fabs((b-a)), b - fabs((b-a)), x_valsh);
    cout << " Матрица для пункта 2:" << endl;
    for(int i = 0; i < 10; i++){
        cout << x[i] << "   |   " << f(x[i]) << "   |   " << f_f_der(x[i]) << "   |   " << f_s_der(x[i]) << endl;
    }
    cout << endl;
    for(int i = 0; i < 3; i ++){
    vector<double> y_valsh = function(x_valsh[i]);
    CubicSpline splineH(x_valsh[i], y_valsh);
    vector<vector<double>> MatrixH;
    MatrixH = create_Matrix(splineH, x);
    cout << "Сетка с шагом = " << fabs((b-a)/(n * pow(2, i))) << endl << endl;
    for( int j = 0; j < x_valsh[i].size() - 1; j++){
        cout << x_valsh[i][j] << "   |    " << y_valsh[j]<< endl;
    }
    cout << endl;
    cout << "Матрица с шагом h:" << fabs((b-a)/(n * pow(2, i))) << endl;     
    print_Matrix_exp(MatrixH);

    cout << "Матрица аппроксимаций для шага h = " << fabs((b-a)/(n * pow(2,i))) << endl;
    print_Matrix(approximate(x, splineH));
    cout << endl;
    }

    vector<double> fd = first_der(dots);
    cout << "Значения первой производной вычисленные по формулам:" << endl;
    print_grid(fd);
    cout << endl;
    vector<double> xy_valsh = regular_grid(dots[0], dots[2], 4 * n);
    vector<double> y_valsh = function(xy_valsh);    
    CubicSpline spline(xy_valsh, y_valsh);
    cout << "Интерполированые значения первой производной:" << endl;
    for(int i = 0; i < 3; i ++){
        cout << spline.first_derivative(dots[i]) << "  ";
    }
    cout<< endl;
    return 0;
}
