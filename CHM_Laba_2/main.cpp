#include <iostream>
#include <vector>
#include <cmath>

// Структура для хранения точки
struct Point {
    double x, y;
};

// Функция для вычисления линейного сплайна между двумя точками
double linearSpline(double x, const Point& p1, const Point& p2) {
    double slope = (p2.y - p1.y) / (p2.x - p1.x);
    return p1.y + slope * (x - p1.x);
}

// Функция для вычисления аппроксимации по методу наименьших квадратов
double leastSquaresApproximation(const std::vector<Point>& points, double x) {
    // Поиск отрезка, где находится x
    for (size_t i = 0; i < points.size() - 1; ++i) {
        if (x >= points[i].x && x <= points[i + 1].x) {
            return linearSpline(x, points[i], points[i + 1]);
        }
    }
    return 0; // Вне диапазона
}

// Функция для вычисления нормы L2 (интеграл ошибки)
double calculateL2Norm(const std::vector<Point>& points) {
    double errorSum = 0.0;
    for (size_t i = 0; i < points.size() - 1; ++i) {
        // Разбиваем отрезок на n маленьких участков и численно интегрируем
        const int n = 1000000; // Количество шагов для интегрирования
        double dx = (points[i + 1].x - points[i].x) / n;
        for (int j = 0; j < n; ++j) {
            double x1 = points[i].x + j * dx;
            double x2 = points[i].x + (j + 1) * dx;
            double f1 = leastSquaresApproximation(points, x1);
            double f2 = leastSquaresApproximation(points, x2);
            double actual1 = points[i].y;  // Табличное значение
            double actual2 = points[i + 1].y; // Табличное значение
            errorSum += ((f1 - actual1) * (f1 - actual1) + (f2 - actual2) * (f2 - actual2)) * dx / 2.0;
        }
    }
    return std::sqrt(errorSum);
}

int main() {
    // Пример табличных данных (x, f(x))
    std::vector<Point> points = {
        {0, 0},
        {1, 2},
        {2, 6}
    };

    // Аппроксимация в точке
    double x;
    std::cin >> x;
    double approxValue = leastSquaresApproximation(points, x);
    std::cout << "Approximated value at x = " << x << " : " << approxValue << std::endl;

    // Вычисление ошибки
    double l2Norm = calculateL2Norm(points);
    std::cout << "L2 Norm of the error: " << l2Norm << std::endl;

    return 0;
}
