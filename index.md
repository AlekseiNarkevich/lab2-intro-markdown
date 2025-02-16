# Первый
## Второй
### Третий
#### Четвёртый
##### Пятый
###### Шестой

# Картинка
# ![Image of Yaktocat](https://octodex.github.com/images/yaktocat.png)

## Л.р. по LDLT разложение матрицы.
```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <iomanip>
using namespace std;
using namespace std::chrono;
// Функция для вывода первых и последних нескольких элементов вектора
void printVector(const vector<double>& v, int limit = 5)
{
    int n = v.size();
    for (int i = 0; i < limit; ++i)
    {
        cout << v[i] << " ";
    }
    cout << "... ";
    for (int i = n - limit; i < n; ++i)
    {
        cout << v[i] << " ";
    }
    cout << endl;
}
// Функция для умножения матрицы на вектор
vector<double> multiply(const vector<vector<double>>& A, const vector<double>& x)
{
    int n = A.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}
// Функция для решения системы уравнений с использованием LDLT разложения
vector<double> ldltSolve(vector<vector<double>>& A, vector<double>& b)
{
    int n = A.size();
    vector<double> x(n);
    vector<double> y(n);
    vector<double> d(n); // Диагональ
    // Разложение A = LDL^T
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < j; ++i)
        {
            A[j][j] -= A[j][i] * A[j][i] * d[i];  // Вычисление диагонального элемента
        }
        d[j] = A[j][j];  // Диагональные элементы
        for (int i = j + 1; i < n; ++i)
        {
            for (int k = 0; k < j; ++k)
            {
                A[i][j] -= A[i][k] * A[j][k] * d[k];  // Вычисление нижнетреугольных элементов
            }
            A[i][j] /= d[j];
        }
    }
    // Прямой ход: решаем Ly = b
    for (int i = 0; i < n; ++i)
    {
        y[i] = b[i];
        for (int j = 0; j < i; ++j)
        {
            y[i] -= A[i][j] * y[j];
        }
    }
    // Решение системы D * z = y (D — диагональная)
    for (int i = 0; i < n; ++i)
    {
        y[i] /= d[i];
    }
    // Обратный ход: решаем L^T * x = z
    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j)
        {
            x[i] -= A[j][i] * x[j];
        }
    }
    return x;
}
int main()
{
    setlocale(LC_ALL, "ru");
    srand(time(0));
    auto start = high_resolution_clock::now();
    int n = 1000;
    int m = 6;
    vector<vector<double>> A(n, vector<double>(n));
    // Генерация недиагональных элементов и создание симметричной матрицы
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < i; ++j) 
        {
            A[i][j] = -100 + rand() % 201;  // Случайное число от -100 до 100
            A[j][i] = A[i][j];  // Симметричность
        }
    }
    // Генерация диагональных элементов
    for (int i = 0; i < n; ++i) 
    {
        double sum_abs = 0.0;
        for (int j = 0; j < n; ++j)
        {
            if (i != j) 
            {
                sum_abs += fabs(A[i][j]);  // Сумма модулей недиагональных элементов
            }
        }
        // Диагональный элемент выбирается случайным образом из указанного диапазона
        double lower_bound = sum_abs + m;
        double upper_bound = sum_abs + 10 * m;
        A[i][i] = lower_bound + static_cast<double>(rand()) / RAND_MAX * (upper_bound - lower_bound);
    }
    // Точное решение x_exact (вектор от m до n+m)
    vector<double> x_exact(n);
    for (int i = 0; i < n; ++i)
    {
        x_exact[i] = m + i;
    }
    // Вычисляем правую часть b = A * x_exact
    vector<double> b = multiply(A, x_exact);
    // Решаем систему уравнений с использованием LDLT разложения
    vector<double> x_approx = ldltSolve(A, b);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end - start);
    cout << fixed << setprecision(16);
    cout << "Первые и последние 5 координат точного решения x:" << endl << endl;
    printVector(x_exact);
    cout << endl;
    cout << "Первые и последние 5 координат приближенного решения x:" << endl << endl;
    printVector(x_approx);
    cout << endl;
    double norm_exact = 0.0;
    double norm_diff = 0.0;
    for (int i = 0; i < n; ++i)
    {
        norm_exact = max(norm_exact, fabs(x_exact[i]));
        norm_diff = max(norm_diff, fabs(x_exact[i] - x_approx[i]));
    }
    double relative_error = (norm_diff / norm_exact) * 100;
    cout << scientific << "Относительная погрешность: " << relative_error << "%" << endl << endl;
    cout << fixed << "Время выполнения: " << duration.count() << " секунд" << endl << endl;
    return 0;
}
```
