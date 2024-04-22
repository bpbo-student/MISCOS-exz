#include <iostream>
#include <locale>

// Константа для точности вычислений
const double eps = 1e-15; // Точность вычислений

// Коэффициенты фильтра
double b0, b1, b2, a1, a2;

// Класс для представления комплексных чисел
class ComplexNumber {
private:
    double real;      // Действительная часть комплексного числа
    double imaginary; // Мнимая часть комплексного числа

public:
    // Конструктор для инициализации комплексного числа
    ComplexNumber(double realPart, double imaginaryPart) : real(realPart), imaginary(imaginaryPart) {}

    // Методы для получения действительной и мнимой частей комплексного числа
    double getReal() const { return real; }
    double getImaginary() const { return imaginary; }

    // Перегрузка оператора сложения для комплексных чисел
    ComplexNumber operator+(const ComplexNumber& other) const {
        double sumReal = real + other.real;
        double sumImaginary = imaginary + other.imaginary;
        return ComplexNumber(sumReal, sumImaginary);
    }

    // Перегрузка оператора умножения для комплексных чисел
    ComplexNumber operator*(const ComplexNumber& other) const {
        double realPart = real * other.real - imaginary * other.imaginary;
        double imaginaryPart = real * other.imaginary + imaginary * other.real;
        return ComplexNumber(realPart, imaginaryPart);
    }
};

// Функция для вычисления абсолютного значения числа
double dmAbs(double x) {
    return (x < 0) ? -x : x;
}

// Функция для вычисления квадратного корня числа
double dmSqrt(double x) {
    if (x < 0) {
        return 0.0; // Квадратный корень от отрицательного числа возвращает 0
    }

    // Начальное предположение
    double guess = x;
    double newGuess;
    int maxIterations = 1000; // Максимальное количество итераций для метода Ньютона

    // Итеративный метод Ньютона для вычисления квадратного корня
    while (maxIterations--) {
        newGuess = 0.5 * (guess + x / guess); // Улучшение предположения
        if (dmAbs(newGuess - guess) < eps) {
            return newGuess; // Возвращение при достижении необходимой точности
        }
        guess = newGuess; // Обновление предположения для следующей итерации
    }

    // Возвращение последнего предположения после превышения максимального числа итераций
    return guess;
}

// Функция для решения квадратного уравнения и вывода его корней
void solveQuadraticEquation(double a, double b, double c) {
    double discriminant = b * b - 4 * a * c; // Вычисление дискриминанта

    // Если дискриминант больше нуля, выводим два действительных корня
    if (discriminant > 0) {
        double root1 = (-b + dmSqrt(discriminant)) / (2 * a);
        double root2 = (-b - dmSqrt(discriminant)) / (2 * a);
        std::cout << root1 << " и " << root2 << std::endl;
    }
    // Если дискриминант равен нулю, выводим один действительный корень
    else if (discriminant == 0) {
        double root = -b / (2 * a);
        std::cout << root << std::endl;
    }
    // Если дискриминант отрицательный, выводим два комплексных корня
    else {
        double realPart = -b / (2 * a);
        double imaginaryPart = dmSqrt(-discriminant) / (2 * a);
        ComplexNumber root1(realPart, imaginaryPart);
        ComplexNumber root2(realPart, imaginaryPart);
        std::cout << root1.getReal() << " + " << root1.getImaginary() << "i"
            << " и " << root2.getReal() << " - " << root2.getImaginary() << "i" << std::endl;
    }
}

// Функция для ввода коэффициентов фильтра
void get_coef() {
    std::cout << "Введите коэффициенты фильтра (со знаком .)" << std::endl;

    std::cout << "b(0) = ";
    std::cin >> b0;

    std::cout << "b(1) = ";
    std::cin >> b1;

    std::cout << "b(2) = ";
    std::cout << "a(1) = ";
    std::cin >> a1;

    std::cout << "a(2) = ";
    std::cin >> a2;

    std::cout << "\n";
}

// Функция для вывода функции H(z) фильтра
void Hz() {
    std::cout << "Функция H(z):" << std::endl;
    std::cout << "         " << b0 << " + " << b1 << " + " << b2 << std::endl;
    std::cout << "H(z) = ———————————————————" << std::endl;
    std::cout << "         " << 1 << " - (" << a1 << " + " << a2 << ") " << std::endl;
}

// Функция для вывода функции y(n) фильтра
void yn() {
    std::cout << "Функция y(n):" << std::endl;
    std::cout << "y(n) = " << b0 << "x(n) + " << b1 << "x(n-1) + " << b2 << "x(n-2) + " << a1 << "y(n-1) + " << a2 << "y(n-2)" << std::endl;

    std::cout << "\n";
}

// Функция для вывода нулей и полюсов фильтра
void zeroplot() {
    std::cout << "Нули: ";
    solveQuadraticEquation(b0, b1, b2);
    std::cout << "Полюса: ";
    solveQuadraticEquation(1, -a1, -a2);

    std::cout << "\n";
}

// Функция для определения импульсной характеристики фильтра
void ih() {
    std::cout << "Импульсная характеристика:" << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::cout << "y(" << i << ") = (" << b0 << " * " << xi(i) << ") + (" << b1 << " * " << xi(i - 1) << ") + (" << b2 << " * " << xi(i - 3) << ") + (" << a1 << " * " << yi(i - 1, b0, b1, b2, a1, a2) << ") + (" << a2 << " * " << yi(i - 2, b0, b1, b2, a1, a2) << ")" << std::endl;
        std::cout << "y(" << i << ") = " << yi(i, b0, b1, b2, a1, a2) << std::endl;
    }
    std::cout << "\n";
}

// Функция для определения полной характеристики фильтра
void ph() {
    std::cout << "Полная характеристика:" << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::cout << "y(" << i << ") = (" << b0 << " * " << xp(i) << ") + (" << b1 << " * " << xp(i - 1) << ") + (" << b2 << " * " << xp(i - 3) << ") + (" << a1 << " * " << yp(i - 1, b0, b1, b2, a1, a2) << ") + (" << a2 << " * " << yp(i - 2, b0, b1, b2, a1, a2) << ")" << std::endl;
        std::cout << "y(" << i << ") = " << yp(i, b0, b1, b2, a1, a2) << std::endl;
    }
    std::cout << "\n";
}

// Главная функция программы
int main(int argc, char* argv[]) {
    // Установка локали на русский язык
    setlocale(LC_ALL, "Russian");

    // Ввод коэффициентов фильтра
    get_coef();

    // Вывод функций фильтра
    Hz();
    yn();

    // Вывод нулей и полюсов фильтра
    zeroplot();

    // Вывод импульсной и полной характеристики фильтра
    ih();
    ph();

    return 0;
}
