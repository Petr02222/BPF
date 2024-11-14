#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

class BPF
{
    std::vector<std::complex<double>> data_;
public:


    BPF(std::vector<std::complex<double>>& data)//конструктор для инициализации свойств, который добивает до степени двойки
    {
        int N = data.size();
        int m = 1;
        while (m < N) m <<= 1;
        data_.resize(m, 0);
        for (int i = 0; i < N; i++)
        {
            data_[i] = data[i];
        }
    }

    void fft(std::vector<std::complex<double>>& x)//Функция прямого преобразования фурье через рекурсивные вызовы
    {
        int N = x.size();
        if (N <= 1) return;

        std::vector<std::complex<double>> chet(N / 2);
        std::vector<std::complex<double>> nechet(N / 2);

        for (int i = 0; i < N / 2; i++)//по формулам
        {
            chet[i] = x[i * 2];
            nechet[i] = x[i * 2 + 1];
        }

        fft(chet);
        fft(nechet);

        for (int i = 0; i < N / 2; i++)
        {
            double pex = -2 * M_PI * i / N;
            std::complex<double> t = std::polar(1.0, pex) * nechet[i];
            x[i] = chet[i] + t;
            x[i + N / 2] = chet[i] - t;
        }
    }

    void ifft(std::vector<std::complex<double>>& x)//обратное фурье
    {
        double N = x.size();
        for (int i = 0; i < N; i++)
        {
            x[i] = std::conj(x[i]);//меняем знак у мнимой части
        }
        fft(x);//просиходит обычное преобразование фурье
        for (int i = 0; i < N; i++)
        {
            x[i] = std::conj(x[i]) / N;//опять меняем знак у мнимой части 
        }
    }

    void transform(bool k)//метод преобразования
    {
        if (k == 1) fft(data_);//если установленная 1, то делается прямое
        else ifft(data_);//если установлен 0, то обратное 
    }
    std::vector<std::complex<double>> getResult()//метод вывода
    {
        return data_;
    }
};

int main()
{
    setlocale(LC_CTYPE, "Ru");
    std::cout << "Введите размер последовательности: ";
    int N;
    std::cin >> N;
    srand(time(0));
    std::complex<double> alfa = 0;
    std::vector<std::complex<double>> data(N);

    for (int i = 0; i < data.size(); i++)
    {
        int real = rand() % 21 - 10;
        int imag = rand() % 21 - 10;
        data[i] = std::complex<int>(real, imag);
    }

    BPF fft(data);
    fft.transform(1);//прямое

    std::vector<std::complex<double>> result = fft.getResult();
    std::cout << "Результат преобразования фурье: " << std::endl;
    for (int i = 0; i < data.size(); i++)
    {
        std::cout << result[i] << std::endl;
    }
    std::cout << std::endl;
    fft.transform(0);//обратное
    std::vector<std::complex<double>> inverse_result = fft.getResult();
    std::cout << "Результат обратного преобразования фурье: " << std::endl;
    for (int i = 0; i < data.size(); i++)
    {
        std::cout << inverse_result[i] << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < data.size(); i++)
    {
        alfa += (data[i] - inverse_result[i]) * (data[i] - inverse_result[i]);
    }
    std::cout << "Средняя квадратичная ошибка " << alfa << std::endl;
}