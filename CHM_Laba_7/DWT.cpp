#include "DWT.h"
DWT::DWT(const std::vector<double> &input) : signals(input), input_signals(input)
{
    size = input.size();
    approx.resize(size / 2);
    details.resize(size / 2);
    restored_signals.resize(size * 2);
}
void DWT::decompose(Wavelet type)
{
    const std::vector<double> *low, *high;
    if (type == D6)
    {
        low = &D6_low;
        high = &D6_high;
    }
    else
    {
        low = &shannon_low;
        high = &shannon_high;
    }

    approx.resize(size / 2);
    details.resize(size / 2);

    for (int i = 0; i < size / 2; ++i)
    {
        double approxSum = 0.0;
        double detailSum = 0.0;
        for (int j = 0; j < low->size(); ++j)
        {
            int index = (2 * i + j) % size;
            approxSum += signals[index] * (*low)[j];
            detailSum += signals[index] * (*high)[j];
        }
        approx[i] = approxSum;
        details[i] = detailSum;
    }
}
void DWT::restore(Wavelet type)
{
    const std::vector<double> *low, *high;
    if (type == D6)
    {
        low = &D6_low;
        high = &D6_high;
    }
    else
    {
        low = &shannon_low;
        high = &shannon_high;
    }

    restored_signals.resize(size * 2);
    fill(restored_signals.begin(), restored_signals.end(), 0.0);

    for (int i = 0; i < size / 2; ++i)
    {
        for (int j = 0; j < low->size(); ++j)
        {
            int index = (2 * i + j) % (size * 2);
            restored_signals[index] += approx[i] * (*low)[j];
            restored_signals[index] += details[i] * (*high)[j];
        }
    }
}
void DWT::save(const std::vector<double> &vec, const std::string &filename)
{
    std::ofstream file(filename);
    for (double v : vec)
    {
        file << std::fixed << std::setprecision(6) << v << "\n";
    }
    file.close();
}
double DWT::computeError()
{
    double error = 0.0;
    for (int i = 0; i < size; ++i)
    {
        error += pow(input_signals[i] - restored_signals[i], 2);
    }
    return sqrt(error / size);
}
const std::vector<double> &DWT::getApproximation() const { return approx; }

const std::vector<double> &DWT::getDetail() const { return details; }

const std::vector<double> &DWT::getRestoredSignal() const { return restored_signals; }

void DWT::setSignals(const std::vector<double> &input)
{
    signals = input;
    size = signals.size();
}