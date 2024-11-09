#include "DWT.h"
int main()
{
    std::vector<double> signal(512, 0.0);
    for (int n = 128; n < 256; ++n)
    {
        signal[n] = sin(pow(n - 128, 1.66) / 128.0);
    }

    for (int n = 384; n < 448; ++n)
    {
        signal[n] = sin(pow(n - 128, 2.89) / 128.0);
    }

    DWT transformer(signal);
    std::vector<double> currentSignal = signal;
    transformer.save(signal, "inputed_signal.txt");
    Wavelet type = SHANNON;
    for (int level = 1; level <= 4; ++level)
    {
        transformer.decompose(type);
        transformer.save(transformer.getApproximation(), "approximation_level_" + std::to_string(level) + ".txt");
        transformer.save(transformer.getDetail(), "detail_" + std::to_string(level) + ".txt");
        transformer.restore(type);
        transformer.save(transformer.getRestoredSignal(), "reconstructed_level_" + std::to_string(level) + ".txt");
        currentSignal = transformer.getApproximation();
        transformer.setSignals(currentSignal);
    }
    transformer.setSignals(currentSignal);
    transformer.decompose(type);
    transformer.restore(type);
    double error = transformer.computeError();
    std::cout << "Reconstruction Error: " << error << std::endl;
    return 0;
}
