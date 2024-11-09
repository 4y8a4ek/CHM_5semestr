#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <functional>

enum Wavelet
{
    D6,
    SHANNON
};

class DWT
{
private:
    const std::vector<double> D6_low = {0.3326705529500826, 0.8068915093110928, 0.4598775021184915, -0.1350110200102546, -0.0854412738822415, 0.0352262918857095};
    const std::vector<double> D6_high = {-0.0352262918857095, -0.0854412738822415, 0.1350110200102546, 0.4598775021184915, -0.8068915093110928, 0.3326705529500826};
    const std::vector<double> shannon_low = {0.3, 0.6, 0.9, 1.0, 0.9, 0.6, 0.3};
    const std::vector<double> shannon_high = {-0.3, -0.6, 0.9, -1.0, 0.9, -0.6, 0.3};
    std::vector<double> signals, input_signals, approx, details, restored_signals;
    int size;
public:
    DWT(const std::vector<double> &input);
    void decompose(Wavelet type);
    void restore(Wavelet type);
    void save(const std::vector<double> &vec, const std::string &filename);
    double computeError();
    const std::vector<double>& getApproximation() const;
    const std::vector<double>& getDetail() const;
    const std::vector<double>& getRestoredSignal() const;
    void setSignals(const std::vector<double> &input);
};
