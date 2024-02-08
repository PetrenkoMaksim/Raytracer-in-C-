#include <array>
#include <vector.h>

class Mat {
public:
    Mat() = default;
    Mat(const Vector& row1, const Vector& row2, const Vector& row3, const Vector& row4) {
        data_[0] = row1[0];
        data_[1] = row1[1];
        data_[2] = row1[2];
        data_[3] = 0;

        data_[4] = row2[0];
        data_[5] = row2[1];
        data_[6] = row2[2];
        data_[7] = 0;

        data_[8] = row3[0];
        data_[9] = row3[1];
        data_[10] = row3[2];
        data_[11] = 0;

        data_[12] = row4[0];
        data_[13] = row4[1];
        data_[14] = row4[2];
        data_[15] = 1;
    }
    double operator()(size_t row, size_t col) const {
        return data_[4 * row + col];
    }
    double& operator()(size_t row, size_t col) {
        return data_[4 * row + col];
    }
    friend Vector MatrixDotProduct(const Vector& vec, const Mat& matrix) {
        std::array<double, 4> res{0, 0, 0, 0}, vec4 = {vec[0], vec[1], vec[2], 1.};

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                res[i] += vec4[j] * matrix(j, i);
            }
        }

        return {res[0] / res[3], res[1] / res[3], res[2] / res[3]};
    }

private:
    std::array<double, 16> data_;
};
