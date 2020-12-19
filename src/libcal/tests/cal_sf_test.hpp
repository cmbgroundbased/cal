
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

#include <cal_test.hpp>


const int CALsfTest::size = 1000;


void compare_fast(double out, double expected) {
    double f32eps = std::numeric_limits <float>::epsilon();
    if ((fabs(out) < f32eps) && (fabs(expected) < f32eps)) {
        // Both are zero within expected error
        return;
    }
    EXPECT_FLOAT_EQ(out, expected);
    return;
}

void CALsfTest::SetUp() {
    angin.resize(size);
    sinout.resize(size);
    cosout.resize(size);
    xin.resize(size);
    yin.resize(size);
    atanout.resize(size);
    sqin.resize(size);
    sqout.resize(size);
    rsqin.resize(size);
    rsqout.resize(size);
    expin.resize(size);
    expout.resize(size);
    login.resize(size);
    logout.resize(size);

    double eps = 0.5 * std::numeric_limits <double>::epsilon();

    // Make sure to sample around zero and 2PI
    angin[0] = -3.0 * eps;
    angin[1] = -2.0 * eps;
    angin[2] = -eps;
    angin[3] = 0.0;
    angin[4] = eps;
    angin[5] = 2.0 * eps;
    angin[6] = 3.0 * eps;
    angin[size - 7] = 2.0 * cal::PI - 3.0 * eps;
    angin[size - 6] = 2.0 * cal::PI - 2.0 * eps;
    angin[size - 5] = 2.0 * cal::PI - eps;
    angin[size - 4] = 2.0 * cal::PI;
    angin[size - 3] = 2.0 * cal::PI + eps;
    angin[size - 2] = 2.0 * cal::PI + 2.0 * eps;
    angin[size - 1] = 2.0 * cal::PI + 3.0 * eps;
    for (int i = 7; i < size - 7; ++i) {
        angin[i] = (double)i * (2.0 * cal::PI / (double)(size - 1));
    }

    for (int i = 0; i < size; ++i) {
        sinout[i] = ::sin(angin[i]);
        cosout[i] = ::cos(angin[i]);
        xin[i] = cosout[i];
        yin[i] = sinout[i];
        atanout[i] = ::atan2(yin[i], xin[i]);
    }

    sqin[0] = 0.0;
    sqin[1] = eps;
    sqin[2] = 2.0 * eps;

    rsqin[0] = eps;
    rsqin[1] = 2.0 * eps;
    rsqin[2] = 3.0 * eps;

    for (int i = 3; i < size; ++i) {
        sqin[i] = (double)i / (double)size;
        rsqin[i] = sqin[i];
    }

    for (int i = 0; i < size; ++i) {
        sqout[i] = ::sqrt(sqin[i]);
        rsqout[i] = 1.0 / ::sqrt(rsqin[i]);
    }

    expin[0] = -1.0e300;
    expin[1] = -1.0e200;
    expin[2] = -1.0e100;

    login[0] = eps;
    login[1] = 2.0 * eps;
    login[2] = 3.0 * eps;

    for (int i = 3; i < size; ++i) {
        expin[i] = 0.005 * (double)i;
        login[i] = expout[i];
    }

    for (int i = 0; i < size; ++i) {
        expout[i] = ::exp(expin[i]);
        logout[i] = ::log(login[i]);
    }

    return;
}

void CALsfTest::TearDown() {
    return;
}

TEST_F(CALsfTest, trig) {
    cal::AlignedVector <double> comp1(size);
    cal::AlignedVector <double> comp2(size);

    cal::vsin(size, angin.data(), comp1.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_DOUBLE_EQ(sinout[i], comp1[i]);
    }

    cal::vcos(size, angin.data(), comp2.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_DOUBLE_EQ(cosout[i], comp2[i]);
    }

    cal::vsincos(size, angin.data(), comp1.data(), comp2.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_DOUBLE_EQ(sinout[i], comp1[i]);
        EXPECT_DOUBLE_EQ(cosout[i], comp2[i]);
    }

    cal::vatan2(size, yin.data(), xin.data(), comp1.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_DOUBLE_EQ(atanout[i], comp1[i]);
    }
}


TEST_F(CALsfTest, fasttrig) {
    cal::AlignedVector <double> comp1(size);
    cal::AlignedVector <double> comp2(size);

    cal::vfast_sin(size, angin.data(), comp1.data());
    for (int i = 0; i < size; ++i) {
        // std::cout << "sin ang = " << angin[i] << ", libm = " << sinout[i] <<
        // ", fast = " << comp1[i] << std::endl;
        compare_fast(comp1[i], sinout[i]);
    }

    cal::vfast_cos(size, angin.data(), comp2.data());
    for (int i = 0; i < size; ++i) {
        // std::cout << "cos ang = " << angin[i] << ", libm = " <<
        //     cosout[i] << ", fast = " << comp2[i] << std::endl;
        compare_fast(comp2[i], cosout[i]);
    }

    cal::vfast_sincos(size, angin.data(), comp1.data(), comp2.data());
    for (int i = 0; i < size; ++i) {
        // std::cout << "sincos ang = " << angin[i] << ", libm = " << sinout[i] << " "
        // <<
        //     cosout[i] << ", fast = " << comp1[i] << " " << comp2[i] << std::endl;
        compare_fast(comp1[i], sinout[i]);
        compare_fast(comp2[i], cosout[i]);
    }

    cal::vfast_atan2(size, yin.data(), xin.data(), comp1.data());
    for (int i = 0; i < size; ++i) {
        // std::cout << "atan in = " << yin[i] << " " << xin[i] << ", libm = " <<
        // atanout[i] << ", fast = " << comp1[i] << std::endl;
        compare_fast(comp1[i], atanout[i]);
    }
}


TEST_F(CALsfTest, sqrtlog) {
    cal::AlignedVector <double> comp(size);

    cal::vsqrt(size, sqin.data(), comp.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_DOUBLE_EQ(sqout[i], comp[i]);
    }

    cal::vrsqrt(size, rsqin.data(), comp.data());
    for (int i = 0; i < size; ++i) {
        if (std::isfinite(comp[i]) || std::isfinite(rsqout[i])) {
            EXPECT_DOUBLE_EQ(rsqout[i], comp[i]);
        }
    }

    cal::vexp(size, expin.data(), comp.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_DOUBLE_EQ(expout[i], comp[i]);
    }

    cal::vlog(size, login.data(), comp.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_DOUBLE_EQ(logout[i], comp[i]);
    }
}


TEST_F(CALsfTest, fast_sqrtlog) {
    cal::AlignedVector <double> comp(size);

    cal::vfast_sqrt(size, sqin.data(), comp.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_FLOAT_EQ(sqout[i], comp[i]);
    }

    cal::vfast_rsqrt(size, rsqin.data(), comp.data());
    for (int i = 0; i < size; ++i) {
        if (std::isfinite(comp[i]) || std::isfinite(rsqout[i])) {
            EXPECT_FLOAT_EQ(rsqout[i], comp[i]);
        }
    }

    cal::vfast_exp(size, expin.data(), comp.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_FLOAT_EQ(expout[i], comp[i]);
    }

    cal::vfast_log(size, login.data(), comp.data());
    for (int i = 0; i < size; ++i) {
        EXPECT_FLOAT_EQ(logout[i], comp[i]);
    }
}


TEST_F(CALsfTest, fast_erfinv) {
    cal::AlignedVector <double> in = {
        -9.990000e-01,
        -7.770000e-01,
        -5.550000e-01,
        -3.330000e-01,
        -1.110000e-01,
        1.110000e-01,
        3.330000e-01,
        5.550000e-01,
        7.770000e-01,
        9.990000e-01
    };

    cal::AlignedVector <double> check = {
        -2.326753765513524e+00,
        -8.616729665092674e-01,
        -5.400720684419686e-01,
        -3.042461029341061e-01,
        -9.869066534119145e-02,
        9.869066534119164e-02,
        3.042461029341063e-01,
        5.400720684419690e-01,
        8.616729665092677e-01,
        2.326753765513546e+00
    };

    cal::AlignedVector <double> out(10);

    cal::vfast_erfinv(10, in.data(), out.data());

    for (int i = 0; i < 10; ++i) {
        EXPECT_FLOAT_EQ(out[i], check[i]);
    }
}
