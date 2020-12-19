#ifndef CAL_TEST_TEST_HPP
#define CAL_TEST_TEST_HPP

#include <cal.hpp>
#include <gtest/gtest.h>

class CALenvTest : public testing::Test {
    public:

        CALenvTest() {}

        ~CALenvTest() {}

        virtual void SetUp() {}

        virtual void TearDown() {}
};



class CALutilsTest : public ::testing::Test {
    public:

        CALutilsTest() {}

        ~CALutilsTest() {}

        virtual void SetUp() {}

        virtual void TearDown() {}
};


class CALqarrayTest : public ::testing::Test {
    public:

        CALqarrayTest() {}

        ~CALqarrayTest() {}

        virtual void SetUp();
        virtual void TearDown() {}

        cal::AlignedVector <double> q1;
        cal::AlignedVector <double> q1inv;
        cal::AlignedVector <double> q2;
        cal::AlignedVector <double> qtonormalize;
        cal::AlignedVector <double> qnormalized;
        cal::AlignedVector <double> vec;
        cal::AlignedVector <double> vec2;
        cal::AlignedVector <double> qeasy;
        cal::AlignedVector <double> mult_result;
        cal::AlignedVector <double> rot_by_q1;
        cal::AlignedVector <double> rot_by_q2;
};


class CALrngTest : public ::testing::Test {
    public:

        CALrngTest() {}

        ~CALrngTest() {}

        virtual void SetUp();
        virtual void TearDown() {}

        static const int64_t size;
        static const uint64_t counter[];
        static const uint64_t key[];
        static const uint64_t counter00[];
        static const uint64_t key00[];

        static const double array_gaussian[];
        static const double array_m11[];
        static const double array_01[];
        static const uint64_t array_uint64[];

        static const double array00_gaussian[];
        static const double array00_m11[];
        static const double array00_01[];
        static const uint64_t array00_uint64[];
};


class CALsfTest : public ::testing::Test {
    public:

        CALsfTest() {}

        ~CALsfTest() {}

        virtual void SetUp();
        virtual void TearDown();

        static const int size;
        cal::AlignedVector <double> angin;
        cal::AlignedVector <double> sinout;
        cal::AlignedVector <double> cosout;
        cal::AlignedVector <double> xin;
        cal::AlignedVector <double> yin;
        cal::AlignedVector <double> atanout;
        cal::AlignedVector <double> sqin;
        cal::AlignedVector <double> sqout;
        cal::AlignedVector <double> rsqin;
        cal::AlignedVector <double> rsqout;
        cal::AlignedVector <double> expin;
        cal::AlignedVector <double> expout;
        cal::AlignedVector <double> login;
        cal::AlignedVector <double> logout;
};


class CALhealpixTest : public ::testing::Test {
    public:

        CALhealpixTest() {}

        ~CALhealpixTest() {}

        virtual void SetUp() {}

        virtual void TearDown() {}
};



#endif // ifndef CAL_TEST_HPP
