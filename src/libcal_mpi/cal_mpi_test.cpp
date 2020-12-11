#include <cal_mpi.hpp>
#include <gtest/gtest.h>
#include <tests/cal_mpi_test.hpp>

#include <tests/cal_mpi_test_shmem.hpp>

int main(int argc, char * argv[]) {

    testing::InitGoogleTest(&argc, argv);

    cal::mpi_init(argc, argv);

    // Disable result printing from all processes except the root one.

    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    ::testing::TestEventListeners & listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    if (rank != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }

    // FIXME:  rank 0 of MPI_COMM_WORLD should create the test
    // output directory here if it does not exist.  Currently none of the
    // C++ unit tests write or read data, so this is not yet an issue.

    int result = RUN_ALL_TESTS();

    return result;
      
}
