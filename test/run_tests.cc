/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh 1
 *  Dependencies: None
 *
 * Runs all Google Test fixtures in the mht/test folder.
 *************************************************************************/
#include "gtest/gtest.h"

/**
 * Runs all the Google Test fixtures in the mht/test folder.
 *
 * @author SCJ Robertson
 * @since 12/10/16 
 */
int main(int argc, char **argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
