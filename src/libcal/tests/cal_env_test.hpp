
// Copyright (c) 2015-2020 by the parties listed in the AUTHORS file.
// All rights reserved.  Use of this source code is governed by
// a BSD-style license that can be found in the LICENSE file.

TEST_F(CALenvTest, print) {
    auto & env = cal::Environment::get();
    env.print();
}


TEST_F(CALenvTest, setlog) {
    auto & env = cal::Environment::get();
    std::string check = env.log_level();
    ASSERT_STREQ(check.c_str(), "INFO");
    env.set_log_level("CRITICAL");
    check = env.log_level();
    ASSERT_STREQ(check.c_str(), "CRITICAL");
    env.set_log_level("INFO");
}
