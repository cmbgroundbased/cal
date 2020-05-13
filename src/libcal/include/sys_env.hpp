/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#ifndef CAL_SYS_ENVIRONMENT_HPP
#define CAL_SYS_ENVIRONMENT_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>

#include <sys_utils.hpp>


namespace cal {
class Environment {
    // Singleton containing runtime settings.

    public:

        // Singleton access
        static Environment & get();

        std::string log_level() const;
        void set_log_level(char const * level);
        std::vector <std::string> signals() const;
        std::vector <std::string> info() const;
        void print() const;
        bool use_mpi() const;
        bool function_timers() const;
        int max_threads() const;
        int current_threads() const;
        void set_threads(int nthread);
        std::string version() const;
        int64_t tod_buffer_length() const;

    private:

        // This class is a singleton- constructor is private.
        Environment();

        std::string loglvl_;
        std::vector <std::string> signals_avail_;
        std::map <std::string, int> signals_value_;
        std::map <std::string, bool> signals_enabled_;
        bool have_mpi_;
        bool use_mpi_;
        bool func_timers_;
        bool at_nersc_;
        bool in_slurm_;
        int max_threads_;
        int cur_threads_;
        std::string git_version_;
        std::string release_version_;
        std::string version_;
        int64_t tod_buffer_length_;
};
}

#endif // ifndef CAL_ENVIRONMENT_HPP
