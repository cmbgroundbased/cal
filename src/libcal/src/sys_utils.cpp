/*
   Copyright (c) 2015-2018 by the parties listed in the AUTHORS file.
   All rights reserved.  Use of this source code is governed by
   a BSD-style license that can be found in the LICENSE file.
 */

#include <sys_env.hpp>
#include <sys_utils.hpp>

#include <cstring>
#include <sstream>

#include <vector>
#include <algorithm>


std::string cal::format_here(std::pair <std::string, int> const & here) {
    std::ostringstream h;
    h << "file \"" << here.first << "\", line " << here.second;
    return std::string(h.str());
}

void * cal::aligned_alloc(size_t size, size_t align) {
    void * mem = NULL;
    int ret = posix_memalign(&mem, align, size);
    if (ret != 0) {
        auto here = cal_HERE();
        auto log = cal::Logger::get();
        std::ostringstream o;
        o << "cannot allocate " << size
          << " bytes of memory with alignment " << align;
        log.error(o.str().c_str(), here);
        throw std::runtime_error(o.str().c_str());
    }
    memset(mem, 0, size);
    return mem;
}

void cal::aligned_free(void * ptr) {
    free(ptr);
    return;
}

cal::Timer::Timer() {
    clear();
}

cal::Timer::Timer(double init_time, size_t init_calls) {
    start_ = time_point();
    stop_ = time_point();
    running_ = false;
    total_ = init_time;
    calls_ = init_calls;
}

void cal::Timer::start() {
    if (!running_) {
        start_ = std::chrono::high_resolution_clock::now();
        running_ = true;
        calls_++;
    }
    return;
}

void cal::Timer::stop() {
    if (running_) {
        stop_ = std::chrono::high_resolution_clock::now();
        std::chrono::duration <double> elapsed =
            std::chrono::duration_cast <std::chrono::duration <double> >
                (stop_ - start_);
        total_ += elapsed.count();
        running_ = false;
    }
    return;
}

void cal::Timer::clear() {
    start_ = time_point();
    stop_ = time_point();
    running_ = false;
    calls_ = 0;
    total_ = 0.0;
    return;
}

double cal::Timer::seconds() const {
    if (running_) {
        auto here = cal_HERE();
        auto log = cal::Logger::get();
        std::string msg("Timer is still running!");
        log.error(msg.c_str(), here);
        throw std::runtime_error(msg.c_str());
    }
    return total_;
}

double cal::Timer::elapsed_seconds() const {
    /* Return the current reading on the timer without incrementing
       the calls_ counter or stopping the timer
     */
    if (not running_) {
        auto here = cal_HERE();
        auto log = cal::Logger::get();
        std::string msg("Timer is not running!");
        log.error(msg.c_str(), here);
        throw std::runtime_error(msg.c_str());
    }
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration <double> elapsed =
        std::chrono::duration_cast <std::chrono::duration <double> >
            (now - start_);
    return total_ + elapsed.count();
}

size_t cal::Timer::calls() const {
    if (running_) {
        auto here = cal_HERE();
        auto log = cal::Logger::get();
        std::string msg("Timer is still running!");
        log.error(msg.c_str(), here);
        throw std::runtime_error(msg.c_str());
    }
    return calls_;
}

bool cal::Timer::is_running() const {
    return running_;
}

void cal::Timer::report(char const * message) {
    double t = seconds();
    cal::Logger & logger = cal::Logger::get();
    std::ostringstream msg;

    msg.precision(2);
    msg << std::fixed << message << ":  " << t << " seconds ("
        << calls_ << " calls)";
    logger.info(msg.str().c_str());
    return;
}

void cal::Timer::report_clear(char const * message) {
    bool was_running = running_;
    if (was_running) {
        stop();
    }
    report(message);
    clear();
    if (was_running) {
        start();
    }
    return;
}

void cal::Timer::report_elapsed(char const * message) {
    /* Report elapsed time from a running timer without changing its state
     */
    double t = elapsed_seconds();
    cal::Logger & logger = cal::Logger::get();
    std::ostringstream msg;

    msg.precision(2);
    msg << std::fixed << message << ":  " << t << " seconds ("
        << calls_ << " calls)";
    logger.info(msg.str().c_str());
    return;
}

cal::GlobalTimers::GlobalTimers() {
    data.clear();
}

cal::GlobalTimers & cal::GlobalTimers::get() {
    static cal::GlobalTimers instance;

    return instance;
}

std::vector <std::string> cal::GlobalTimers::names() const {
    std::vector <std::string> ret;
    for (auto const & it : data) {
        ret.push_back(it.first);
    }
    std::stable_sort(ret.begin(), ret.end());
    return ret;
}

void cal::GlobalTimers::start(std::string const & name) {
    if (data.count(name) == 0) {
        data[name].clear();
    }
    data.at(name).start();
    return;
}

void cal::GlobalTimers::clear(std::string const & name) {
    data[name].clear();
    return;
}

void cal::GlobalTimers::stop(std::string const & name) {
    if (data.count(name) == 0) {
        auto here = cal_HERE();
        auto log = cal::Logger::get();
        std::ostringstream o;
        o << "Cannot stop timer " << name << " which does not exist";
        log.error(o.str().c_str(), here);
        throw std::runtime_error(o.str().c_str());
    }
    data.at(name).stop();
    return;
}

double cal::GlobalTimers::seconds(std::string const & name) const {
    if (data.count(name) == 0) {
        auto here = cal_HERE();
        auto log = cal::Logger::get();
        std::ostringstream o;
        o << "Cannot get seconds for timer " << name
          << " which does not exist";
        log.error(o.str().c_str(), here);
        throw std::runtime_error(o.str().c_str());
    }
    return data.at(name).seconds();
}

size_t cal::GlobalTimers::calls(std::string const & name) const {
    if (data.count(name) == 0) {
        auto here = cal_HERE();
        auto log = cal::Logger::get();
        std::ostringstream o;
        o << "Cannot get seconds for timer " << name
          << " which does not exist";
        log.error(o.str().c_str(), here);
        throw std::runtime_error(o.str().c_str());
    }
    return data.at(name).calls();
}

bool cal::GlobalTimers::is_running(std::string const & name) const {
    if (data.count(name) == 0) {
        return false;
    }
    return data.at(name).is_running();
}

void cal::GlobalTimers::stop_all() {
    for (auto & tm : data) {
        tm.second.stop();
    }
    return;
}

void cal::GlobalTimers::clear_all() {
    for (auto & tm : data) {
        tm.second.clear();
    }
    return;
}

void cal::GlobalTimers::report() {
    stop_all();
    std::vector <std::string> names;
    for (auto const & tm : data) {
        names.push_back(tm.first);
    }
    std::stable_sort(names.begin(), names.end());
    std::ostringstream msg;
    for (auto const & nm : names) {
        msg.str("");
        msg << "Global timer: " << nm;
        data.at(nm).report(msg.str().c_str());
    }
    return;
}

cal::Logger::Logger() {
    // Prefix for messages
    prefix_ = std::string("cal ");
    return;
}

cal::Logger & cal::Logger::get() {
    static cal::Logger instance;

    // Check the level every time we get a reference to the singleton,
    // in case the level has been changed manually during runtime.
    instance.check_level();
    return instance;
}

void cal::Logger::check_level() {
    auto & env = cal::Environment::get();
    std::string val = env.log_level();
    if (strncmp(val.c_str(), "DEBUG", 5) == 0) {
        level_ = log_level::debug;
    } else if (strncmp(val.c_str(), "INFO", 4) == 0) {
        level_ = log_level::info;
    } else if (strncmp(val.c_str(), "WARNING", 7) == 0) {
        level_ = log_level::warning;
    } else if (strncmp(val.c_str(), "ERROR", 5) == 0) {
        level_ = log_level::error;
    } else if (strncmp(val.c_str(), "CRITICAL", 8) == 0) {
        level_ = log_level::critical;
    } else {
        level_ = log_level::none;
    }
    return;
}

void cal::Logger::debug(char const * msg) {
    if (level_ <= log_level::debug) {
        fprintf(stdout, "%sDEBUG: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void cal::Logger::debug(char const * msg,
                          std::pair <std::string, int> const & here) {
    if (level_ <= log_level::debug) {
        std::string hstr = cal::format_here(here);
        fprintf(stdout, "%sDEBUG: %s (%s)\n", prefix_.c_str(), msg,
                hstr.c_str());
        fflush(stdout);
    }
    return;
}

void cal::Logger::info(char const * msg) {
    if (level_ <= log_level::info) {
        fprintf(stdout, "%sINFO: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void cal::Logger::info(char const * msg,
                         std::pair <std::string, int> const & here) {
    if (level_ <= log_level::info) {
        std::string hstr = cal::format_here(here);
        fprintf(stdout, "%sINFO: %s (%s)\n", prefix_.c_str(), msg,
                hstr.c_str());
        fflush(stdout);
    }
    return;
}

void cal::Logger::warning(char const * msg) {
    if (level_ <= log_level::warning) {
        fprintf(stdout, "%sWARNING: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void cal::Logger::warning(char const * msg,
                            std::pair <std::string, int> const & here) {
    if (level_ <= log_level::warning) {
        std::string hstr = cal::format_here(here);
        fprintf(stdout, "%sWARNING: %s (%s)\n", prefix_.c_str(), msg,
                hstr.c_str());
        fflush(stdout);
    }
    return;
}

void cal::Logger::error(char const * msg) {
    if (level_ <= log_level::error) {
        fprintf(stdout, "%sERROR: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void cal::Logger::error(char const * msg,
                          std::pair <std::string, int> const & here) {
    if (level_ <= log_level::error) {
        std::string hstr = cal::format_here(here);
        fprintf(stdout, "%sERROR: %s (%s)\n", prefix_.c_str(), msg,
                hstr.c_str());
        fflush(stdout);
    }
    return;
}

void cal::Logger::critical(char const * msg) {
    if (level_ <= log_level::critical) {
        fprintf(stdout, "%sCRITICAL: %s\n", prefix_.c_str(), msg);
        fflush(stdout);
    }
    return;
}

void cal::Logger::critical(char const * msg,
                             std::pair <std::string, int> const & here) {
    if (level_ <= log_level::critical) {
        std::string hstr = cal::format_here(here);
        fprintf(stdout, "%sCRITICAL: %s (%s)\n", prefix_.c_str(), msg,
                hstr.c_str());
        fflush(stdout);
    }
    return;
}
