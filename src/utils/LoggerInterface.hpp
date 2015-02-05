#ifndef LOGGERINTERFACE_HPP
#define LOGGERINTERFACE_HPP

#ifdef HAS_CXX11

#include "Logger.hpp"
#include "Timer.hpp"

#define LOG_INIT logging::logger<logging::FileLogPolicy>::Instance().initialize
#define LOG      logging::logger<logging::FileLogPolicy>::Instance().print<logging::coarse>
#define LOG_FINE logging::logger<logging::FileLogPolicy>::Instance().print<logging::fine>
#define LOG_ALL  logging::logger<logging::FileLogPolicy>::Instance().print<logging::everything>
#define LOG_TIME logging::logger<logging::FileLogPolicy>::Instance().print<logging::timings>(Timer::TheTimer())

#else // HAS_CXX11

#define LOG_INIT(...)
#define LOG(...)
#define LOG_FINE(...)
#define LOG_ALL(...)
#define LOG_TIME

#endif // HAS_CXX11
 
#endif // LOGGERINTERFACE_HPP
