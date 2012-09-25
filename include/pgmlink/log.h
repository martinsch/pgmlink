/*
 * Log.h
 *
 * Copyright (c) 2009 Bernhard Kausler <bernhard.kausler@iwr.uni-karlsruhe.de>
 * Initial version based on ideas of Petru Marginea <petru.marginean@gmail.com>
 * (see: http://www.ddj.com/cpp/201804215)
 *
 */

#ifndef LOG_H_
#define LOG_H_

#include <ctime>
#include <cstdio>

#include <string>
#include <sstream>
#include <iostream>



/**
 * @page log Logging
 *
 *
 *
 * @section usage Usage
 * Logging is provided to @c stderr via the @c LOG macro:
 * @code
 * LOG(logINFO) << "Hello" << username << "no endl, will be appended automatically";
 * @endcode
 *
 * You can set a plateu logging level:
 * @code
 * #undef FILELOG_MAX_LEVEL
 * #define FILELOG_MAX_LEVEL tracking::logWARNING
 * @endcode
 *
 * This can be used to set a global logging level or to construct a more fine grained
 * structure (for example, different logging levels per file or code segment).
 *
 * The logging macro LOG ensures, that logging code will only be compiled into the final
 * binary, if it corresponds to the set plateu logging level.
 *
 * <em>
 * As a consequence, don't use the Log API directly, but only via the provided macros.
 * </em>
 *
 *
 *
 * @section logginglevels Provided logging levels
 * These are the available logging levels ordered from the highest to the deepest:
 * @code
 * logERROR, logWARNING, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4
 * @endcode
 *
 *
 *
 * @section disablelogging Turning off logging
 * If logging and the accompanying code is not desired in the final binary, a special logging
 * level is provided:
 * @code
 * logNO_LOGGING
 * @endcode
 * Set it as the plateu logging level to disable logging.
 *
 * <em>
 * Don't use logNO_LOGGING as a paramter to the logging macro.
 * </em>
 *
 *
 *
 * @author Bernhard X. Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-04-02
 */



namespace Tracking
{



// nowTime()
/**
 * The current time as a string.
 *
 * The exact presentation of the string is locale dependent, but should be the following
 * in most cases:
 * @code
 * 15:50:57.979
 * @endcode
 * The precision is up to milliseconds.
 */
inline std::string nowTime();



// LogLevel
/**
 * The logging levels, which can be used.
 *
 * Assign one of these logging levels to each message.
 * Don't use the highest level 'logNO_LOGGING'. It is needed to turn off logging and would
 * default to logINFO if used.
 */
enum LogLevel {logNO_LOGGING, logERROR, logWARNING, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4};



// Output2FILE
/**
 * Redirector of the logging stream to a file handle.
 *
 * At the moment, the file handle is hardcoded to 'stderr'.
 *
 * Use it in conjunction with the Log<T> class: Log<Output2FILE>.
 *
 * @author Bernhard X. Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-04-02
 */
class Output2FILE
{
public:
    // getRedirect()
    /**
     * The file handle to which the logging stream is redirected.
     */
    static FILE*& getRedirect();

    // output()
    /**
     * Writes message to file handle.
     *
     * This function is used in the Log<T> and mandatory for every Redirector.
     */
    static void output(const std::string& msg);
};



// getRedirect()
inline FILE*& Tracking::Output2FILE::getRedirect()
{
    static FILE* pStream = stderr;
    return pStream;
}



// output()
inline void Tracking::Output2FILE::output(const std::string& msg)
{
    FILE* pStream = getRedirect();
    if (!pStream) {
        return;
    }
    fprintf(pStream, "%s", msg.c_str());
    fflush(pStream);
}



//Log<T>
/**
 * A thread-safe logging tool
 *
 * The Log<T> passes its internal logging stream to a Redirector T.
 *
 * You have to provide a Redirector with the following interface:
 * @code
 * void T::output(const std::string& msg)
 * @endcode
 *
 * @author Bernhard X. Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-04-02
 */
template <typename T>
class Log
{
public:
    Log();
    virtual ~Log();

    // get()
    /**
     * The internal logging stream.
     *
     * Calling this function writes something like the following into the internal logging
     * stream (depending on the chosen logging level and the current time):
     * @code
     * '- 16:17:23.714 WARNING: '
     * @endcode
     * Afterwards, it returns a reference to that internal stream.
     */
    std::ostringstream& get(LogLevel level = logINFO);

    // getReportingLevel
    /**
     * Returns the deepest logging level available.
     *
     * This function can be used as a safeguard in the logging macros to defend against illegal
     * user defined global logging levels.
     *
     * @return The deepest Tracking::LogLevel available.
     */
    static LogLevel& getReportingLevel();

    // toString()
    /**
     * Converts enumerated logging level to a string.
     *
     * For example: logINFO to "INFO".
     */
    static std::string toString(LogLevel level);

    // fromString()
    /**
     * Converts from string to a enumerated logging level
     *
     * For example: "INFO" to logINFO.
     */
    static LogLevel fromString(const std::string& level);

protected:
    // os_
    /**
     * Internal string stream used for logging messages.
     *
     * In the destructor, this stream is written to the Redirector T.
     */
    std::ostringstream os_;

private:
    // Declare copy constructor etc. as private, since we don't want them to be used.
    Log(const Log&);
    Log& operator =(const Log&);
};



// Log()
template <typename T>
Tracking::Log<T>::Log()
{
}



// get()
template <typename T>
std::ostringstream& Tracking::Log<T>::get(LogLevel level)
{
    // check for valid logging level
    LogLevel ll = level;
    if ( ll <= logNO_LOGGING || ll > getReportingLevel() ) {
        Tracking::Log<T>().get(Tracking::logWARNING) << "Log<T>::get(): Invalid logging level '" << ll << "'. Using INFO level as default.";
        ll = logINFO;
    }

    // print standard logging preambel to logging stream
    os_ << "- " << Tracking::nowTime();
    os_ << " " << toString(ll) << ": ";
    os_ << std::string(ll > Tracking::logDEBUG ? ll - Tracking::logDEBUG : 0, '\t');

    return os_;
}



// ~Log()
template <typename T>
Tracking::Log<T>::~Log()
{
    os_ << std::endl;
    T::output(os_.str());
}



// getReportingLevel()
template <typename T>
Tracking::LogLevel& Tracking::Log<T>::getReportingLevel()
{
    static Tracking::LogLevel reportingLevel = Tracking::logDEBUG4;
    return reportingLevel;
}



// toString()
template <typename T>
std::string Tracking::Log<T>::toString(LogLevel level)
{
    if (level > getReportingLevel() || level < logNO_LOGGING) {
        Tracking::Log<T>().get(Tracking::logWARNING) << "Log<T>::toString(): Unknown logging level '" << level << "'. Using INFO level as default.";
        return "INFO";
    }

    static const char* const buffer[] = {"NO_LOGGING", "ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"};
    return buffer[level];
}



// fromString()
template <typename T>
Tracking::LogLevel Tracking::Log<T>::fromString(const std::string& level)
{
    if (level == "DEBUG4")
        return Tracking::logDEBUG4;
    if (level == "DEBUG3")
        return Tracking::logDEBUG3;
    if (level == "DEBUG2")
        return Tracking::logDEBUG2;
    if (level == "DEBUG1")
        return Tracking::logDEBUG1;
    if (level == "DEBUG")
        return Tracking::logDEBUG;
    if (level == "INFO")
        return Tracking::logINFO;
    if (level == "WARNING")
        return Tracking::logWARNING;
    if (level == "ERROR")
        return Tracking::logERROR;
    if (level == "NO_LOGGING")
        return Tracking::logNO_LOGGING;

    // else
    Tracking::Log<T>().get(Tracking::logWARNING) << "Log<T>::fromString(): Unknown logging level '" << level << "'. Using INFO level as default.";
    return Tracking::logINFO;
}



// FILELOG_DECLSPEC
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#   if defined (BUILDING_FILELOG_DLL)
#       define FILELOG_DECLSPEC   __declspec (dllexport)
#   elif defined (USING_FILELOG_DLL)
#       define FILELOG_DECLSPEC   __declspec (dllimport)
#   else
#       define FILELOG_DECLSPEC
#   endif // BUILDING_DBSIMPLE_DLL
#else
#   define FILELOG_DECLSPEC
#endif // _WIN32



// FILELog
/**
 * An instance of LOG<T> which is writing to a FILE.
 */
class FILELOG_DECLSPEC FILELog : public Log<Output2FILE> {};
//typedef Log<Output2FILE> FILELog;



// FILELOG_MAX_LEVEL
#ifndef FILELOG_MAX_LEVEL
/**
 * The deepest logging level to be compiled into the code.
 *
 * Every logging message depper than that level will not be compiled into the code.
 */
#define FILELOG_MAX_LEVEL Tracking::logDEBUG4
#endif



// LOG()
/**
 * Logs to a file handle.
 *
 * This macro checks, if the logging level should be compiled. After that, it creates an
 * anonymous instance of FILELog and writes to its logging stream. Afterwards, the
 * anonymous object is destroyed and the logging stream flushed out to the FILE.
 *
 * Use it like this:
 * @code
 * LOG(logINFO) << "some logging" << 1224 << "no endl, will be appended automatically";
 * @endcode
 */
#define LOG(level) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else if (level > Tracking::FILELog::getReportingLevel() || !Tracking::Output2FILE::getRedirect()) ; \
    else Tracking::FILELog().get(level)



// nowTime()
// We have to do the following yaketiyak, because the standard <ctime> is not thread safe.
// (It is using static internal buffers in some functions like ctime() .)
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
    } // Temporarily close Tracking namespace to include the external Windows headers.

    // winsocks2.h has always to be included BEFORE windows.h
	// We don't use winsocks2 here, but it may be used in a file including this header.
	#include <winsock2.h>
	#include <windows.h>

// Reopen the Tracking namespace.
namespace Tracking {
inline std::string nowTime()
{
    const int MAX_LEN = 200;
    char buffer[MAX_LEN];

    // get time
    if (GetTimeFormatA(
                LOCALE_USER_DEFAULT,    // locale
                0,                      // time format flags
                0,                      // optional ptr to systemtime structure
                "HH':'mm':'ss",         // format
                buffer,                 // ptr to output buffer
                MAX_LEN)                // size of output buffer
            == 0) {
        return "Error in nowTime()";
    }

    // format time according to our format: "hh:mm:ss.ms"
    static DWORD first = GetTickCount();
    char result[100] = {0};
    std::sprintf(result, "%s.%03ld", buffer, (long)(GetTickCount() - first) % 1000);

    return result;
}

#else
} // Temporarily close Tracking namespace to inclue header files.
#include <sys/time.h>

// Reopen namespace Tracking.
namespace Tracking {
inline std::string nowTime()
{
    // get time
    time_t t;
    t = time(NULL);
    if ( t == static_cast<std::time_t>(-1) ) {
        return "Error_in_nowTime().time";
    }

    // convert time to local time
    tm r;
    if (localtime_r(&t, &r) == NULL) {
        return "Error_in_nowTime().localtime_r";
    }

    // convert localtime to a string
    char buffer[101];
    if (strftime(buffer, sizeof(buffer), "%X", &r) == 0) {
        return "Error_in_nowTime().strftime";
    }

    // format the string according to our format: "hh:mm:ss.ms"
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[101] = {0};
    std::sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);

    return result;
}

#endif //WIN32



} /* namespace Tracking */
#endif /* LOG_H__ */
