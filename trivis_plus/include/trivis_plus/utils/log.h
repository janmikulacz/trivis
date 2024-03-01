/**
 * File:   logging.h
 *
 * Date:   18.11.2021
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_UTILS_LOG_H_
#define TRIVIS_PLUS_UTILS_LOG_H_

#include <boost/log/trivial.hpp>

#ifndef LOG_TRC
#define LOG_TRC BOOST_LOG_TRIVIAL(trace)
#endif

#ifndef LOG_DBG
#define LOG_DBG BOOST_LOG_TRIVIAL(debug)
#endif

#ifndef LOG_INF
#define LOG_INF BOOST_LOG_TRIVIAL(info)
#endif

#ifndef LOG_WRN
#define LOG_WRN BOOST_LOG_TRIVIAL(warning)
#endif

#ifndef LOG_ERR
#define LOG_ERR BOOST_LOG_TRIVIAL(error)
#endif

#ifndef LOG_FTL
#define LOG_FTL BOOST_LOG_TRIVIAL(fatal)
#endif

#ifndef LOGF_TRC
#define LOGF_TRC(msg) LOG_TRC << trivis_plus::utils::Color(trivis_plus::utils::kBlue) << msg << trivis_plus::utils::kReset
#endif

#ifndef LOGF_DBG
#define LOGF_DBG(msg) LOG_DBG << trivis_plus::utils::Color(trivis_plus::utils::kGreen) << msg << trivis_plus::utils::kReset
#endif

#ifndef LOGF_INF
#define LOGF_INF(msg) LOG_INF << msg << trivis_plus::utils::kReset
#endif

#ifndef LOGF_WRN
#define LOGF_WRN(msg) LOG_WRN << trivis_plus::utils::Color(trivis_plus::utils::kYellow) << msg << trivis_plus::utils::kReset
#endif

#ifndef LOGF_ERR
#define LOGF_ERR(msg) LOG_ERR << trivis_plus::utils::Color(trivis_plus::utils::kRed) << msg << trivis_plus::utils::kReset
#endif

#ifndef LOGF_FTL
#define LOGF_FTL(msg) LOG_FTL << trivis_plus::utils::Style(trivis_plus::utils::kBold, trivis_plus::utils::kRed) << msg << trivis_plus::utils::kReset
#endif

#ifndef LOG_NAME
#define LOG_NAME ""
#endif

namespace trivis_plus::utils {

using severity_level = boost::log::trivial::severity_level;

static constexpr const char *kLibName = LOG_NAME;
static constexpr const char *kReset = "\u001b[0m";
static constexpr const char *kEndL = "\u001b[0m\n";

enum Style : int {
    kBold = 1,
    kUnderline = 4,
    kReversed = 7
};

enum Color : int {
    kBlack = 0,
    kRed = 1,
    kGreen = 2,
    kYellow = 3,
    kBlue = 4,
    kMagenta = 5,
    kCyan = 6,
    kWhite = 7,
};

inline std::string Color(enum Color txt_color) {
    return "\u001b[3" + std::to_string(txt_color) + "m";
}

inline std::string BgColor(enum Color bg_color) {
    return "\u001b[4" + std::to_string(bg_color) + "m";
}

inline std::string BgColor(enum Color bg_color, enum Color txt_color) {
    return "\u001b[4" + std::to_string(bg_color) + ";3" + std::to_string(txt_color) + "m";
}

inline std::string BgColor(enum Color bg_color, enum Style style) {
    return "\u001b[4" + std::to_string(bg_color) + ";" + std::to_string(style) + "m";
}

inline std::string Style(enum Style style) {
    return "\u001b[" + std::to_string(style) + "m";
}

inline std::string Style(enum Style style, enum Color txt_color) {
    return "\u001b[" + std::to_string(style) + ";3" + std::to_string(txt_color) + "m";
}

inline std::string Style(enum Style style, enum Color txt_color, enum Color bg_color) {
    return "\u001b[" + std::to_string(style) + ";3" + std::to_string(txt_color) + ";4" + std::to_string(bg_color) + "m";
}

void InitLogging(const boost::log::trivial::severity_level &level);

}

#endif //TRIVIS_PLUS_UTILS_LOG_H_
