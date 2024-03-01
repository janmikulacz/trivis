/**
 * File:   logging.cc
 *
 * Date:   18.11.2021
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_plus/utils/log.h"

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>

using namespace trivis_plus::utils;

namespace bl = boost::log;

void trivis_plus::utils::InitLogging(const severity_level &level) {
    bl::add_common_attributes();
    bl::register_simple_formatter_factory<severity_level, char>("Severity");
    std::string log_format = std::string("") + "[%TimeStamp%]" + (std::string(kLibName).empty() ? "" : " [" + std::string(kLibName) + "]") + " [%Severity%] >> %Message%";
    bl::add_console_log(std::cout, bl::keywords::format = log_format);
    bl::core::get()->set_filter(bl::trivial::severity >= bl::trivial::trace);
    bl::core::get()->set_filter(bl::trivial::severity >= level);
}