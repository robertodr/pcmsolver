# Copyright (C) 2015-2016 Jonathan Müller <jonathanmueller.dev@gmail.com>
# This file is subject to the license terms in the LICENSE file
# found in the top-level directory of this distribution.

if(NOT COMP_API_VERSION)
    message(FATAL_ERROR "needs newer comp_base.cmake version")
endif()
comp_api_version(1)

comp_feature(decltype "int main() {int i; decltype(i) j;}" COMP_CPP11_FLAG)
comp_workaround(decltype
"
#ifndef ${COMP_PREFIX}DECLTYPE
    #if ${COMP_PREFIX}HAS_DECLTYPE
        #define ${COMP_PREFIX}DECLTYPE(x) decltype(x)
    #elif defined(__GNUC__)
        #define ${COMP_PREFIX}DECLTYPE(X) __typeof__(x)
    #else
        #error \"no decltype replacement available\"
    #endif
#endif" COMP_CPP98_FLAG)

if(COMP_API_VERSION VERSION_GREATER 1.0)
    comp_sd6_macro(decltype __cpp_decltype 200707)
endif()
