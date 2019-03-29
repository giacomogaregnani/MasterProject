file(GLOB SUBSOURCES SubLibrary/*.cpp SubLibrary/*.hpp)

add_library(Class
        Class.hpp
        Class.cpp
        ${SUBSOURCES}
        )

