file(GLOB SOURCES CONFIGURE_DEPENDS *.c)

add_library(timing STATIC ${SOURCES})
