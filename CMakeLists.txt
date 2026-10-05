check_include_file_cxx(stdint.h HAVE_STDINT_H)
if(HAVE_STDINT_H)
    add_definitions(-DHAVE_STDINT_H)
endif()

set(examples_as_tests_sources)
if(${ENABLE_EXAMPLES})
    set(examples_as_tests_sources
        #test/trace-matrix-propagation-model-examples-test-suite.cc
        )
endif()

build_lib(
    LIBNAME trace-matrix-propagation-model
    SOURCE_FILES model/trace-matrix-propagation-model.cc
                 helper/trace-matrix-propagation-model-helper.cc
    HEADER_FILES model/trace-matrix-propagation-model.h
                 helper/trace-matrix-propagation-model-helper.h
    LIBRARIES_TO_LINK ${libcore}
    TEST_SOURCES test/trace-matrix-propagation-model-test-suite.cc
                 ${examples_as_tests_sources}
)
