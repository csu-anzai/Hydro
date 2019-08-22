INCLUDE_DIRECTORIES(${LAVAflow_INCLUDES})

set(LAVAUTIL_SOURCES
	ScalarList.cxx
	IndexOperations.cxx
	SortOperations.cxx
    ArrayOperations.cxx
    FlashIndexing.cxx
)

ADD_LIBRARY(LAVAUTIL STATIC ${LAVAUTIL_SOURCES})

list(APPEND LAVAflow_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
set(LAVAflow_INCLUDES ${LAVAflow_INCLUDES} PARENT_SCOPE)
