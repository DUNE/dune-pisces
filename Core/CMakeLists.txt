set(LIBRARY ${PROJECT_NAME}Core)

file(GLOB HEADER_FILES *.h)
file(GLOB SOURCES *.cxx)

add_library(${LIBRARY} SHARED
            ${HEADER_FILES}
	    ${SOURCES}
)

target_link_directories(${LIBRARY}
        PUBLIC $ENV{CAFANACORE_LIB}
        PUBLIC $ENV{CAFANA_LIB}
)
target_link_libraries(${LIBRARY}
        PUBLIC CAFAnaCoreExt
        PUBLIC CAFAnaAnalysis
        PUBLIC CAFAnaCore
        PUBLIC CAFAnaCuts
        PUBLIC CAFAnaDecomp
        PUBLIC CAFAnaExperiment
        PUBLIC CAFAnaExtrap
        PUBLIC CAFAnaFit
        PUBLIC CAFAnaPrediction
        PUBLIC CAFAnaSysts
        PUBLIC CAFAnaVars
)
link_root(${LIBRARY})

install(TARGETS ${LIBRARY} LIBRARY DESTINATION lib)
install(FILES ${HEADER_FILES} DESTINATION include/${PROJECT_NAME}/Core)
