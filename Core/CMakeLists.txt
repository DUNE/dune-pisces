set(LIBRARY ${PROJECT_NAME}Core)

file(GLOB HEADER_FILES *.h)
file(GLOB SOURCES *.cxx)

add_library(${LIBRARY} SHARED
            ${HEADER_FILES}
	    ${SOURCES}
)

# currently depends on full CAFAna... at some point hope this won't be the case
target_link_directories(${LIBRARY}
        PUBLIC $ENV{CAFANACORE_LIB}
        PUBLIC $ENV{CAFANA_LIB}
        PUBLIC $ENV{DUNEANAOBJ_LIB}
)
# may not need all these libraries -- will tidy it up as the dependency on CAFAna slowly unravels over time (well, that's the idea...)
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
        PUBLIC duneanaobj_StandardRecordProxy
)
link_root(${LIBRARY})

install(TARGETS ${LIBRARY} LIBRARY DESTINATION lib)
install(FILES ${HEADER_FILES} DESTINATION include/${PROJECT_NAME}/Core)
