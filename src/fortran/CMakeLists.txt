add_executable(${ExampleName} ${ExampleName}.F90)

target_link_libraries(${ExampleName} PRIVATE opencmisslibs)

install(TARGETS ${ExampleName} DESTINATION bin)

configure_file(visualise.cmgui visualise.cmgui)
