
if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_IMAGE_ANALYZE_IO__)
endif ()


set(image_io_header_files
        ImageIOExport.h
        ImageIOBase.h
        ImageIOAnalyze.h)

set(image_io_src_files
        ImageIOBase.cpp
        ImageIOAnalyze.cpp)

add_library(gadgetron_toolbox_image_analyze_io SHARED ${image_io_header_files} ${image_io_src_files})
set_target_properties(gadgetron_toolbox_image_analyze_io PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
target_link_libraries(gadgetron_toolbox_image_analyze_io gadgetron_toolbox_log gadgetron_toolbox_cpucore boost)
target_include_directories(gadgetron_toolbox_image_analyze_io
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})


install(TARGETS gadgetron_toolbox_image_analyze_io
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        ${image_io_header_files}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
