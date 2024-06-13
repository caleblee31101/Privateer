
find_path(ZSTD_INCLUDE_DIRS
  NAMES zstd.h
  HINTS ../build_debug/zstd/lib)
  #HINTS ${zstd_ROOT_DIR}/lib)

find_library(ZSTD_LIBRARIES
  NAMES zstd
  HINTS ../build_debug/zstd/lib)
  #HINTS ${zstd_ROOT_DIR}/lib)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(zstd
  DEFAULT_MSG
  ZSTD_LIBRARIES
  ZSTD_INCLUDE_DIRS
)

mark_as_advanced(
  ZSTD_LIBRARIES
  ZSTD_INCLUDE_DIRS)

if(ZSTD_FOUND AND NOT (TARGET zstd::zstd))
  add_library (zstd::zstd UNKNOWN IMPORTED)
  set_target_properties(zstd::zstd
    PROPERTIES
      IMPORTED_LOCATION ${ZSTD_LIBRARIES}
      INTERFACE_INCLUDE_DIRECTORIES ${ZSTD_INCLUDE_DIRS})
endif()

if (NOT ZSTD_FOUND)
    message(STATUS "Zstd not found!!!")
endif()