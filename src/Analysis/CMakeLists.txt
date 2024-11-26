if(enable_GLOG)
  set(compile_options_tmp ${base_compile_options} " -D_GLOG")
elseif()
  set(compile_options_tmp ${base_compile_options})
endif()

string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

set(main_target "Analysis_AD")

add_executable(${main_target} ${main_src})
target_sources(${main_target}
  PRIVATE
    main.cpp
)
target_compile_options(${main_target} PUBLIC ${compile_options_tmp})

if(compiler STREQUAL intel)
  target_link_libraries(${main_target} PRIVATE
    lib
  )
elseif(compiler STREQUAL gcc)
  target_link_libraries(${main_target} PRIVATE
    lib
  )
endif()

target_include_directories(${main_target} PRIVATE
  ${PROJECT_SOURCE_DIR}/src/lib
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})