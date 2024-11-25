set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(lib STATIC)
target_sources(lib
  PRIVATE
    initialize.cpp
    main_Analysis.cpp
    main_Explicit.cpp
    main_Implicit.cpp
    main_DG.cpp
)



target_compile_options(lib PUBLIC ${compile_options_tmp})