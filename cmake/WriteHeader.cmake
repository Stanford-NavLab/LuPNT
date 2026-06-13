# Create a single header file lupnt.h with all includes
function(write_header dir)
  file(WRITE ${dir}/lupnt.h "#pragma once\n")
  file(
    GLOB SUBDIRS
    LIST_DIRECTORIES true
    RELATIVE ${dir}
    ${dir}/*
  )
  list(FILTER SUBDIRS EXCLUDE REGEX "^\\.|\\.\\./")
  foreach(subdir ${SUBDIRS})
    if(IS_DIRECTORY ${dir}/${subdir})
      # Format the directory name for a comment
      string(REGEX REPLACE "^${dir}/" "" subdir_name ${subdir})
      file(APPEND ${dir}/lupnt.h "\n// ${subdir_name}\n")
      # Find all .h files in the current subdirectory
      file(
        GLOB HEADER_FILES
        RELATIVE ${dir}/${subdir}
        ${dir}/${subdir}/*.h
      )
      list(SORT HEADER_FILES)
      foreach(file ${HEADER_FILES})
        # Create relative include path
        file(APPEND ${dir}/lupnt.h "#include \"lupnt/${subdir}/${file}\"\n")
      endforeach()
    endif()
  endforeach()
endfunction()
