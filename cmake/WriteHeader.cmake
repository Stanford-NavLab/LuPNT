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
      # Find all .h files under this top-level subdirectory, recursing into nested subdirectories
      # (e.g. numerics/filters/) so their headers are grouped under the top-level category.
      file(
        GLOB_RECURSE HEADER_FILES
        RELATIVE ${dir}/${subdir}
        ${dir}/${subdir}/*.h
      )
      # `plasma` is a self-contained, vendored module included directly by its consumers via
      # explicit paths (and defines symbols that clash with the top-level core/ headers); keep it
      # out of the umbrella header.
      list(FILTER HEADER_FILES EXCLUDE REGEX "(^|/)plasma/")
      list(SORT HEADER_FILES)
      foreach(file ${HEADER_FILES})
        # Create relative include path
        file(APPEND ${dir}/lupnt.h "#include \"lupnt/${subdir}/${file}\"\n")
      endforeach()
    endif()
  endforeach()
endfunction()
