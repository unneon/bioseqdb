cmake_minimum_required(VERSION 3.14)

project(bioseqdb LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(PostgreSQL REQUIRED COMPONENTS Server)
find_program(PG_CONFIG pg_config)
execute_process(COMMAND ${PG_CONFIG} --pkglibdir OUTPUT_VARIABLE PG_CONFIG_PKGLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PG_CONFIG} --sharedir OUTPUT_VARIABLE PG_CONFIG_SHAREDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

find_package(BZip2 REQUIRED)

find_library(BWA_LIBRARIES bwa REQUIRED)
find_library(HTS_LIBRARIES hts REQUIRED)

add_library(bioseqdb-pg SHARED
        bioseqdb-pg/bwa.cpp
        bioseqdb-pg/extension.cpp
        bioseqdb-pg/sequence.cpp
        )
add_executable(bioseqdb-import
        bioseqdb-import/main.cpp
        )

target_include_directories(bioseqdb-pg PRIVATE ${PostgreSQL_TYPE_INCLUDE_DIR})
target_link_libraries(bioseqdb-pg PRIVATE ${PostgreSQL_LIBRARIES})
target_link_libraries(bioseqdb-pg PRIVATE ${BZIP2_LIBRARIES})
target_link_libraries(bioseqdb-pg PRIVATE ${HTS_LIBRARIES} ${BWA_LIBRARIES})
target_include_directories(bioseqdb-import PRIVATE ${PostgreSQL_INCLUDE_DIRS})
target_link_libraries(bioseqdb-import PRIVATE ${PostgreSQL_LIBRARIES})

install(TARGETS bioseqdb-pg DESTINATION ${PG_CONFIG_PKGLIBDIR})
install(FILES bioseqdb-pg/bioseqdb.control DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
install(FILES bioseqdb-pg/bioseqdb--0.0.0.sql DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
