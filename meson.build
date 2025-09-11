project(
  'dualization',
  'cpp',
  default_options : [
    'warning_level=2',
	'c_std=c17',
	'cpp_std=c++17',
#	'b_sanitize=address',
	'buildtype=release',
	'default_library=shared'
  ],
  version : '1.2.0'
)
major_version = '1'

main_source_file = 'src/experiment.cpp'
gen_source_file = 'src/generate_hypergraph.cpp'
project_source_files = [
  'src/bipartite_hypergraph.cpp',
  'src/bipartite_02.cpp'
]

project_include_dir = include_directories('include')

cc = meson.get_compiler('cpp')

lgmp_dep = dependency('libgmp', required : false)
if not lgmp_dep.found()
  lgmp_dep = cc.find_library('gmp', required : true)
endif

project_dependencies = [
  dependency('openmp'),
  lgmp_dep
]

build_args = [
#  '-fsanitize=address',
#  '-fno-omit-frame-pointer',
  '-flto',
  '-march=native'
]

linker_args = [
#  '-fsanitize=address',
#  '-static-libasan',
  '-flto',
  '-lgmp',
  '-lgmpxx'
]


# ===================================================================

# ======
# Target
# ======

build_args += [
  '-DPROJECT_NAME=' + meson.project_name(),
  '-DPROJECT_VERSION=' + meson.project_version(),
]

libdual = library(
  'libdual',
  project_source_files,
  version : meson.project_version(),
  soversion : major_version,
  include_directories: project_include_dir,
  dependencies: project_dependencies,
  c_args : build_args,
  link_args: linker_args)

dualization_target = executable(
  meson.project_name(),
  main_source_file,
  link_with : libdual,
  include_directories: project_include_dir,
  dependencies: project_dependencies,
  install : true,
  c_args : build_args,
  link_args: linker_args
)

generator_target = executable(
  'hypergraph_generator',
  gen_source_file,
  link_with : libdual,
  include_directories: project_include_dir,
  dependencies: project_dependencies,
  install : true,
  c_args : build_args,
  link_args: linker_args
)

#test('basic', project_target)
