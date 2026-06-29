# meson.build of the external library "blake3"
# Paul Boersma, 27 June 2026

blake3_objects = '''
	blake3.cpp blake3_dispatch.cpp blake3_portable.cpp'''.split()
	

sources = whisper_and_ggml_objects + ggml_cpu_objects + model_objects

libblake3 = static_library (
	'blake3',
	sources : blake3_objects,
	include_directories : [melder_inc]
)

libblake3_dep = declare_dependency (
	include_directories : include_directories('.'), link_with: libblake3
)
