venv := justfile_directory() / ".venv"
python := venv / "bin/python"
nanobind_dir := venv / "lib/python3.14/site-packages/nanobind/cmake"

# Configure CMake for compile_commands.json for clangd
configure:
    cmake -S . -B build \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
        -DPython_EXECUTABLE={{python}} \
        -Dnanobind_DIR={{nanobind_dir}}
