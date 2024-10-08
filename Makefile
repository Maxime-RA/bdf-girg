LIB_BUILD_DIR=girgs_cpplib/build
LIB_TARGET_DIR=girgs_cpplib/install
LIB_GIRGS=girg_sampling/libgirgs.so.1
LIB_GIRGS_WRAP=girg_sampling/_libgirgs_wrapper.so

PYBIND11_INCLUDE=/home/maxime/.local/lib/python3.8/site-packages/pybind11/include
PYTHON_INCLUDE=$(shell poetry run python -c "from sysconfig import get_paths; print(get_paths()['include'])")

.PHONY: all clean build submodules test

all: build

clean:
	rm -rf dist/ build/ $(LIB_BUILD_DIR) $(LIB_TARGET_DIR)
	rm -rf girg_sampling/*.so girg_sampling/*.so.*
	rm -rf girg_sampling.egg-info/ __pycache__/

$(LIB_GIRGS) $(LIB_HYPERGIRGS):
	cmake girgs_cpplib -B $(LIB_BUILD_DIR)
	cmake --build $(LIB_BUILD_DIR)
	cmake --install $(LIB_BUILD_DIR) --prefix $(LIB_TARGET_DIR)
	# NOTE: renaming to prevent conflicts with other installs of lib(hyper)girgs
	cp $(LIB_TARGET_DIR)/lib/libgirgs.so $(LIB_GIRGS)

$(LIB_GIRGS_WRAP): girg_sampling/_libgirgs_wrapper.cpp $(LIB_GIRGS)
	$(CC) -fPIC -Wall -shared -g -o $@ $< \
        -I$(PYBIND11_INCLUDE) -I$(LIB_TARGET_DIR)/include/ -I$(PYTHON_INCLUDE) \
        -lstdc++ $(LIB_GIRGS) -Wl,-rpath=.

build: submodules $(LIB_GIRGS_WRAP) $(LIB_HYPERGIRGS_WRAP)

submodules:
	git submodule init
	git submodule update

test:
	poetry install
	poetry run pytest .
