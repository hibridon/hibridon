SHELL:=/bin/bash 
TEMP_PATH:=$(shell echo ~/work/hibridon/tmp)
HIBRIDON_ROOT_PATH:=$(shell pwd)
CONFIG_ID:=fr.univ-rennes1.ipr.physix.gfortran

CONFIG_FILE=CONFIG
LIB_SRC_FILES=$(shell find ./src -name "*.f")
LIB_BUILD_LOG_FILE=obj.log

.PHONY: all
all: test

.PHONY: config
config: $(CONFIG_FILE)

$(CONFIG_FILE):
	./ci/make-config.bash "$(HIBRIDON_ROOT_PATH)" "$(TEMP_PATH)" "$(CONFIG_ID)"

.PHONY: build
build: $(LIB_BUILD_LOG_FILE)

$(LIB_BUILD_LOG_FILE): CONFIG $(LIB_SRC_FILES)
	./ci/build.bash "$(HIBRIDON_ROOT_PATH)" "$(TEMP_PATH)" "$(CONFIG_ID)"

.PHONY: test
test: build
	./ci/test-hibridon.bash "$(HIBRIDON_ROOT_PATH)" "$(TEMP_PATH)" "$(CONFIG_ID)"




