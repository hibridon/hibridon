SHELL:=/bin/bash 
TEMP_PATH:=$(shell echo ~/work/hibridon/tmp)
HIBRIDON_ROOT_PATH:=$(shell pwd)

CONFIG:
	./ci/make-config.bash $(HIBRIDON_ROOT_PATH) $(TEMP_PATH)

.PHONY: build
build: CONFIG
	./ci/build.bash $(HIBRIDON_ROOT_PATH) $(TEMP_PATH)

.PHONY: test
test: build
	./ci/test-hibridon.bash $(HIBRIDON_ROOT_PATH) $(TEMP_PATH)




