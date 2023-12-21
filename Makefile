OKTOBERFEST_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
# include $(OKTOBERFEST_DIR)MakefileShared # sourcing this is not needed because it is already sourced in MakefileCe
include $(OKTOBERFEST_DIR)MakefileShared

dependencies:
	git config --global credential.helper cache

registry:
	docker login gitlab.lrz.de:5005
	docker build -t gitlab.lrz.de:5005/proteomics/github/oktoberfest .
	docker push gitlab.lrz.de:5005/proteomics/github/oktoberfest

jump:
	$(DOCKER_CMD) \
		$(IMAGE) bash

# --no-cache
build: dependencies
	git describe --long --dirty --always > hash.file
	docker build -f Dockerfile -t $(IMAGE) . || (exit 1)

bootstrap: DATA=/root/data
bootstrap:
	bash -c "cp /root/Makefile* $(LOCAL_DIR)"

run_oktoberfest: rm_err_file
	$(DOCKER_CMD) \
		$(IMAGE) python3 -u -m oktoberfest --config_path $(LOCAL_DIR)/config.json || (echo "2" > $(DATA)err.out; exit 2)

compress: run_oktoberfest
	zip -j -r -9 "$(DATA)/results.zip" "$(DATA)/results/" || (echo "3" > $(DATA)err.out; exit 3)

all: compress


run_local:
	python3 -u -m oktoberfest --config_path $(DATA)/config.json

clean_data_folder:
	bash -c "rm -rf $(DATA)/{proc,msms,results,mzML,msms.prosit,err.out,results.zip}"
