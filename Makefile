OKTOBERFEST_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
# include $(OKTOBERFEST_DIR)MakefileShared # sourcing this is not needed because it is already sourced in MakefileCe
include $(OKTOBERFEST_DIR)MakefileShared

dependencies:
	git config --global credential.helper cache

registry:
	docker login gitlab.lrz.de:5005
	docker build -t gitlab.lrz.de:5005/proteomics/prosit_tools/oktoberfest .
	docker push gitlab.lrz.de:5005/proteomics/prosit_tools/oktoberfest

jump: 
	$(DOCKER_CMD) \
		$(IMAGE) bash

# --no-cache
build: dependencies
	docker build -f Dockerfile -t $(IMAGE) . || (exit 1)


run_oktoberfest: rm_err_file
	$(DOCKER_CMD) \
		$(IMAGE) python3 -u -m oktoberfest /root/data || (echo "2" > $(DATA)err.out; exit 2)

compress: run_oktoberfest
	zip -j -r -9 "$(DATA)/results.zip" "$(DATA)/percolator/"  --exclude '*_prosit.tab' '*_andromeda.tab' || (echo "3" > $(DATA)err.out; exit 3)

all: compress


run_local: 
	python3 -u -m oktoberfest "$(DATA)"

clean_data_folder: 
	rm -rf "$(DATA)/{proc,msms,percolator,mzML,msms.prosit}"
