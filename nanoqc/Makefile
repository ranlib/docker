#VER=latest
VER=v0.9.4

build:
	docker build -t dbest/nanoqc:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/nanoqc:$(VER) .

push:
	docker push dbest/nanoqc:$(VER)
