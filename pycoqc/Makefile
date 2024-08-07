VER=latest

build:
	docker build -t dbest/pycoqc:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/pycoqc:$(VER) .

push:
	docker push dbest/pycoqc:$(VER)
