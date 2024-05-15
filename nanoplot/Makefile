VER=v1.42.0

build:
	docker build -t dbest/nanoplot:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/nanoplot:$(VER) .

push:
	docker push dbest/nanoplot:$(VER)
