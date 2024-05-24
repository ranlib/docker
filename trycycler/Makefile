VER=v0.5.5

build:
	docker build -t dbest/trycycler:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/trycycler:$(VER) .

push:
	docker push dbest/trycycler:$(VER)
