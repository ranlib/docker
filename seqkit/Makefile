VER=v2.10.0

build:
	docker build -t dbest/seqkit:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/seqkit:$(VER) .

push:
	docker push dbest/seqkit:$(VER)
