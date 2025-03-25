VER=v2.9.3

build:
	docker build -t dbest/flye:${VER} .

build_no_cache:
	docker build --no-cache -t dbest/flye:${VER} .

push:
	docker push dbest/flye:${VER}

clean:
	rm *~

