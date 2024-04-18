VER=v2.2

build:
	docker build -t dbest/canu:${VER} .

build_no_cache:
	docker build --no-cache -t dbest/canu:${VER} .

push:
	docker push dbest/canu:${VER}

clean:
	rm *~

