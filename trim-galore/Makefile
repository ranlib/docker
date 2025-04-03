VER=v0.6.10

build:
	docker build -t dbest/trim_galore:$(VER) -f Dockerfile .

build_no_cache:
	docker build --no-cache -t dbest/trim_galore:$(VER) -f Dockerfile .

push:
	docker push  dbest/trim_galore:$(VER)

clean:
	rm *~
