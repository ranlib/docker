VER=v0.24.0

build:
	docker build -t dbest/fastp:$(VER) -f Dockerfile .

build_no_cache:
	docker build --no-cache -t dbest/fastp:$(VER) -f Dockerfile .

push:
	docker push  dbest/fastp:$(VER)

clean:
	rm *~
