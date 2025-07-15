VER=v1.22.1

build:
	docker build -t dbest/samtools:$(VER) . 

build_no_cache:
	docker build --no-cache -t dbest/samtools:$(VER) . 

push:
	docker push dbest/samtools:$(VER)

clean:
	rm *~
