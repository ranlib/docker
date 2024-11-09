VER=v1.21

build:
	docker build -t dbest/samtools:$(VER) . 

build_no_cache:
	docker build --no-cache -t dbest/samtools:$(VER) . 

push:
	docker push dbest/samtools:$(VER)

clean:
	rm *~
