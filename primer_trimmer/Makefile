VER=v0.1

build:
	docker build -t dbest/primer_trimmer:$(VER) . 

build_no_cache:
	docker build --no-cache -t dbest/primer_trimmer:$(VER) . 

push:
	docker push dbest/primer_trimmer:$(VER)

clean:
	rm *~
