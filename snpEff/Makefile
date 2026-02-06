VER=v5.4a

build:
	docker build -t dbest/snpeff:$(VER) . 

build_no_cache:
	docker build --no-cache -t dbest/snpeff:$(VER) . 

push:
	docker push dbest/snpeff:$(VER)

clean:
	rm *~
