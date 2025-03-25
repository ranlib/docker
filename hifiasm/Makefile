VER=0.25.0
build:
	docker build -t dbest/hifiasm:$(VER) .

push: build
	docker push dbest/hifiasm:$(VER)

clean:
	rm *~

